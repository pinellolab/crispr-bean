# `bean run`: Quantify variant effects
BEAN uses Bayesian network to incorporate gRNA editing outcome to provide posterior estimate of variant phenotype. The Bayesian network reflects data generation process. Briefly,  
1. Cellular phenotype (either for cells are sorted upon for sorting screen, or log(proliferation rate)) is modeled as the Gaussian mixture distribution of wild-type phenotype and variant phenotype.
2. The weight of the mixture components are inferred from the reporter editing outcome and the chromatin accessibility of the loci.
3. Cells with each gRNA, formulated as the mixture distribution, is sorted by the phenotypic quantile to produce the gRNA counts.

For the full detail, see the method section of the [BEAN manuscript](https://www.medrxiv.org/content/10.1101/2023.09.08.23295253v1).

<img src="/crispr-bean/assets/bean.gif" alt="model" width="700"/>   
  
<br></br>

# Usage example
```bash
bean run sorting[survival] variant[tiling] my_sorting_screen_filtered.h5ad \
[--uniform-edit, --scale-by-acc [--acc-bw-path accessibility_signal.bw, --acc-col accessibility]] \
-o output_prefix/ \
--fit-negctrl
```
See full list of parameters [below](#full-parameters).


# Input
`my_sorting_screen_filtered.h5ad` can be produced by one of the following:  
1. [`bean count-samples`]((#bean-count-samples-count-reporter-screen-data)) when you have raw `.fastq` file
2. (Limited to `bean run variant` mode)  `bean create-screen` when you have flat `.csv` tables of gRNA metadata table, sample metadata table, gRNA counts table (# guides x # samples), and optionally # edits table.   
    ```bash
    bean create-screen gRNA_info_table.csv sample_info_table.csv gRNA_counts_table.csv \
    [--edits edit_counts_table.csv -o output.h5ad] 
    ```  
    * `gRNA_info_table.csv` should have following columns.
      * `name`: gRNA ID column
      * `target`: This column denotes which target variant/element of each gRNA.
      * `target_group [Optional]`: If negative control gRNA will be used, specify as "NegCtrl" in this column. 
    * `sample_info_table.csv` should have following columns.
      * `sample_id`: ID of sequencing sample
      * `replicate`: Replicate # of this sample
      * `bin`: Name of the sorting bin
      * `upper_quantile`: FACS sorting upper quantile
      * `lower_quantile`: FACS sorting lower quantile  
    * `gRNA_counts_table.csv` should be formatted as follows.
      * Columns include one of `sample_id` columns in `sample_info_table.csv` file.
      * 1st row (row index) follows `name` (gRNA ID) in `gRNA_info_table.csv` file.
3. You can manually create the `AnnData` object with more annotations including allele counts: see [API tutorial](#using-bean-as-python-module) for full detail.


# Output
<img src="/crispr-bean/assets/model_output.png" alt="model" width="700"/>

Above command produces
* `output_prefix/bean_element_result.[model_type].csv` with following columns:
  * Estimated variant effect sizes
    * `mu` (Effect size): Mean of variant phenotype, given the wild type has standard normal phenotype distribution of `mu = 0, sd = 1`.
    * `mu_sd`: Mean of variant phenotype `mu` is modeled as normal distribution. The column shows fitted standard deviation of `mu` that quantify the uncertainty of the variant effect.
    * `mu_z`: z-score of `mu`
    * `sd`: Standard deviation of variant phenotype, given the wild type has standard normal phenotype distribution of `mu = 0, sd = 1`.
    * `CI[0.025`, `0.975]`: Credible interval of `mu`
    * When negative control is provided, above columns with `_adj` suffix are provided, which are the corresponding values adjusted for negative control.  
  * Metrics on per-variant evidence provided in input (provided in `tiling` mode)
    * `effective_edit_rate`: Sum of per-variant editing rates over all alleles observed in the input. Allele-level editing rate is divided by the number of variants observed in the allele prior to summing up.
    * `n_guides`: # of guides covering the variant.
    * `n_coocc`: # of cooccurring variants with a given variant in any alleles observed in the input.
* `output_prefix/bean_sgRNA_result.[model_type].csv`: 
  * `edit_rate`: Estimated editing rate at the target loci.

See the example output [here](https://github.com/pinellolab/crispr-bean/tree/main/docs/example_run_output).