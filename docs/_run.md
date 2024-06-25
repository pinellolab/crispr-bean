# `bean run`: Quantify variant effects
BEAN uses Bayesian network to incorporate gRNA editing outcome to provide posterior estimate of variant phenotype. The Bayesian network reflects data generation process. Briefly,  
1. Cellular phenotype (either for cells are sorted upon for sorting screen, or log(proliferation rate)) is modeled as the Gaussian mixture distribution of wild-type phenotype and variant phenotype.
2. The weight of the mixture components are inferred from the reporter editing outcome and the chromatin accessibility of the loci.
3. Cells with each gRNA, formulated as the mixture distribution, is sorted by the phenotypic quantile to produce the gRNA counts.

For the full detail on modeling, see the [model description](https://pinellolab.github.io/crispr-bean/model.html).

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
## `bean_element_result.[model_type].csv`
- Variant ID / grouping
  - `edit`: Variant ID.
  - `group`: The grouping of the coding variants, assigned as one of nonsense/missense/synonymous.
  - `int_pos`: The integer position of the noncoding variants.
  - `chrom`: The chromosome of the variant.
  - `pos`: The position of the variant. If coding variant, starts with `A` and the position specified 1-based amino acid position. If noncodig variant, numeric genomic position.
  - `ref`: The reference base/amino acid of the variant.
  - `alt`: The alternative base/amino acid of the variant.
  - `coding`: A flag indicating if the element is coding variant or not.

- Per-variant summary of variant-producing guides (`tiling` mode)
  - `guide_target_group`: Aggregated `target_group` column in the input sgRNA_info.csv file. All unique values of the guides that produced (filtered) edited alleles that includes this variant is listed.
  - `effective_edit_rate`: The effective editing rate of the element. Calculated as `sum_over_guides(sum_over_alleles(per_guide_allele_editing_rate / # variants in the allele))`.
  - `editing_guides`: List of guides that edited the variant.
  - `per_guide_editing_rates`: The per-guide editing rates of the variant.
  - `n_guides`: The number of guides that edited the variant.
  - `n_coocc`: The number of unique co-occurring variants that appeared together in any alleles that contains the variant.
  
- Variant effect size: Use `mu_z_adj` whenever available, otherwise `mu_z_scaled`, otherwise `mu_z`.
  - `mu`: The mean value of the variant effect size.
  - `mu_sd`: The standard deviation of the mean value of the variant effect size.
  - `mu_z`: The z-score of the mean value of the variant effect size.
  - `sd`: The standard deviation of the phenotype induced by the variant.
  - `CI[0.025,0.975]`: The 95% credible interval of the mean value of the variant effect size. Corresponds to `mu_z_adj` when available, otherwise `mu_z_scaled`, otherwise `mu_z`. 
  - `[]_scaled`: Above values scaled by negative control variants.
  - `[]_adj`: Above values scaled by synonymous variants.
  
## `bean_sgRNA_result.[model_type].csv`
- `name`: sgRNA ID provided in the `name` column of the input.
- `edit_rate`: Effective editing rates
- `accessibility`: (Only if you have used `--scale-by-acc`) Accessibility signal that is used for scaling of the editing rate.
- `scaled_edit_rate`: (Only if you have used `--scale-by-acc`) Endogenous editing rate used for modeling, estimated by scaling reporter editing rate by accessibility signal
- `[cond1]_[cond2].median_lfc`: Raw LFC with pseudocount fed in with `--guide-lfc-pseudocount` argument (default 5).
- For `tiling` mode
  - `variants`: Variants generated by this gRNA
  - `variant_edit_rates`: Editing rate of this gRNA for each variant it creates.

See the full output file description and example output [here](https://github.com/pinellolab/crispr-bean/tree/main/docs/example_run_output).