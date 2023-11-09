# <img src="imgs/bean_title.svg" alt="crispr-bean" width="300"/>

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/crispr-bean)](https://pypi.org/project/crispr-bean/)
[![PyPI version](https://img.shields.io/pypi/v/crispr-bean)](https://pypi.org/project/crispr-bean/)
[![Code style](https://img.shields.io/badge/code%20style-black-black)](https://github.com/psf/black)

**CRISPR** **B**ayesian **E**stimation of variant effect (from **B**ase **E**diting reporter screens) with guide **A**ctivity **N**ormalization  
This is an analysis toolkit for the pooled CRISPR reporter or sensor data. The reporter technique transfects cells with plasmid with not only sgRNA but with the **target sequence surrogate** which we call **reporter** or **sensor**.  


<img src="imgs/reporter.jpg" alt="Reporter construct" width="700"/>

## Overview
`crispr-bean` supports end-to-end analysis of pooled sorting screens, with or without reporter.  

<img src="imgs/dag_bean.png" alt="dag_bean.svg" width="700"/>  

1. [`bean-count-sample`](#bean-count-samples-count-reporter-screen-data): Base-editing-aware **mapping** of guide, optionally with reporter from `.fastq` files.  
2. [`bean-qc`](#bean-qc-qc-of-reporter-screen-data): Quality control report and filtering out / masking of aberrant sample and guides  
3. [`bean-filter`](#bean-filter-filtering-and-optionally-translating-alleles): Filter reporter alleles; essential for `tiling` mode that allows for all alleles generated from gRNA.
4. [`bean-run`](#bean-run-quantify-variant-effects): Quantify targeted variants' effect sizes from screen data.  

### Data structure
BEAN stores mapped gRNA and allele counts in `ReporterScreen` object which is compatible with [AnnData](https://anndata.readthedocs.io/en/latest/index.html). See [Data Structure](#data-structure) section for more information.

### Examples
We provide example scripts in `tests/`. Running `pytest --sparse-ordering` generates example input/output files from running 1 and 2-4 sequentially.

### Pipeline run options by library design
The `bean-filter` and `bean-run` steps depend on the type of gRNA library design, where BEAN supports two modes of running.
1. `variant` library: Several gRNAs tile each of the targeted variants  
  Ex)  
  <img src="imgs/variant.png" alt="variant library design" width="700"/>  

2. `tiling` library: gRNA densely tiles a long region (e.g. gene(s), exon(s), coding sequence(s))  
  Ex)  
  <img src="imgs/tiling.png" alt="tiling library design" width="450"/>  

<br/><br/>

## Installation 
### Full installation
First install [pytorch](https://pytorch.org/) >=12.1,<2.
Then download from PyPI:
```
pip install crispr-bean[model]
```

### Mapping and data wrangling, without variant effect quantification
```
pip install crispr-bean
```
This wouldn't have variant effect size quantification (`bean-run`) functionality.  

<br/><br/>

## `bean-count-samples`: Count (reporter) screen data  
`bean-count-samples` (or `bean-count` for a single sample) maps guide into guide counts, **allowing for base transition in spacer sequence**. When the matched reporter information is provided, it can count the **target site edits** and **alleles produced by each guide**. Mapping is efficiently done based on [CRISPResso2](https://github.com/pinellolab/CRISPResso2) modified for base-edit-aware mapping.



```python
bean-count-samples \
  --input sample_list.csv   `# sample with lines 'R1_filepath,R2_filepath,sample_name\n'` \
  -b A                      `# base that is being edited (A/G)` \
  -f sgRNA_info_table.csv   `# sgRNA information` \
  -o .                      `# output directory` \
  -r                        `# read edit/allele information from reporter` \
  -t 12                     `# number of threads` \
  --name my_sorting_screen  `# name of this sample run`
```

### Input file format
#### 1. gRNA_library.csv
File should contain following columns. 
* `name`: gRNA ID column
* `sequence`: gRNA sequence
* `barcode`: R2 barcode to help match reporter to gRNA
* In order to use accessibility in the [variant effect quantification](#bean-run-quantify-variant-effects), provide accessibility information in one of two options. (For non-targeting guides, provide NA values (empty cell).)   
  * Option 1: `chrom` & `genomic_pos`: Chromosome (ex. `chr19`) and genomic position of guide sequence. You will have to provide the path to the bigwig file with matching reference version in `bean-run`. 
  * Option 2: `accessibility_signal`: ATAC-seq signal value of the target loci of each guide.  
* For variant library (gRNAs are designed to target specific variants and ignores bystander edits)
  * `target`: This column denotes which target variant/element of each gRNA. This is not used in `bean-count[-samples]` but required to run `bean-run` in later steps.
  * `target_group [Optional]`: If negative/positive control gRNA will be considered in `bean-qc` and/or `bean-run`, specify as "NegCtrl"/"PosCtrl" in this column. 
  * `target_pos [Optional]`: If `--match_target_pos` flag is used, input file needs `target_pos` which specifies 0-based relative position of targeted base within Reporter sequence.
* For tiling library (gRNAs tile coding / noncoding sequences)
  * `strand`: Specifies gRNA strand information relative to the reference genome.
  * `chrom`: Chromosome of gRNA targeted locus.
  * `start_pos`: gRNA starting position in the genome. Required when you provide `strand` column. Should specify the smaller coordinate value among start and end position regardless of gRNA strandedness.

Also see examples for [variant library](tests/data/test_guide_info.csv) and [tiling library](tests/data/test_guide_info_tiling.csv).

#### 2. sample_list.csv
File should contain following columns with header.
* `R1_filepath`: Path to read 1 `.fastq[.gz]` file
* `R2_filepath`: Path to read 1 `.fastq[.gz]` file
* `sample_id`: ID of sequencing sample
* `rep [Optional]`: Replicate # of this sample
* `bin [Optional]`: Name of the sorting bin
* `upper_quantile [Optional]`: FACS sorting upper quantile
* `lower_quantile [Optional]`: FACS sorting lower quantile  

Optional columns are not required but can be provided for compatibility with `bean-qc` and `bean-run`. See [example](tests/data/sample_list.csv).
  
### Output file format
`count` or `count-samples` produces `.h5ad` and `.xlsx` file with guide and per-guide allele counts.  
* `.h5ad`: This output file follows annotated matrix format compatible with `AnnData` and is based on `Screen` object in [purturb_tools](https://github.com/pinellolab/perturb-tools). See [Data Structure](#data-structure) section for more information.  
* `.xlsx`: This output file contains `.guides`, `.samples`, `.X[_bcmatch,_edits]`. (`allele_tables` are often too large to write into an Excel!)  

### Parameters
* `-b`, `--edited-base` (`str`, default: `None`): For base editors, the base that should be ignored when matching the gRNA sequence 
* `-f`, `--sgRNA-filename` (`str`, default: `None`): sgRNA description file. (See above)
* `--guide-start-seq`: Guide starts after this sequence in R1 (default: '')
* `--guide-end-seq`: Guide ends after this sequence in R1 (default: '')
* `-r`, `--count-reporter` (default: `False`): Count edited alleles in reporters.
* `-q`, `--min-average-read-quality` (default: `30`): Minimum average quality score (phred33) to keep a read
* `-s`, `--min-single-bp-quality` (default: `0`): Minimum single bp score (phred33) to keep a read (default: 0)
* `-n`, `--name`: Name of the output file will be `bean_count_{name}.h5ad`.
* `-o`, `--output-folder`: Output folder
* `-l`, `--reporter-length` (default: `32`): Length of the reporter sequence.
* `--keep-intermediate` (default: `False`): Keep all intermediate files generated from filtering.
* `--qstart-R1` (default: `0`): Start position of the read when filtering for quality score of the read 1 
* `--qend-R1` (default: `47`): End position of the read when filtering for quality score of the read 1 
* `--qstart-R2` (default: `0`): Start position of the read when filtering for quality score of the read 2 
* `--qend-R2` (default: `36`): End position of the read when filtering for quality score of the read 2 
* `--gstart-reporter` (default: `6`): Start position of the guide sequence in the reporter
* `--match-target-pos` (default: `False`): Count the edit in the exact target position.
* `--target-pos-col` (default: `target_pos`): Column name specifying the relative target position within *reporter* sequence. 
* `--guide-bc` (default: `True`): Construct has guide barcode 
* `--guide-bc-len` (default: `4`): Guide barcode sequence length at the beginning of the R2
* `--offset` (default: `False`): Guide file has `offest` column that will be added to the relative position of reporters. 
* `--align-fasta` (default: ''): gRNA is aligned to this sequence to infer the offset. Can be used when the exact offset is not provided.
* `--string-allele` (default: `False`): Store allele as quality filtered string instead of Allele object
* `-g`, `--count-guide-edits` (default: `False`): count the self editing of guides (default: False)
* `-m`, `--count-guide-reporter-alleles`: count the matched allele of guide and reporter edit
* `--tiling`: Specify that the guide library is tiling library 
* `-i`, `--input`: List of fastq and sample ids. See above for the format.
* `-t`, `--threads` (default: `10`): Number of threads to use
* `--guide-start-seqs-file` (default: `None`): CSV file path with per-sample `guide_start_seq` to be used, if provided. Formatted as `sample_id, guide_start_seq` 
* `--guide-end-seqs-file` (default: `None`): CSV file path with per-sample `guide_end_seq` to be used, if provided. Formatted as `sample_id,guide_end_seq` 
* `--rerun` (default: `False`): Recount each sample. If `False`, existing count for each sample is taken.




<br/><br/>

## `bean-qc`: QC of reporter screen data
```bash
bean-qc \
  my_sorting_screen.h5ad    `# Input ReporterScreen .h5ad file path` \
  -o my_sorting_screen_masked.h5ad   `# Output ReporterScreen .h5ad file path` \
  -r qc_report_my_sorting_screen   `# Prefix for QC report` 
```

`bean-qc` supports following quality control and masks samples with low quality. Specifically:
* Plots guide coverage and the uniformity of coverage
* Guide count correlation between samples
* Log fold change correlation when positive controls are provided
* Plots editing rate distribution
* Identify samples with low guide coverage/guide count correlation/editing rate and mask the sample in `bdata.samples.mask`
* Identify outlier guides to filter out

#### Output
Above command produces 
* `my_sorting_screen_masked.h5ad` without problematic replicate and guides and with sample masks, and  
* `qc_report_my_sorting_screen.[html,ipynb]` as QC report.  


#### Additional Parameters
* `--tiling` (default: `None`): If set as `True` or `False`, it sets the screen object to be tiling (`True`) or variant (`False`)-targeting screen when calculating editing rate. 
* `--replicate-label` (default: `"rep"`): Label of column in `bdata.samples` that describes replicate ID.
* `--condition-label` (default: `"bin"`)": Label of column in `bdata.samples` that describes experimental condition. (sorting bin, time, etc.).
* `--target-pos-col` (default: `"target_pos"`): Target position column in `bdata.guides` specifying target edit position in reporter.
* `--rel-pos-is-reporter` (default: `False`): Specifies whether `edit_start_pos` and `edit_end_pos` are relative to reporter position. If `False`, those are relative to spacer position.
* `--edit-start-pos` (default: `2`): Edit start position to quantify editing rate on, 0-based inclusive.
* `--edit-end-pos` (default: `7`): Edit end position to quantify editing rate on, 0-based exclusive.
* `--count-correlation-thres` (default: `0.8`): Threshold of guide count correlation to mask out.
* `--edit-rate-thres` (default: `0.1`): Median editing rate threshold per sample to mask out.
* `--posctrl-col` (default: `group`): Column name in .h5ad.guides DataFrame that specifies guide category.
* `--posctrl-val` (default: `PosCtrl`): Value in .h5ad.guides[`posctrl_col`] that specifies guide will be used as the positive control in calculating log fold change.
* `--lfc-thres` (default: `0.1`): Positive guides' correlation threshold to filter out.
* `--lfc-conds` (default: `"top,bot"`): Values in of column in `ReporterScreen.samples[condition_label]` for LFC will be calculated between, delimited by comma
* `--recalculate-edits` (default: `False`): Even when `ReporterScreen.layers['edit_count']` exists, recalculate the edit counts from `ReporterScreen.uns['allele_count']`."

<br/><br/>


## `bean-filter`: Filtering (and optionally translating) alleles
As `tiling` mode of `bean-run` accounts for any robustly observed alleles, `bean-filter` filters for such alleles.
```bash
bean-filter my_sorting_screen_masked.h5ad \
-o my_sorting_screen_filtered.h5ad  `# Output file path` \
```

### Output
Above command produces 
* `my_sorting_screen_filtered.h5ad` with filtered alleles stored in `.uns`,   
* `my_sorting_screen_filtered.filtered_allele_stats.pdf`, and `my_sorting_screen_filtered.filter_log.txt` that report allele count stats in each filtering step.  

You may want to adjust the flitering parameters to obtain optimal balance between # guides per variant & # variants that are scored. See example outputs of filtering step [here](docs/example_filtering_outputs/).


### Translating alleles
If you want to obtain **amino acid level variant** for coding sequence tiling screens, provide coding sequence positions which variants occuring within the coding sequence will be translated. *This is optional, but **highly recommended** to increase per-(coding)variant support.*  

<img src="imgs/translation.png" alt="Allele translation" width="500"/>  
  

```bash
bean-filter my_sorting_screen.h5ad \
-o my_sorting_screen_masked.h5ad \
--translate   `# Translate coding variants` \
[ --translate-gene-name GENE_SYMBOL OR
  --translate-genes-list path_to_gene_names_file.txt OR
  --translate-fasta gene_exon.fa, OR
  --translate-fastas-csv gene_exon_fas.csv]
```
* When library covers a single gene, do either of the following:
  1. Feed `--translate-gene-name GENE_SYMBOL` if your `genomic_pos` column of `sgRNA_info_tbl` is compatible with [MANE transcript](https://useast.ensembl.org/info/genome/genebuild/mane.html)'s reference genome. (Per 10/23/2023, GRCh38). This will automatically load the exon positions based on MANE transcript annotation.
  2. To use your custom coding sequence and exon positions, feed `--translate-fasta gene_exon.fa` argument where `gene_exon.fa` is the FASTA file with entries of exons. [See full details here](docs/exon_fa_format.md).
* When library covers multiple genes, do either of the following:  
  1. Feed `--translate-genes-list path_to_gene_names_file.txt` where `path_to_gene_names_file.txt` is file with one gene symbol per line.
  2. Feed `--translate-fastas-csv gene_exon_fas.csv` where `gene_exon_fas.csv` is the csv file with lines `gene_id,gene_exon_fasta_path` without header. Each FASTA file in `gene_exon_fasta_path` is formatted [as the single-gene FASTA file](docs/exon_fa_format.md).
* Translation will keep the variants outside the coding sequence as nucleotide-level variants, while aggregating variants leading to the same coding sequence variants.

### Full list of parameters
* `-o`, `--output-prefix` (default: `None`): Output prefix for log and ReporterScreen file with allele assignment.
* `-p`, `--plasmid-path` (default: `None`): Plasmid `ReporterScreen` object path. 
  * If provided, alleles are filtered based on if a nucleotide edit is more significantly enriched in sample compared to the plasmid data. 
  * Negative control data where no edit is expected can be fed in instead of plasmid library.
* `-w`, `--filter-window` (default: `False`): "Only consider edit within window provided by (`edit_start_pos`, `edit_end_pos`). If this flag is not provided, `--edit-start-pos` and `--edit-end-pos` flags are ignored.
* `-s`, `--edit-start-pos` (default: `2`): 0-based start posiiton (inclusive) of edit relative to the start of guide spacer.
* `-e`, `--edit-end-pos` (default: `7`): 0-based end position (exclusive) of edit relative to the start of guide spacer.
* `-b`, `--filter-target-basechange` (default: `False`): Filter for target (intended) base change of edits (stored in `bdata.uns['target_base_change']`)
* `-j`, `--jaccard-threshold` (default: `0.3`): Jaccard Index threshold when the alleles are mapped to the most similar alleles. In each filtering step, allele counts of filtered out alleles will be mapped to the most similar allele only if they have Jaccard Index of shared edit higher than this threshold.
* `--translate` (default: `False`): Translate nucleotide-level variants prior to allele proportion filtering.
* `-f`, `--translate-fasta` (defulat: `None`): fasta file path with exon positions. If not provided and `--translate` flag is provided, LDLR hg19 coordinates will be used.
* `-fs`, `--translate-fastas-csv` (defulat: `None`): .csv with two columns with gene IDs and FASTA file path corresponding to each gene.
* `-g`, `--translate-gene-name` (default: `None`): Gene symbol for translation
* `-gs`, `--translate-genes-list` (default: `None`): Path to the text file with gene symbols in each line
* `-ap`, `--filter-allele-proportion` (default: `0.05`): If provided, only the alleles that exceed `filter_allele_proportion` in `filter-sample-proportion` will be retained.
* `-ac`, `--filter-allele-count` (default: `5`): If provided, alleles that exceed `filter_allele_proportion` AND `filter_allele_count` in `filter-sample-proportion` will be retained.
* `sp`, `--filter-sample-proportion` (default: `0.2`): "If `filter_allele_proportion` is provided, alleles that exceed `filter_allele_proportion` in `filter-sample-proportion` will be retained.


<br/><br/>

## `bean-run`: Quantify variant effects
BEAN uses Bayesian network to incorporate gRNA editing outcome to provide posterior estimate of variant phenotype. The Bayesian network reflects data generation process. Briefly,  
1. Cellular phenotype (either for cells are sorted upon for sorting screen, or log(proliferation rate)) is modeled as the Gaussian mixture distribution of wild-type phenotype and variant phenotype.
2. The weight of the mixture components are inferred from the reporter editing outcome and the chromatin accessibility of the loci.
3. Cells with each gRNA, formulated as the mixture distribution, is sorted by the phenotypic quantile to produce the gRNA counts.

For the full detail, see the method section of the [BEAN manuscript](https://www.medrxiv.org/content/10.1101/2023.09.08.23295253v1).

<img src="imgs/bean.gif" alt="model" width="700"/>   
  
<br></br>
```bash
bean-run sorting[survival] variant[tiling] my_sorting_screen_filtered.h5ad \
[--uniform-edit, --scale-by-acc [--acc-bw-path accessibility_signal.bw, --acc-col accessibility]] \
-o output_prefix/ \
--fit-negctrl
```

### Input
`my_sorting_screen_filtered.h5ad` can be produced by one of the following:  
1. [`bean-count-samples`]((#bean-count-samples-count-reporter-screen-data)) when you have raw `.fastq` file
2. (Limited to `bean-run variant` mode)  `bean-create-screen` when you have flat `.csv` tables of gRNA metadata table, sample metadata table, gRNA counts table (# guides x # samples), and optionally # edits table.   
    ```bash
    bean-create-screen gRNA_info_table.csv sample_info_table.csv gRNA_counts_table.csv \
    [--edits edit_counts_table.csv -o output.h5ad] 
    ```  
    * `gRNA_info_table.csv` should have following columns.
      * `name`: gRNA ID column
      * `target`: This column denotes which target variant/element of each gRNA.
      * `target_group [Optional]`: If negative control gRNA will be used, specify as "NegCtrl" in this column. 
    * `sample_info_table.csv` should have following columns.
      * `sample_id`: ID of sequencing sample
      * `rep`: Replicate # of this sample
      * `bin`: Name of the sorting bin
      * `upper_quantile`: FACS sorting upper quantile
      * `lower_quantile`: FACS sorting lower quantile  
    * `gRNA_counts_table.csv` should be formatted as follows.
      * Columns include one of `sample_id` columns in `sample_info_table.csv` file.
      * 1st row (row index) follows `name` (gRNA ID) in `gRNA_info_table.csv` file.
3. You can manually create the `AnnData` object with more annotations including allele counts: see [API tutorial](#using-bean-as-python-module) for full detail.


### Output
<img src="imgs/model_output.png" alt="model" width="700"/>

Above command produces
* `output_prefix/bean_element_result.[model_type].csv` with following columns:
  * `mu` (Effect size): Mean of variant phenotype, given the wild type has standard normal phenotype distribution of `mu = 0, sd = 1`.
  * `mu_sd`: Mean of variant phenotype `mu` is modeled as normal distribution. The column shows fitted standard deviation of `mu` that quantify the uncertainty of the variant effect.
  * `mu_z`: z-score of `mu`
  * `sd`: Standard deviation of variant phenotype, given the wild type has standard normal phenotype distribution of `mu = 0, sd = 1`.
  * `CI[0.025`, `0.975]`: Credible interval of `mu`
  * When negative control is provided, above columns with `_adj` suffix are provided, which are the corresponding values adjusted for negative control.  
* `output_prefix/bean_sgRNA_result.[model_type].csv`: 
  * `edit_rate`: Estimated editing rate at the target loci.

### Optional parameters
* Run options
  * `-o`, `--outdir`: Directory to save the run result.
  * `--result-suffix` (default: ""): Suffix of the output files.
  * `--uniform-edit` (default: `False`): Ignores variable editing
  * `--scale-by-acc` (default: `False`): Scale guide editing efficiency by the target loci accessibility. Needs one of `acc_bw_path` or `acc_col` that provides per-guide accessibility information.
    * `--acc-bw-path` : Accessibility .bigWig file to be used to assign accessibility of guides. To read accessibility information from the .bigWig file, input `.h5ad` file needs `chr` and `genomic_pos` column in `.guides` DataFrame. 
  * `--fit-negctrl` (default: `False`): Fit the shared negative control distribution to normalize the fitted parameters.
    * `--negctrl-col` (default: `target_group`): Column in bdata.obs specifying if a guide is negative control. If the `bdata.guides[negctrl_col].tolower() == negctrl_col_value`, it is treated as negative control guide.
    * `--negctrl-col-value` (default: `negctrl`): Column value in bdata.guides specifying if a guide is negative control. If the `bdata.guides[negctrl_col].tolower() == negctrl_col_value`, it is treated as negative control guide.
  * `--repguide-mask` (default: `repguide_mask`): n_replicate x n_guide mask to mask the outlier guides. bdata.uns[repguide_mask] will be used. This is calculated with `bean-qc`. 
  * `--cuda` (default: `False`): Run on GPU
    * `--device`: Optionally use GPU if provided valid GPU device name (ex. cuda:0)
  * `--ignore-bcmatch` (default: `False`): If provided, even if the screen object has .X_bcmatch, ignore the count when fitting.
  * `--allele-df-key` (default: `allele_counts`): bdata.uns[allele_df_key] will be used as the allele count.
  * `--control-guide-tag` (default: `CONTROL`): (Relevant in bean-run *tiling* mode) If this string is in guide name, treat each guide separately not to mix the position. Used for non-targeting negative controls.
* Guide annotations (`bdata.guides` column keys)
  * `--acc-col`: Column name in bdata.guides that specify raw ATAC-seq signal.
  * `--target-column` (default: `target`): Column key in `bdata.guides` that describes the target element of each guide.
  * `--guide-activity-col`: Column in `bdata.guides` DataFrame showing the editing rate estimated via external tools.
* Sample annotations (`bdata.samples` column keys)
  * `--condition-column` (default: `bin`): Column key in `bdata.samples` that describes experimental condition.
  * `-uq`, `--sorting-bin-upper-quantile-column` (default: `upper_quantile`): Column name with upper quantile values of each sorting bin in bdata.samples
  * `-lq`, `--sorting-bin-lower-quantile-column` (default: `lower_quantile`): Column name with lower quantile values of each sorting bin in bdata.samples

<br></br>

## Data Structure
### ReporterScreen object
BEAN stores mapped gRNA and allele counts in `ReporterScreen` object which is compatible with [AnnData](https://anndata.readthedocs.io/en/latest/index.html).  

<img src="imgs/data_structure_v2.png" alt="ReporterScreen object structure" width="700" />

  * `.guides`: guide information provided in input (`gRNA_library.csv` in above example)
  * `.samples`: sample information provided in input (`sample_list.csv` in above example)
  * `.X`: Main guide count matrix, where row corresponds to each guide in `.guides` and columns correspond to samples in `.samples`.
Following attributes are included if matched reporter is provided and you chose to read edit/allele information from the reporter using `-r` option.
  * `.X_bcmatch [Optional]`: Contains information about number of barcode-matched reads. Information about R2 barcode should be specified as `barcode` column in your `gRNA_library.csv` file.
  * `.X_edits [Optional]`: If target position of each guide is specified as `target_pos` in input `gRNA_library.csv` file and `--match-target-position` option is provided, the result has the matrix with the number of target edit at the specified positions.
  * `.allele_tables [Optional]`: Dictionary with a single allele count table that counts per guide and allele combination, what is the count per sample. 

### Using BEAN as Python module
```
import bean as be
cdata = be.read_h5ad("bean_counts_sample.h5ad")
```
Python package `bean` supports multiple data wrangling functionalities for `ReporterScreen` objects. See the [**ReporterScreen API tutorial**](docs/ReporterScreen_api.ipynb) for more detail.

<br></br>

## Run time
* Installation takes 14.4 mins after pytorch installation with pytorch in Dell XPS 13 Ubuntu WSL.
* `bean-run` takes 4.6 mins with `--scale-by-acc` tag in Dell XPS 13 Ubuntu WSL for variant screen dataset with 3455 guides and 6 replicates with 4 sorting bins.
* Full pipeline takes 90.1s in GitHub Action for toy dataset of 2 replicates and 30 guides.

## Citation
If you have used BEAN, please cite:  
Ryu, J. et al. Joint genotypic and phenotypic outcome modeling improves base editing variant effect quantification. medRxiv (2023) doi:10.1101/2023.09.08.23295253
