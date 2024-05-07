# `bean create-screen`: Create ReporterScreen object from flat files
```bash
bean create-screen gRNA_library.csv sample_list.csv gRNA_counts_table.csv
```
## Input
  * gRNA_library.csv (`variant` screen [example](https://github.com/pinellolab/crispr-bean/blob/main/tests/data/var_mini_guides.csv), `tiling` screen [example](https://github.com/pinellolab/crispr-bean/blob/main/tests/data/test_guide_info_tiling_chrom.csv))
  * sample_list.csv (`sorting` screen [example](https://github.com/pinellolab/crispr-bean/blob/main/tests/data/sample_list_survival.csv), `variant` screen [example](https://github.com/pinellolab/crispr-bean/blob/main/tests/data/var_mini_samples.csv))
  * gRNA_counts_table.csv: Table with gRNA ID in the first column and sample IDs as the column names (first row)
`gRNA_library.csv` and `sample_list.csv` should be formatted as :ref:`input`. ([example](https://github.com/pinellolab/crispr-bean/blob/main/tests/data/var_mini_counts.csv))