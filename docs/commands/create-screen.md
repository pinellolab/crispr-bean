## `bean create-screen`: Create ReporterScreen object from flat files
```bash
bean create-screen gRNA_library.csv sample_list.csv gRNA_counts_table.csv
```
### Input
  * [gRNA_library.csv](#1-gRNA_librarycsv)
  * [sample_list.csv](#2-sample_listcsv)
  * gRNA_counts_table.csv: Table with gRNA ID in the first column and sample IDs as the column names (first row)

### Full Parameters
  * `-e`, `--edits` (default: `None`): Path to edit counts .csv table, with index at first column and column names at the first row.
  * `-o`, `--output-prefix` (default: `None`): Output file prefix (output will be saved as `output_prefix.h5ad`). If not provided, `gRNA_counts_table_csv` file prefix is used.