BEAN stores mapped gRNA and allele counts in `ReporterScreen` object which is compatible with [AnnData](https://anndata.readthedocs.io/en/latest/index.html).  

<img src="../imgs/data_structure_v2.png" alt="ReporterScreen object structure" width="700" />

  * `.guides`: guide information provided in input (`gRNA_library.csv` in above example)
  * `.samples`: sample information provided in input (`sample_list.csv` in above example)
  * `.X`: Main guide count matrix, where row corresponds to each guide in `.guides` and columns correspond to samples in `.samples`.
Following attributes are included if matched reporter is provided and you chose to read edit/allele information from the reporter using `-r` option.
  * `.X_bcmatch [Optional]`: Contains information about number of barcode-matched reads. Information about R2 barcode should be specified as `barcode` column in your `gRNA_library.csv` file.
  * `.X_edits [Optional]`: If target position of each guide is specified as `target_pos` in input `gRNA_library.csv` file and `--match-target-position` option is provided, the result has the matrix with the number of target edit at the specified positions.
  * `.allele_tables [Optional]`: Dictionary with a single allele count table that counts per guide and allele combination, what is the count per sample. 