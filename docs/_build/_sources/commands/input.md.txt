This document describes the input files of :ref:`count_samples`.
## sgRNA_info_table.csv
File should contain following columns. 
* `name`: gRNA ID column
* `sequence`: gRNA sequence
* `barcode`: R2 barcode to help match reporter to gRNA, written in the sense direction (as in R1)
* In order to use accessibility in the [variant effect quantification](#bean-run-quantify-variant-effects), provide accessibility information in one of two options. (For non-targeting guides, provide NA values (empty cell).)   
  * Option 1: `chrom` & `genomic_pos`: Chromosome (ex. `chr19`) and genomic position of guide sequence. You will have to provide the path to the bigwig file with matching reference version in `bean run`. 
  * Option 2: `accessibility_signal`: ATAC-seq signal value of the target loci of each guide.  
* For variant library (gRNAs are designed to target specific variants and ignores bystander edits)
  * `target`: This column denotes which target variant/element of each gRNA. This is not used in `bean count[-samples]` but required to run `bean run` in later steps.
  * `target_group`: If negative/positive control gRNA will be considered in `bean qc` and/or `bean run`, specify as "NegCtrl"/"PosCtrl" in this column. 
  * `target_pos`: If `--match_target_pos` flag is used, input file needs `target_pos` which specifies 0-based relative position of targeted base within Reporter sequence.
* For tiling library (gRNAs tile coding / noncoding sequences)
  * `strand`: Specifies gRNA strand information relative to the reference genome.
  * `chrom`: Chromosome of gRNA targeted locus.
  * `start_pos`: gRNA starting position in the genome. Required when you provide `strand` column. Should specify the smaller coordinate value among start and end position regardless of gRNA strandedness.

Also see examples for [variant library](tests/data/test_guide_info.csv) and [tiling library](tests/data/test_guide_info_tiling.csv).

## sample_list.csv
File should contain following columns with header.
* `R1_filepath`: Path to read 1 `.fastq[.gz]` file
* `R2_filepath`: Path to read 1 `.fastq[.gz]` file
* `sample_id`: ID of sequencing sample
* `replicate`: Replicate # of this sample (Should NOT contain `.`)
* `condition`: Name of the sorting bin (ex. `top`, `bot`), or label of timepoint (ex. `D5`, `D18`)  

For FACS sorting screens:
* `upper_quantile`: FACS sorting upper quantile
* `lower_quantile`: FACS sorting lower quantile  

For proliferation / survival screens:
* `time`: Numeric time following the base editing of each sample.


Also see examples for [FACS sorting screen](tests/data/sample_list.csv).