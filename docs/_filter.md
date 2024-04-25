# `filter`: Filtering (and optionally translating) alleles
As `tiling` mode of `bean run` accounts for any robustly observed alleles, `bean filter` filters for such alleles.
```bash
bean filter my_sorting_screen_masked.h5ad \
-o my_sorting_screen_filtered.h5ad  `# Output file path` \
```

# Output
Above command produces 
* `my_sorting_screen_filtered.h5ad` with filtered alleles stored in `.uns`,   
* `my_sorting_screen_filtered.filtered_allele_stats.pdf`, and `my_sorting_screen_filtered.filter_log.txt` that report allele count stats in each filtering step.  

You may want to adjust the flitering parameters to obtain optimal balance between # guides per variant & # variants that are scored. See example outputs of filtering step [here](docs/example_filtering_output/).


# Translating alleles
If you want to obtain **amino acid level variant** for coding sequence tiling screens, provide coding sequence positions which variants occuring within the coding sequence will be translated. *This is optional, but **highly recommended** to increase per-(coding)variant support.*  

<img src="/crispr-bean/assets/translation.png" alt="Allele translation" width="500"/>  
  

```bash
bean filter my_sorting_screen.h5ad \
-o my_sorting_screen_masked.h5ad \
--translate   `# Translate coding variants` \
[ --translate-gene-name GENE_SYMBOL OR
  --translate-genes-list path_to_gene_names_file.txt OR
  --translate-fasta gene_exon.fa, OR
  --translate-fastas-csv gene_exon_fas.csv]
```
* When library covers a single gene, do either of the following:
  1. Feed `--translate-gene-name GENE_SYMBOL` if your `genomic_pos` column of `sgRNA_info_tbl` is compatible with [MANE transcript](https://useast.ensembl.org/info/genome/genebuild/mane.html)'s reference genome. (Per 10/23/2023, GRCh38). This will automatically load the exon positions based on MANE transcript annotation.
  2. To use your custom coding sequence and exon positions, feed `--translate-fasta gene_exon.fa` argument where `gene_exon.fa` is the FASTA file with entries of exons. [See full details here](https://github.com/pinellolab/crispr-bean/blob/main/docs/exon_fa_format.md).
* When library covers multiple genes, do either of the following:  
  1. Feed `--translate-genes-list path_to_gene_names_file.txt` where `path_to_gene_names_file.txt` is file with one gene symbol per line.
  2. Feed `--translate-fastas-csv gene_exon_fas.csv` where `gene_exon_fas.csv` is the csv file with lines `gene_id,gene_exon_fasta_path` without header. Each FASTA file in `gene_exon_fasta_path` is formatted [as the single-gene FASTA file](https://github.com/pinellolab/crispr-bean/blob/main/docs/exon_fa_format.md).
* Translation will keep the variants outside the coding sequence as nucleotide-level variants, while aggregating variants leading to the same coding sequence variants.

# Getting splice sites
We provide the utility script to obtain the splice sites if you use the MANE transcript with gene symbol (`--translate-gene-name GENE_SYMBOL` or `--translate-genes-list path_to_gene_names_file.txt`).
```bash
bean get-splice-sites LDLR A LDLR_splice_sites.csv --gene-name
```