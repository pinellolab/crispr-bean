# Tiling sorting screen tutorial
Tiling screen that tiles gRNA densely across locus or multiple loci, selected based on FACS signal quantiles.  

<table>
  <tr>
    <th>Library design</th>
    <td>Tiling (gRNAs tile each locus densely)   <br> <img src="../../imgs/tiling.png" alt="tiling library design" width="300"/> </td>
  </tr>
  <tr>
    <th>Selection</th>
    <td>Cells are sorted based on FACS signal quantiles  <br>  <img src="../../imgs/sorting_bins@8x.png" alt="variant library design" width="300"/></td>
  </tr>
</table>

<br></br>

## 1. Count gRNA & reporter ([`bean-count-samples`](../../README#bean-count-samples-count-reporter-screen-data))
```
bean-count-samples \
--input tests/data/sample_list_tiling.csv          `# Contains fastq file path; see test file for example.`\
-b A                                               `# Base A is edited (into G)` \
-f tests/data/test_guide_info_tiling_chrom.csv     `# Contains gRNA metadata; see test file for example.`\
-o tests/test_res/ \
-r                                                 `# Quantify reporter edits` \
--tiling
```
Make sure you follow the [input file format](../../README#input-file-format) for seamless downstream steps. This will produce `tests/test_res/bean_count_sample_list.h5ad`. 

## 2. QC ([`bean-qc`](../../README#bean-qc-qc-of-reporter-screen-data))
Base editing data will include QC about editing efficiency. As QC uses predefined column names and values, beware to follow the [input file guideline](../../README#input-file-format), but you can change the parameters with the full argument list of [`bean-qc`](../../README#bean-qc-qc-of-reporter-screen-data). (Common factors you may want to tweak is `--ctrl-cond=bulk` and `--lfc-conds=top,bot` if you have different sample condition labels.)
```
bean-qc \
  my_sorting_screen.h5ad           `# Input ReporterScreen .h5ad file path` \
  -o my_sorting_screen_masked.h5ad `# Output ReporterScreen .h5ad file path` \
  -r qc_report_my_sorting_screen   `# Prefix for QC report` \
  [--tiling]                       `# Not required if you have passed --tiling in counting step`
```



If the data does not include reporter editing data, you can provide `--no-editing` flag to omit the editing rate QC.

## 3. Filter alleles ([`bean-filter`](../../README#bean-filter-filtering-and-optionally-translating-alleles))
As tiling library doesn't have designated per-gRNA target variant, any base edit observed in reporter may be the candidate variant, while having too many variants with very low editing rate significantly decreases the power. Variants are filtered based on multiple criteria in `bean-fitler`.  

If the screen targets coding sequence, it's beneficial to translate edits into coding varaints whenever possible for better power. For translation, provide `--translate` and one of the following:
```
[ --translate-gene-name GENE_SYMBOL OR
  --translate-genes-list path_to_gene_names_file.txt OR
  --translate-fasta gene_exon.fa, OR
  --translate-fastas-csv gene_exon_fas.csv]
```
where `path_to_gene_names_file.txt` has one gene symbol per line, and gene symbol uses its MANE transcript (hg38) coordinates of exons. In order to use other reference versions or transcript ID, you'll need to feed in fasta file. See detailed formatting of fasta file [here](../../README#translating-alleles).

Example allele filtering given we're translating based on MANE transcript exons of multiple gene symbols:

```bash
bean-filter tests/data/tiling_mini_screen_masked.h5ad \
-o tests/data/tiling_mini_screen_annotated \
--filter-target-basechange                             `# Filter based on intended base changes. If -b A was provided in bean-count, filters for A>G edit. If -b C was provided, filters for C>T edit.`\
--filter-window --edit-start-pos 0 --edit-end-pos 19   `# Filter based on editing window in spacer position within reporter.`\
--filter-allele-proportion 0.1 --filter-sample-proportion 0.3 `#Filter based on allele proportion larger than 0.1 in at least 0.3 (30%) of the control samples.` \
--translate --translate-genes-list tests/data/gene_symbols.txt
```

Ouptut file `` shows number of alleles per guide and number of guides per variant, where we want high enough values for the latter. See the typical output for dataset with good editing coverage & filtering result [here](../example_filtering_ouptut/).

## 4. Quantify variant effect ([`bean-run`](../../README#bean-run-quantify-variant-effects))
By default, `bean-run [sorting,survival] tiling` uses most filtered allele counts table for variant identification and quantification of their effects. **Check [allele filtering output](../example_filtering_ouptut/)** and choose alternative filtered allele counts table if necessary.   

`bean-run` can take 3 run options to quantify editing rate:  
1. From **reporter + accessibility**  
    If your gRNA metadata table (`tests/data/test_guide_info.csv` above) included per-gRNA accessibility score, 
    ```
    bean-run sorting tiling \
    tests/data/tiling_mini_screen_annotated.h5ad \
    -o tests/test_res/var/ \
    --fit-negctrl \
    --scale-by-acc \
    --accessibility-col accessibility
    ```
    If your gRNA metadata table (`tests/data/test_guide_info.csv` above) included per-gRNA chromosome & position and you have bigWig file with accessibility signal, 
    ```
    bean-run sorting tiling \
    tests/data/tiling_mini_screen_annotated.h5ad \
    -o tests/test_res/var/ \
    --fit-negctrl \
    --scale-by-acc \
    --accessibility-bw accessibility.bw
    ```

2. From **reporter**
    ```
    bean-run sorting tiling \
    tests/data/tiling_mini_screen_annotated.h5ad \
    -o tests/test_res/var/ \
    --fit-negctrl 
    ```
3. No reporter information, assume the same editing efficiency of all gRNAs.  
    Use this option if your data don't have editing rate information.
    ```
    bean-run sorting tiling \
    tests/data/tiling_mini_screen_annotated.h5ad \
    -o tests/test_res/var/ \
    --fit-negctrl \
    --uniform-edit
    ```

See [full argument list](../../README#optional-parameters) to accommodate different input sample & guide metadata columns/values and run options.