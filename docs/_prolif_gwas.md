## Variant survival screen tutorial
GWAS variant screen with per-variant gRNA tiling design, selected based on FACS signal quantiles.  

<table>
  <tr>
    <th>Library design</th>
    <td>Variant (gRNAs tile each target variant)   <br> <img src="/crispr-bean/assets/variant.png" alt="variant library design" width="600"/></td>
  </tr>
  <tr>
    <th>Selection</th>
    <td>Cells are sorted based on FACS signal quantiles  <br>  <img src="/crispr-bean/assets/proliferation.png" alt="variant library design" width="300"/></td>
  </tr>
</table>

<br></br>

## Example workflow
```bash
screen_id=my_sorting_tiling_screen
working_dir=my_workdir

# 1. Count gRNA & reporter
bean count-samples \
--input ${working_dir}/sample_list.csv    `# Contains fastq file path; see test file for example.`\
-b A                                  `# Base A is edited (into G)` \
-f ${working_dir}/test_guide_info.csv     `# Contains gRNA metadata; see test file for example.`\
-o ${working_dir}                                 `# Output directory` \
-r                                    `# Quantify reporter edits` \
-n ${screen_id}                          `# ID of the screen to be counted`   

# 2. QC samples & guides
bean qc \
  ${working_dir}/bean_count_${screen_id}.h5ad             `# Input ReporterScreen .h5ad file path` \
  -o ${working_dir}/bean_count_${screen_id}_masked.h5ad   `# Output ReporterScreen .h5ad file path` \
  -r ${working_dir}/qc_report_${screen_id}                `# Prefix for QC report` \
  --lfc-conds D0,D14                `# Conditions to calculate LFC of positive controls` \
  -b                                       ` # Remove replicates with no good samples.

# 3. Quantify variant effect
bean run survival variant \
    ${working_dir}/bean_count_${screen_id}_masked.h5ad \
    -o ${working_dir}/ \
    --fit-negctrl \
    --scale-by-acc \
    --accessibility-col accessibility
```
See more details below.

## 1. Count gRNA & reporter (:ref:`count_samples`)
```bash
bean count-samples \
--input ${working_dir}/sample_list.csv    `# Contains fastq file path; see test file for example.`\
-b A                                  `# Base A is edited (into G)` \
-f ${working_dir}/test_guide_info.csv     `# Contains gRNA metadata; see test file for example.`\
-o ${working_dir}                                 `# Output directory` \
-r                                    `# Quantify reporter edits` \
-n ${screen_id}                          `# ID of the screen to be counted`  
```
Make sure you follow the [input file format](../../README#input-file-format) for seamless downstream steps. This will produce `./bean_count_${screen_id}.h5ad`. 

## 2. QC samples & guides (:ref:`qc`)
Base editing data will include QC about editing efficiency. As QC uses predefined column names and values, beware to follow the [input file guideline](../../README#input-file-format), but you can change the parameters with the full argument list of [`bean qc`](../../README#bean qc-qc-of-reporter-screen-data). (Common factors you may want to tweak is `--ctrl-cond=bulk` and `--lfc-conds=top,bot` if you have different sample condition labels.)
```
bean qc \
  ${working_dir}/bean_count_${screen_id}.h5ad             `# Input ReporterScreen .h5ad file path` \
  -o ${working_dir}/bean_count_${screen_id}_masked.h5ad   `# Output ReporterScreen .h5ad file path` \
  -r ${working_dir}/qc_report_${screen_id}                `# Prefix for QC report` \
  --lfc-conds D0,D14                `# Conditions to calculate LFC of positive controls` \
  -b                                       ` # Remove replicates with no good samples.

```



If the data does not include reporter editing data, you can provide `--no-editing` flag to omit the editing rate QC.


## 3. Quantify variant effect (:ref:`run`)

`bean run` can take 3 run options to quantify editing rate:  
1. From **reporter + accessibility**  
    If your gRNA metadata table (`${working_dir}/test_guide_info.csv` above) included per-gRNA accessibility score, 
    ```
    bean run sorting variant \
    ${working_dir}/bean_count_${screen_id}_masked.h5ad \
    -o $working_dir \
    --fit-negctrl \
    --scale-by-acc \
    --accessibility-col accessibility
    ```
    If your gRNA metadata table (`${working_dir}/test_guide_info.csv` above) included per-gRNA chromosome & position and you have bigWig file with accessibility signal, 
    ```
    bean run sorting variant \
    ${working_dir}/bean_count_${screen_id}_masked.h5ad \
    -o $working_dir \
    --fit-negctrl \
    --scale-by-acc \
    --accessibility-bw accessibility.bw
    ```

2. From **reporter**, without accessibility

    This assumes the all target sites have the uniform chromatin accessibility.
    ```
    bean run sorting variant \
    ${working_dir}/bean_count_${screen_id}_masked.h5ad \
    -o $working_dir \
    --fit-negctrl 
    ```
3. No reporter information, assume the same editing efficiency of all gRNAs.  
    Use this option if your data don't have editing outcome information.
    ```
    bean run sorting variant \
    ${working_dir}/bean_count_${screen_id}_masked.h5ad \
    -o $working_dir \
    --fit-negctrl \
    --uniform-edit
    ```