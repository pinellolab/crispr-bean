# Variant sorting screen tutorial
GWAS variant screen with per-variant gRNA tiling design, selected based on FACS signal quantiles.  

<table>
  <tr>
    <th>Library design</th>
    <td>Variant (gRNAs tile each target variant)   <br> <img src="/crispr-bean/assets/variant.png" alt="variant library design" width="600"/></td>
  </tr>
  <tr>
    <th>Selection</th>
    <td>Cells are sorted based on FACS signal quantiles  <br>  <img src="/crispr-bean/assets/sorting_bins@8x.png" alt="variant library design" width="300"/></td>
  </tr>
</table>

<br></br>

## Example workflow
```bash
screen_id=var_mini_screen
working_dir=tests/data/

# 1. Count gRNA & reporter
bean count-samples \
--input ${working_dir}/sample_list.csv    `# Contains fastq file path; see test file for example.`\
-b A                                  `# Base A is edited (into G)` \
-f ${working_dir}/test_guide_info.csv     `# Contains gRNA metadata; see test file for example.`\
-o ./                                 `# Output directory` \
-r                                    `# Quantify reporter edits` \
-n ${screen_id}                          `# ID of the screen to be counted`   
# count-samples output from above test run is too low in read depth. Downstream processes can be run with test file included in the Github repo.

# (Optional) Profile editing patterns
bean profile tests/data/${screen_id}.h5ad --pam-col '5-nt PAM'

# 2. QC samples & guides
bean qc \
  ${working_dir}/bean_count_${screen_id}.h5ad             `# Input ReporterScreen .h5ad file path` \
  -o ${working_dir}/bean_count_${screen_id}_masked.h5ad   `# Output ReporterScreen .h5ad file path` \
  -r ${working_dir}/qc_report_${screen_id}                `# Prefix for QC report` \
  -b                                                      ` # Remove replicates with no good samples.

# 3. Quantify variant effect
bean run sorting variant \
    ${working_dir}/${screen_id}_masked.h5ad \
    -o ${working_dir}/ \
    --fit-negctrl \
    --scale-by-acc \
    --acc-bw-path tests/data/accessibility_signal_chr6.bw
```

See more details below.

## 1. Count gRNA & reporter (:ref:`count_samples`)
```bash
screen_id=my_sorting_tiling_screen

# 1. Count gRNA & reporter
bean count-samples \
--input ${working_dir}/sample_list.csv    `# Contains fastq file path; see test file for example.`\
-b A                                  `# Base A is edited (into G)` \
-f ${working_dir}/test_guide_info.csv     `# Contains gRNA metadata; see test file for example.`\
-o ./                                 `# Output directory` \
-r                                    `# Quantify reporter edits` \
-n ${screen_id}                          `# ID of the screen to be counted`   
```

Make sure you follow the [input file format](https://pinellolab.github.io/crispr-bean/input.html) for seamless downstream steps. This will produce `./bean_count_${screen_id}.h5ad`. 


## (Optional) Profile editing pattern (:ref:`profile`)
You can profile the pattern of base editing based on the allele counts. 

```bash
bean profile tests/data/${screen_id}.h5ad --pam-col '5-nt PAM'
```


## 2. QC samples & guides (:ref:`qc`)
Base editing data will include QC about editing efficiency. As QC uses predefined column names and values, beware to follow the [input file guideline](https://pinellolab.github.io/crispr-bean/input.html), but you can change the parameters with the full argument list of [bean qc](https://pinellolab.github.io/crispr-bean/qc.html). (Common factors you may want to tweak is `--ctrl-cond=bulk` and `--lfc-conds=top,bot` if you have different sample condition labels.)

```bash
bean qc \
  ${screen_id}.h5ad    `# Input ReporterScreen .h5ad file path` \
  -o ${screen_id}_masked.h5ad   `# Output ReporterScreen .h5ad file path` \
  -r qc_report_${screen_id}   `# Prefix for QC report` 
```


If the data does not include reporter editing data, you can provide `--no-editing` flag to omit the editing rate QC.


## 3. Quantify variant effect (:ref:`run`)

`bean run` can take 3 run options to quantify editing rate:  
1. From **reporter + accessibility**  
  If your gRNA metadata table (`${working_dir}/test_guide_info.csv` above) included per-gRNA accessibility score, 
    ```bash
    bean run sorting variant \
    ${working_dir}/${screen_id}_masked.h5ad \
    -o ${working_dir}/ \
    --fit-negctrl \
    --scale-by-acc \
    --accessibility-col accessibility
    ```

  If your gRNA metadata table (`${working_dir}/test_guide_info.csv` above) included per-gRNA chromosome & position and you have bigWig file with accessibility signal, 
    
    ```bash
    bean run sorting variant \
    ${working_dir}/${screen_id}_masked.h5ad \
    -o ${working_dir}/ \
    --fit-negctrl \
    --scale-by-acc \
    --accessibility-bw accessibility.bw
    ```

2. From **reporter**, without accessibility

  This assumes the all target sites have the uniform chromatin accessibility.

    ```bash
    bean run sorting variant \
    ${working_dir}/${screen_id}_masked.h5ad \
    -o ${working_dir}/ \
    --fit-negctrl 
    ```

3. No reporter information, assume the same editing efficiency of all gRNAs.  
  Use this option if your data don't have editing outcome information.

    ```bash
    bean run sorting variant \
    ${working_dir}/${screen_id}_masked.h5ad \
    -o ${working_dir}/ \
    --fit-negctrl \
    --uniform-edit
    ```

See the example output [here](https://github.com/pinellolab/crispr-bean/tree/main/tests/var_res_example/).
  