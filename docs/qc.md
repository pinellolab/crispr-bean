# `bean qc`: QC of reporter screen data
```bash
bean qc \
  my_sorting_screen.h5ad             `# Input ReporterScreen .h5ad file path` \
  -o my_sorting_screen_masked.h5ad   `# Output ReporterScreen .h5ad file path` \
  -r qc_report_my_sorting_screen     `# Prefix for QC report` \
  --ctrl-cond presort                `# "condition" column in the control sample before selection. Mean gRNA editing rates in these samples are reported. ` \
# Inspect the output qc_report_my_sorting_screen.html to tweak QC threshold

bean qc \
  my_sorting_screen.h5ad              \
  -o my_sorting_screen_masked.h5ad    \
  -r qc_report_my_sorting_screen      \
  #[--count-correlation-thres 0.7 ...]\
  -b
```

`bean qc` supports following quality control and masks samples with low quality. Specifically:  

<img src="assets/qc_output.png" alt="Allele translation" width="900"/>  

* Plots guide coverage and the uniformity of coverage
* Guide count correlation between samples
* Log fold change correlation when positive controls are provided
* Plots editing rate distribution
* Identify samples with low guide coverage/guide count correlation/editing rate and mask the sample in `bdata.samples.mask`
* Identify outlier guides to filter out

# Output
Above command produces 
* `my_sorting_screen_masked.h5ad` without problematic replicate and guides and with sample masks, and  
* `qc_report_my_sorting_screen.[html,ipynb]` as QC report.  
##### Optional arguments:
* `-o OUT_SCREEN_PATH`, `--out-screen-path OUT_SCREEN_PATH`
                        Path where quality-filtered ReporterScreen object to be written to
* `-r OUT_REPORT_PREFIX`, `--out-report-prefix OUT_REPORT_PREFIX`
                        Output prefix of qc report (prefix.html, prefix.ipynb)

##### QC thresholds:
* `--count-correlation-thres COUNT_CORRELATION_THRES`
                        Correlation threshold to mask out.
* `--edit-rate-thres EDIT_RATE_THRES`
                        Mean editing rate threshold per sample to mask out.
* `--lfc-thres LFC_THRES`
                        Positive guides' correlation threshold to filter out.

##### Run options:
* `-b`, `--remove-bad-replicates`
                        Remove replicates with at least two of its samples meet the QC threshold (bean run does not support having only one sorting bin sample for a replicate).
* `-i`, `--ignore-missing-samples`
                        If the flag is not provided, if the ReporterScreen object does not contain all condiitons for
                        each replicate, make fake empty samples. If the flag is provided, don't add dummy samples.
* `--no-editing`          Ignore QC about editing. Can be used for QC of other editing modalities.
* `--dont-recalculate-edits`
                        When ReporterScreen.layers['edit_count'] exists, do not recalculate the edit counts from
                        ReporterScreen.uns['allele_count'].

##### Input `.h5ad` formatting:
Note that these arguements will change the way the QC metrics are calculated for guides, samples, or replicates.
* `--tiling TILING`       Specify that the guide library is tiling library without 'n guides per target' design
* `--replicate-label REPLICATE_LABEL`
                        Label of column in `bdata.samples` that describes replicate ID.
* `--sample-covariates SAMPLE_COVARIATES`
                        Comma-separated list of column names in `bdata.samples` that describes non-selective
                        experimental condition. (drug treatment, etc.)
* `--condition-label CONDITION_LABEL`
                        Label of column in `bdata.samples` that describes experimental condition. (sorting bin, time,
                        etc.)
###### Editing rate calculation
  * `--control-condition CTRL_COND`
                        Values in of column in `ReporterScreen.samples[condition_label]` for guide-level editing rate
                        to be calculated. Default is `None`, which considers all samples.
  * `--rel-pos-is-reporter`
                        Specifies whether `edit_start_pos` and `edit_end_pos` are relative to reporter position. If
                        `False`, those are relative to spacer position.
  Editing rate is calculated with following parameters in 
    * Variant screens: 
      * `--target-pos-col TARGET_POS_COL`
                        Target position column in `bdata.guides` specifying target edit position in reporter
    * tiling screens:
      * `--edit-start-pos EDIT_START_POS`
                            Edit start position to quantify editing rate on, 0-based inclusive.
      * `--edit-end-pos EDIT_END_POS`
                            Edit end position to quantify editing rate on, 0-based exclusive.
###### LFC of positive controls
  * `--posctrl-col POSCTRL_COL`
                          Column name in ReporterScreen.guides DataFrame that specifies guide category. To use all
                          gRNAs, feed empty string ''.
  * `--posctrl-val POSCTRL_VAL`
                          Value in ReporterScreen.guides[`posctrl_col`] that specifies guide will be used as the
                          positive control in calculating log fold change.
  * `--lfc-conds LFC_CONDS`
                          Values in of column in `ReporterScreen.samples[condition_label]` for LFC will be calculated
                          between, delimited by comma