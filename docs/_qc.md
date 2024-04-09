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

<img src="/crispr-bean/assets/qc_output.png" alt="Allele translation" width="900"/>  

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
