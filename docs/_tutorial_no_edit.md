# Sorting screen tutorial without reporter edits
Screens without editing outcome, where each gRNA is assigned to a target.  

<table>
  <tr>
    <th>Library design</th>
    <td>Each guide RNA (or other perturbation such as ORF) has a specified target. </td>
  </tr>
  <tr>
    <th>Selection</th>
    <td>Cells are sorted based on FACS signal quantiles  <br>  <img src="/crispr-bean/assets/sorting_bins@8x.png" alt="sorting_bins" width="300"/></td>
  </tr>
</table>

<br></br>
Here, we consider an example where BEAN uses the **external count** without reporter, gRNA and sample information to infer the per-target effect sizes.

## Example workflow
```bash
screen_id=my_sorting_screen
working_dir=my_workdir

# 1. Given that we have gRNA count for each sample, generate ReporterScreen (.h5ad) object for downstream analysis.
bean create-screen ${working_dir}/gRNA_info.csv ${working_dir}/sample_list.csv ${working_dir}/gRNA_counts.csv -o ${working_dir}/bean_count_${screen_id}

# 2. QC samples & guides
bean qc \
  ${working_dir}/bean_count_${screen_id}.h5ad             `# Input ReporterScreen .h5ad file path` \
  -o ${working_dir}/bean_count_${screen_id}_masked.h5ad   `# Output ReporterScreen .h5ad file path` \
  -r ${working_dir}/qc_report_${screen_id}                `# Prefix for QC report` \
  -b                                       ` # Remove replicates with no good samples.

# 3. Quantify variant effect
bean run sorting variant \
    ${working_dir}/bean_count_${screen_id}_masked.h5ad \
    -o ${working_dir}/ \
    --uniform-edit             `# As we have no edit information.` \
    [--no-negative-control]    `# If you don't have the negative control gRNAs.`
```

See the example input files [here](https://pinellolab.github.io/crispr-bean/create_screen.html).