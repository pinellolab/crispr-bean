# `bean run` output files
These are the example output of [`bean run`](https://pinellolab.github.io/crispr-bean/run.html).
`bean run` produces [`bean_element_result.[model_type].csv`] as the main output and will produce log and per-guide information files, and intermediate `.pkl` file that saves ELBO loss for variational inference.

## `bean_element_result.[model_type].csv`
- Variant ID / grouping
  - `edit`: Variant ID.
  - `group`: The grouping of the coding variants, assigned as one of nonsense/missense/synonymous.
  - `int_pos`: The integer position of the noncoding variants.
  - `chrom`: The chromosome of the variant.
  - `pos`: The position of the variant. If coding variant, starts with `A` and the position specified 1-based amino acid position. If noncodig variant, numeric genomic position.
  - `ref`: The reference base/amino acid of the variant.
  - `alt`: The alternative base/amino acid of the variant.
  - `coding`: A flag indicating if the element is coding variant or not.

- Per-variant summary of variant-producing guides
  - `guide_target_group`: Aggregated `target_group` column in the input sgRNA_info.csv file. All unique values of the guides that produced (filtered) edited alleles that includes this variant is listed.
  - `effective_edit_rate`: The effective editing rate of the element. Calculated as `sum_over_guides(sum_over_alleles(per_guide_allele_editing_rate / # variants in the allele))`.
  - `editing_guides`: List of guides that edited the variant.
  - `per_guide_editing_rates`: The per-guide editing rates of the variant.
  - `n_guides`: The number of guides that edited the variant.
  - `n_coocc`: The number of unique co-occurring variants that appeared together in any alleles that contains the variant.
- Variant effect size: Use `mu_z_adj` whenever available, otherwise `mu_z_scaled`, otherwise `mu_z`.
  - `mu`: The mean value of the variant effect size.
  - `mu_sd`: The standard deviation of the mean value of the variant effect size.
  - `mu_z`: The z-score of the mean value of the variant effect size.
  - `sd`: The standard deviation of the phenotype induced by the variant.
  - `CI[0.025,0.975]`: The 95% credible interval of the mean value of the variant effect size. Corresponds to `mu_z_adj` when available, otherwise `mu_z_scaled`, otherwise `mu_z`. 
  - `[]_scaled`: Above values scaled by negative control variants.
  - `[]_adj`: Above values scaled by synonymous variants.
  