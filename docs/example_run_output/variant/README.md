# `bean run` output files
These are the example output of [`bean run`](https://pinellolab.github.io/crispr-bean/run.html).
`bean run` produces [`bean_element_result.[model_type].csv`] as the main output and will produce log and per-guide information files, and intermediate `.pkl` file that saves ELBO loss for variational inference.

## `bean_element_result.[model_type].csv`
- Variant ID / grouping
  - `target`: Target variant / element.
  - `target_group`: The group of the target variant / element, as specified in the input.
- Per-variant summary of variant-producing guides
  - `n_guides`: The number of guides targeting the variant, as described in the input.
  - `edit_rate_mean`: The mean editing rate of the variant.
  - `edit_rate_std`: The standard deviation of the editing rate of the variant.
- Variant effect size: Use `mu_z_adj` whenever available, otherwise `mu_z_scaled`, otherwise `mu_z`.
  - `mu`: The mean value of the variant effect size.
  - `mu_sd`: The standard deviation of the mean value of the variant effect size.
  - `mu_z`: The z-score of the mean value of the variant effect size.
  - `sd`: The standard deviation of the phenotype induced by the variant.
  - `CI[0.025,0.975]`: The 95% credible interval of the mean value of the variant effect size. Corresponds to `mu_z_adj` when available, otherwise `mu_z_scaled`, otherwise `mu_z`. 
  - `[]_scaled`: Above values scaled by negative control variants.
  - `[]_adj`: Above values scaled by negative control variants.

## `bean_sgRNA_result.[model_type].csv`
- `name`: sgRNA ID provided in the `name` column of the input.
- `edit_rate`: Editing rates
- `accessibility`: (Only if you have used `--scale-by-acc`) Accessibility signal that is used for scaling of the editing rate.
- `scaled_edit_rate`: (Only if you have used `--scale-by-acc`) Endogenous editing rate used for modeling, estimated by scaling reporter editing rate by accessibility signal
- `[replicate].[cond1]_[cond2]`: Raw per-replicate LFC with pseudocount fed in with `--guide-lfc-pseudocount` argument (default 5).