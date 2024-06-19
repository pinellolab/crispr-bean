# Model run

```
bean-run [run_mode] screen_data_path [options]
```

## Input file format
### Shared
#### screen.guides
* `sequence`
* `reporter_sequence`
##### optional columns for models with accessibility
* `chr`, `genomic_pos`
#### screen.samples
* `rep`: Replicate column showing experimental replicate.
* `condition`: Condition column showing the conditions that are selected in the screen. Can be substituted for different string than `condition`, and it will be specified during the `bean-run` as `--condition-column your_col_name`.

### Variant screens
