# Feeding custom prior into `bean run`

In this tutorial, we consider a case where user want to specify per-variant prior beliefs into the model using `--prior-params=prior_params.pkl` for `bean run`.  
In particular, we see a case where the library is split into two sublibraries, where guides are disjoint but many guides in two sublibraries may target shared variant. In this particular case, we provide helper command `bean build-prior`.

## Example workflow
```bash
screen_id1=var_mini_screen_sub1
screen_id2=var_mini_screen_sub2
working_dir=tests/data/
output_dir=tests/test_res/

# 1. Run sublibrary 1
# It is important that you feed --save-raw so the run output will contain the input for the next step.
bean run sorting variant $working_dir/${screen_id1}.h5ad --control-condition D14_1 -o $working_dir --fit-negctrl --save-raw

# 2. Build prior
# Feed first and second `bean run` scripts, and output file pickle file path that will store prior_params
bean build-prior \
'bean run sorting variant $working_dir/${screen_id1}.h5ad --control-condition D14_1 -o $working_dir --fit-negctrl --save-raw' \
'bean run sorting variant $working_dir/${screen_id2}.h5ad --control-condition D14_2 -o $working_dir --fit-negctrl' \
$working_dir/bean_run_result.${screen_id1}/MixtureNormal.result.pkl \
$working_dir/prior_params.pkl

# 3. Run sublibrary 2 with the specified prior
# Feed in the prior_param.pkl file from the previous step.
bean run sorting variant $working_dir/${screen_id1}.h5ad --control-condition D14_2 -o BE_part1_variant_masked_sorting --fit-negctrl --prior-param $working_dir/prior_params.pkl
```

## Manually specifying prior
TBD. If this function is desired please open an issue for expedited documentation :)