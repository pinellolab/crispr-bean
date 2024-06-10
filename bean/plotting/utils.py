import argparse
import os
import bean as be


def check_args(args):
    if args.output_prefix is None:
        sample_id = os.path.splitext(os.path.basename(args.bdata_path))[0]
        sample_path = os.path.dirname(args.bdata_path)
        if not sample_path:
            sample_path = "."
        args.output_prefix = f"{sample_path}/bean_profile.{sample_id}/bean_profile"
        os.makedirs(f"{sample_path}/bean_profile.{sample_id}/", exist_ok=True)
    if args.window_length < 1:
        raise ValueError(f"window_length {args.window_length} is too small.")
    cdata = be.read_h5ad(args.bdata_path)
    cdata.samples["replicate"] = cdata.samples[args.replicate_col]
    cdata_bulk = cdata[:, cdata.samples[args.condition_col] == args.control_condition]
    if len(cdata_bulk) == 0:
        raise ValueError(
            f"No samples with bdata.samples['{args.condition_col}'] == {args.control_condition}. Please check your input arguments --condition-col & --control-condition."
        )
    if args.pam_col is not None and args.pam_col not in cdata.guides.columns:
        raise ValueError(
            f"Specified --pam-col `{args.pam_col}` does not exist in ReporterScreen.guides.columns ({cdata.guides.columns}). Please check your input. If you don't want PAM output, please do not provide --pam-col argument.s"
        )
    if args.replicate_col not in cdata.samples.columns:
        raise ValueError(
            f"Specified --replicate-col `{args.replicate_col}` does not exist in ReporterScreen.samples.columns ({cdata.samples.columns}). Please check your input."
        )
    if args.condition_col not in cdata.samples.columns:
        raise ValueError(
            f"Specified --condition-col `{args.condition_col}` does not exist in ReporterScreen.samples.columns ({cdata.samples.columns}). Please check your input."
        )
    if args.control_condition not in cdata.samples[args.condition_col].tolist():
        raise ValueError(
            f"Specified --control-condition `{args.control_condition}` does not exist in ReporterScreen.samples[{args.condition_col}] ({cdata.samples[args.condition_col]}). Please check your input."
        )
    return args
