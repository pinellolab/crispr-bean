#!/usr/bin/env python
import os
import papermill as pm
import bean as be
from bean.qc.utils import check_args


def main(args):
    print("  \n~~~BEANQC~~~")
    print("-Check guide/sample level quality and mask / discard-")
    print(
        r"""
    _ _       
  /  \ '\        ___   ___ 
  |   \  \      / _ \ / __|
   \   \  |    | (_) | (__ 
    `.__|/      \__\_\\___|
    """
    )
    if args.out_screen_path is None:
        args.out_screen_path = f"{args.bdata_path.rsplit('.h5ad', 1)[0]}.masked.h5ad"
    if args.out_report_prefix is None:
        args.out_report_prefix = f"{args.bdata_path.rsplit('.h5ad', 1)[0]}.qc_report"
    args = check_args(args)
    os.system(
        "python -m ipykernel install --user --name bean_python3 --display-name bean_python3"
    )
    pm.execute_notebook(
        f"{os.path.dirname(be.__file__)}/notebooks/sample_quality_report.ipynb",
        f"{args.out_report_prefix}.ipynb",
        parameters=dict(
            bdata_path=args.bdata_path,
            out_bdata_path=args.out_screen_path,
            tiling=args.tiling,
            edit_quantification_start_pos=args.edit_start_pos,
            edit_quantification_end_pos=args.edit_end_pos,
            target_pos_col=args.target_pos_col,
            rel_pos_is_reporter=args.rel_pos_is_reporter,
            count_correlation_thres=args.count_correlation_thres,
            edit_rate_thres=args.edit_rate_thres,
            posctrl_col=args.posctrl_col,
            posctrl_val=args.posctrl_val,
            lfc_thres=args.lfc_thres,
            replicate_label=args.replicate_col,
            condition_label=args.condition_col,
            comp_cond1=args.lfc_cond1,
            comp_cond2=args.lfc_cond2,
            ctrl_cond=args.control_condition,
            exp_id=args.out_report_prefix,
            recalculate_edits=(not args.dont_recalculate_edits),
            base_edit_data=args.base_edit_data,
            remove_bad_replicates=args.remove_bad_replicates,
            reporter_length=args.reporter_length,
            reporter_right_flank_length=args.reporter_right_flank_length,
        ),
        kernel_name="bean_python3",
    )
    os.system(f"jupyter nbconvert --to=html {args.out_report_prefix}.ipynb")
