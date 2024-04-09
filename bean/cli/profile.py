#!/usr/bin/env python
import os
import papermill as pm
import bean as be
from bean.plotting.utils import parse_args, check_args


def main(args):
    print("  \n~~~BEAN Profile~~~")
    print("-Profile editing patterns of your editor-")
    print(
        r"""
    _ _                     __ _ _     
  /  \ '\     _ __ _ _ ___ / _(_) |___ 
  |   \  \   | '_ \ '_/ _ \  _| | / -_)
   \   \  |  | .__/_| \___/_| |_|_\___|
    `.__|/   |_|                       
    """
    )
    args = check_args(args)
    os.system(
        "python -m ipykernel install --user --name bean_python3 --display-name bean_python3"
    )
    pm.execute_notebook(
        f"{os.path.dirname(be.__file__)}/notebooks/profile_editing_preference.ipynb",
        f"{args.output_prefix}.ipynb",
        parameters=dict(
            bdata_path=args.bdata_path,
            output_prefix=args.output_prefix,
            replicate_col=args.replicate_col,
            condition_col=args.condition_col,
            control_condition=args.control_condition,
            max_editing_window_length=args.window_length,
            pam_col=args.pam_col,
        ),
        kernel_name="bean_python3",
    )
    os.system(
        f"jupyter nbconvert --to html {args.output_prefix}_editing_preference.ipynb"
    )
