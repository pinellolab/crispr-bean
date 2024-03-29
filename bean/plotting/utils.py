import argparse


def parse_args():
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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bdata_path", help="Path to the ReporterScreen object to run QC on", type=str
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        help="Output prefix of editing pattern report (prefix.html, prefix.ipynb). If not provided, base name of `bdata_path` is used.",
        type=str,
    )
    parser.add_argument(
        "--replicate-col",
        help="Column name in `bdata.samples` that describes replicate ID.",
        type=str,
        default="replicate",
    )
    parser.add_argument(
        "--condition-col",
        help="Column name in `bdata.samples` that describes experimental condition. (sorting bin, time, etc.)",
        type=str,
        default="bin",
    )
    parser.add_argument(
        "--pam-col",
        help="Column name describing PAM of each gRNA in `bdata.guides`.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--control-condition",
        help="Control condition where editing preference would be profiled at. Pre-filters data where `bdata.samples[condition_col] == control_condition`.",
        type=str,
        default="bulk",
    )
    parser.add_argument(
        "-w",
        "--window-length",
        help="Window length of editing window of maximal editing efficiency to be identified. This window is used to quantify context specificity within the window.",
        type=int,
        default=6,
    )

    args = parser.parse_args()
    if args.output_prefix is None:
        args.output_prefix = f"{args.bdata_path.rsplit('.h5ad', 1)[0]}"
    return args


def check_args(args):
    if args.window_length < 1:
        raise ValueError(f"window_length {args.window_length} is too small.")
    if args.window_length > 20:
        raise ValueError(f"window_length {args.window_length} is too large.")
    return args
