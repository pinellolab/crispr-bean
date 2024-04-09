import argparse
import os


def parse_args(parser=None):
    if parser is None:
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
        default="condition",
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
    parser.add_argument(
        "--save-fig",
        action="store_true",
        help="Save .pdf of the figures included in the report.",
    )

    return parser


def check_args(args):
    if args.output_prefix is None:
        sample_id = os.path.splitext(os.path.basename(args.bdata_path))[0]
        args.output_prefix = (
            f"{os.path.dirname(args.bdata_path)}/bean_profile.{sample_id}/{sample_id}"
        )
        os.makedirs(args.output_prefix, exist_ok=True)
    if args.window_length < 1:
        raise ValueError(f"window_length {args.window_length} is too small.")
    if args.window_length > 20:
        raise ValueError(f"window_length {args.window_length} is too large.")
    return args
