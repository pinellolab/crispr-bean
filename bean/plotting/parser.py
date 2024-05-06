import argparse


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
    parser.add_argument(
        "--reporter-length",
        type=int,
        default=32,
        help="Length of reporter sequence in the construct.",
    )
    parser.add_argument(
        "--reporter-right-flank-length",
        type=int,
        default=6,
        help="Length of the right-flanking nucleotides of protospacer in the reporter.",
    )
    return parser
