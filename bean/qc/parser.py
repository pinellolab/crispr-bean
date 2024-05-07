import argparse


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def parse_args(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument(
        "bdata_path", help="Path to the ReporterScreen object to run QC on", type=str
    )
    thres_parser = parser.add_argument_group("QC thresholds")
    run_parser = parser.add_argument_group("Run options")
    input_parser = parser.add_argument_group("Input .h5ad formatting")

    thres_parser.add_argument(
        "--count-correlation-thres",
        help="Correlation threshold to mask out.",
        type=float,
        default=0.7,
    )
    thres_parser.add_argument(
        "--edit-rate-thres",
        help="Mean editing rate threshold per sample to mask out.",
        type=float,
        default=0.1,
    )
    thres_parser.add_argument(
        "--lfc-thres",
        help="Positive guides' correlation threshold to filter out.",
        type=float,
        default=-0.1,
    )

    parser.add_argument(
        "-o",
        "--out-screen-path",
        help="Path where quality-filtered ReporterScreen object to be written to",
        type=str,
    )
    parser.add_argument(
        "-r",
        "--out-report-prefix",
        help="Output prefix of qc report (prefix.html, prefix.ipynb)",
        type=str,
    )

    run_parser.add_argument(
        "-b",
        "--remove-bad-replicates",
        help="Remove replicates with at least two of its samples meet the QC threshold.",
        action="store_true",
    )
    run_parser.add_argument(
        "-i",
        "--ignore-missing-samples",
        help="If the flag is not provided, if the ReporterScreen object does not contain all condiitons for each replicate, make fake empty samples. If the flag is provided, don't add dummy samples.",
        action="store_true",
    )
    run_parser.add_argument(
        "--no-editing",
        help="Ignore QC about editing. Can be used for QC of other editing modalities.",
        action="store_true",
    )
    run_parser.add_argument(
        "--dont-recalculate-edits",
        help="When ReporterScreen.layers['edit_count'] exists, do not recalculate the edit counts from ReporterScreen.uns['allele_count'].",
        action="store_true",
    )

    input_parser.add_argument(
        "--tiling",
        dest="tiling",
        type=str2bool,
        help="Specify that the guide library is tiling library without 'n guides per target' design",
    )
    input_parser.add_argument(
        "--replicate-col",
        help="Label of column in `bdata.samples` that describes replicate ID.",
        type=str,
        default="replicate",
    )
    input_parser.add_argument(
        "--sample-covariates",
        help="Comma-separated list of column names in `bdata.samples` that describes non-selective experimental condition. (drug treatment, etc.)",
        type=str,
        default=None,
    )
    input_parser.add_argument(
        "--condition-col",
        help="Label of column in `bdata.samples` that describes experimental condition. (sorting bin, time, etc.)",
        type=str,
        default="condition",
    )
    input_parser.add_argument(
        "--target-pos-col",
        help="Target position column in `bdata.guides` specifying target edit position in reporter",
        type=str,
        default="target_pos",
    )
    input_parser.add_argument(
        "--rel-pos-is-reporter",
        help="Specifies whether `edit_start_pos` and `edit_end_pos` are relative to reporter position. If `False`, those are relative to spacer position.",
        action="store_true",
        default=False,
    )
    input_parser.add_argument(
        "--edit-start-pos",
        help="Edit start position to quantify editing rate on, 0-based inclusive.",
        default=2,
        type=int,
    )
    input_parser.add_argument(
        "--edit-end-pos",
        help="Edit end position to quantify editing rate on, 0-based exclusive.",
        default=7,
        type=int,
    )

    input_parser.add_argument(
        "--posctrl-col",
        help="Column name in ReporterScreen.guides DataFrame that specifies guide category. To use all gRNAs, feed empty string ''.",
        type=str,
        default="target_group",
    )
    input_parser.add_argument(
        "--posctrl-val",
        help="Value in ReporterScreen.guides[`posctrl_col`] that specifies guide will be used as the positive control in calculating log fold change.",
        type=str,
        default="PosCtrl",
    )

    input_parser.add_argument(
        "--lfc-conds",
        help="Values in of column in `ReporterScreen.samples[condition_label]` for LFC will be calculated between, delimited by comma",
        type=str,
        default="top,bot",
    )
    input_parser.add_argument(
        "--control-condition",
        help="Values in of column in `ReporterScreen.samples[condition_label]` for guide-level editing rate to be calculated",
        type=str,
        default="bulk",
    )
    parser.add_argument(
        "--reporter-length",
        type=int,
        default=None,
        help="Length of reporter sequence in the construct.",
    )
    parser.add_argument(
        "--reporter-right-flank-length",
        type=int,
        default=None,
        help="Length of the right-flanking nucleotides of protospacer in the reporter.",
    )
    return parser
