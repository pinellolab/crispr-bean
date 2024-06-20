import argparse


def parse_args(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="Generate prior_param.pkl for two batched runs, where two runs have no overlap but where guides that targeting a single edit are present in both libraries."
        )
    parser.add_argument(
        "command1",
        type=str,
        help="bean run command for the first batched run.",
    )
    parser.add_argument(
        "command2",
        type=str,
        help="bean run command for the second batched run.",
    )
    parser.add_argument(
        "raw_run_output1",
        type=str,
        help="bean run output .pkl path for the first batched run, which should be ran with --save-raw",
    )
    parser.add_argument(
        "output_path",
        type=str,
        help="Output path to save prior parameters.",
    )
    return parser
