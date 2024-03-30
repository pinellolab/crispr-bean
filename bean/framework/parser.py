import argparse


def get_input_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description="bean-create-screen parameters",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    parser.add_argument(
        "gRNA_info_table_csv",
        type=str,
        help="Path to gRNA info .csv table, with index at first column and column names at the first row.",
    )
    parser.add_argument(
        "samples_info_table_csv",
        type=str,
        help="Path to sample info .csv table, with index at first column and column names at the first row.",
    )
    parser.add_argument(
        "gRNA_counts_table_csv",
        type=str,
        help="Path to gRNA counts .csv table, with index at first column and column names at the first row.",
    )
    parser.add_argument(
        "-e",
        "--edits",
        type=str,
        help="Path to edit counts .csv table, with index at first column and column names at the first row.",
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        help="Output file prefix (output will be saved as `output_prefix.h5ad`). If not provided, `gRNA_counts_table_csv` file prefix is used.",
        default=None,
    )
    return parser
