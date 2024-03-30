import argparse
from bean.mapping.utils import get_input_parser_count as get_count_parser
from bean.mapping.utils import get_input_parser as get_count_samples_parser
from bean.plotting.utils import parse_args as get_profile_parser
from bean.qc.utils import parse_args as get_qc_parser
from bean.annotate.utils import parse_args as get_filter_parser
from bean.model.run import parse_args as get_run_parser
from bean.framework.read_from_csvs import get_input_parser as get_create_screen_parser
from bean.cli.count import main as count
from bean.cli.count_samples import main as count_samples
from bean.cli.profile import main as profile
from bean.cli.qc import main as qc
from bean.cli.filter import main as filter
from bean.cli.run import main as run
from bean.cli.create_screen import main as create_screen


def get_parser():
    parser = argparse.ArgumentParser(prog="bean")
    subparsers = parser.add_subparsers(help="Subcommands", dest="subcommand")
    count_parser = subparsers.add_parser("count", help="count")
    count_parser = get_count_parser(count_parser)
    count_samples_parser = subparsers.add_parser("count-samples", help="count samples")
    count_samples_parser = get_count_samples_parser(count_samples_parser)
    profile_parser = subparsers.add_parser("profile", help="profile")
    profile_parser = get_profile_parser(profile_parser)
    qc_parser = subparsers.add_parser("qc", help="qc")
    qc_parser = get_qc_parser(qc_parser)
    filter_parser = subparsers.add_parser("filter", help="filter")
    filter_parser = get_filter_parser(filter_parser)
    run_parser = subparsers.add_parser("run", help="run")
    run_parser = get_run_parser(run_parser)
    create_screen_parser = subparsers.add_parser("create_screen", help="create")
    create_screen_parser = get_create_screen_parser(create_screen_parser)
    return parser


global_parser = None


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()
    if args.subcommand == "count":
        count(args)
    elif args.subcommand == "count-samples":
        count_samples(args)
    elif args.subcommand == "profile":
        profile(args)
    elif args.subcommand == "qc":
        qc(args)
    elif args.subcommand == "filter":
        filter(args)
    elif args.subcommand == "run":
        run(args)
    elif args.subcommand == "create-screen":
        create_screen(args)
    else:
        parser.print_help()
