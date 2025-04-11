import json
import casmam
import warnings
import argparse
import pandas as pd

warnings.filterwarnings("ignore", category=pd.io.pytables.PerformanceWarning)


def main():
    parser = argparse.ArgumentParser("casm-alloy-manager")
    subparser = parser.add_subparsers(dest="command")

    # map command
    mapper = subparser.add_parser(
        "map",
        help="Maps relaxed structures of a casm project to a set of parent crystal structures and predicts the best parent crystal structure for every configuration",
    )

    mapper.add_argument(
        "--child",
        "-c",
        type=str,
        required=True,
        help="A file containing a list of where relaxed/child structures can be found",
    )

    mapper.add_argument(
        "--parents",
        "-p",
        nargs="?",
        type=str,
        default="common",
        choices=["all", "common"],
        help="What parent structures to use",
    )

    mapper.add_argument(
        "--settings",
        "-s",
        type=str,
        default=None,
        help="All settings to be passed onto the mapper",
    )

    # outfile name. If outfile name is *.html, results will be written out to html file
    # If outfile name is *.hdf, results will be written to hdf file
    mapper.add_argument(
        "--outfile",
        "-o",
        type=str,
        required=True,
        help="Output file name (a pandas dataframe dumped as a hdf5/html) file",
    )

    args = parser.parse_args()

    if args.command == "map":

        if args.settings is not None:
            with open(args.settings, "r") as f:
                settings = json.load(f)

        else:
            settings = {}

        # read child structure paths
        with open(args.child, "r") as f:
            child_structure_paths = f.read().splitlines()

        # read parent structure paths
        if not (args.parents == "all" or args.parents == "common"):
            with open(args.parents, "r") as f:
                parent_structure_paths = f.read().splitlines()
            mapping_results = casmam.map.map(
                child_structure_paths, parent_structure_paths, **settings
            )
        else:
            mapping_results = casmam.map.map(
                child_structure_paths, args.parents, **settings
            )

        if ".html" in args.outfile:
            mapping_results.to_html(args.outfile)

        if ".hdf" in args.outfile:
            mapping_results.to_hdf(args.outfile, key="mapping_results")


if __name__ == "main":
    main()
