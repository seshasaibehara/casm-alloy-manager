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
    predict = subparser.add_parser(
        "predict",
        help="Maps relaxed structures of a casm project to a set of parent crystal structures and predicts the best parent crystal structure for every configuration",
    )

    # List of configurations in ccasm query json format
    predict.add_argument(
        "--configurations",
        "-c",
        type=str,
        required=True,
        help="List of configurations in ccasm query json format",
    )

    # outfile name. If outfile name is *.html, results will be written out to html file
    # If outfile name is *.hdf, results will be written to hdf file
    predict.add_argument(
        "--outfile",
        "-o",
        type=str,
        required=True,
        help="Output file name (a pandas dataframe dumped as a hdf5/html) file",
    )

    predict.add_argument(
        "--configtype",
        nargs="?",
        type=str,
        default="relaxed",
        choices=["relaxed", "unrelaxed"],
        help="What to read. If relaxed will read properties.calc.json, if unrelaxed will read structure.json",
    )

    # ignored if reading unrelaxed structure.json files
    predict.add_argument(
        "--calctype",
        nargs="?",
        type=str,
        default="default",
        help="calctype from where to read the properties",
    )

    # TODO: Add option for user to give parent crystal structures
    # TODO: Rename top5
    predict.add_argument(
        "--parents",
        "-p",
        nargs="?",
        type=str,
        default="frequent",
        choices=["all", "frequent"],
        help="What parent structures to use",
    )

    args = parser.parse_args()

    # read configurations
    if args.command == "predict":

        with open(args.configurations, "r") as f:
            configs = json.load(f)

        config_names = [config["name"] for config in configs]

        if args.configtype == "relaxed":
            relaxed = True
        if args.configtype == "unrelaxed":
            relaxed = False
        # construct child structures to be used in mapping
        # get child properties paths
        child_paths = casmam.mapping.mapping.get_properties_json_paths(
            config_names, args.calctype, relaxed
        )
        mapping_results = (
            casmam.mapping.mapping.map_configurations_onto_parent_structures(
                child_paths, args.parents
            )
        )

        if ".html" in args.outfile:
            mapping_results.to_html(args.outfile)

        if ".hdf" in args.outfile:
            mapping_results.to_hdf(args.outfile, key="mapping_results")


if __name__ == "main":
    main()
