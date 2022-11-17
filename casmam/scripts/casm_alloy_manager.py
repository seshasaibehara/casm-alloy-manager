import json
import casmam
import warnings
import argparse
import pandas as pd

warnings.filterwarnings("ignore", category=pd.io.pytables.PerformanceWarning)


def function(arg1):
    """TODO: Docstring for function.

    Parameters
    ----------
    arg1 : TODO

    Returns
    -------
    TODO

    """
    pass


def main():
    parser = argparse.ArgumentParser("casm-alloy-manager")
    subparser = parser.add_subparsers(dest="command")

    # map command
    mapper = subparser.add_parser(
        "map",
        help="Maps relaxed structures of a casm project to a set of parent crystal structures and predicts the best parent crystal structure for every configuration",
    )

    # List of configurations in ccasm query json format
    mapper.add_argument(
        "--configurations",
        "-c",
        type=str,
        required=True,
        help="List of configurations in ccasm query json format",
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

    mapper.add_argument(
        "--configtype",
        nargs="?",
        type=str,
        default="relaxed",
        choices=["relaxed", "unrelaxed"],
        help="What to read. If relaxed will read properties.calc.json, if unrelaxed will read structure.json",
    )

    mapper.add_argument(
        "--calctype",
        nargs="?",
        type=str,
        default="default",
        help="calctype from where to read the properties",
    )

    # TODO: Add option for user to give parent crystal structures
    mapper.add_argument(
        "--parents",
        "-p",
        nargs="?",
        type=str,
        default="common",
        choices=["all", "common"],
        help="What parent structures to use",
    )
    # TODO: Add input settings to mapping arguments
    # TODO: Add input settings to orgainizing mapping results

    # analyze command
    analyze = subparser.add_parser(
        "analyze",
        help="Analyzes the given mapping results and finds best map to each of the configurations",
    )

    # TODO: what to do if it's html
    analyze.add_argument(
        "--infile", "-i", type=str, required=True, help="Mapping results as hdf5 file"
    )

    # TODO: need a html argument?
    analyze.add_argument(
        "--outfile", "-o", type=str, required=True, help="Output file name"
    )

    args = parser.parse_args()

    # read configurations
    if args.command == "map":

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

    if args.command == "analyze":
        mapping_results = pd.read_hdf(args.infile)
        best_maps = casmam.mapping.mapping.analyze_mapping_data(mapping_results)

        best_maps.to_hdf(args.outfile, key="best_maps")


if __name__ == "main":
    main()
