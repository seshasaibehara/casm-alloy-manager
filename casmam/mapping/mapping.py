import os
import json
import casm.xtal
import numpy as np
import pandas as pd
import importlib.resources
import casmam.xtal.xtal as casmamxtal
import casm.mapping.info as mapperinfo
import casm.mapping.methods as mappermethods


class MappingResult:

    """An object containing all the mapping result
    info and which can be dumped as a picke

    """

    def __init__(self):
        """TODO: to be defined."""
        self.parent_path = "not available"
        self.child_path = "not available"

        # mapping costs
        self.atomic_cost = np.nan
        self.lattice_cost = np.nan
        self.total_cost = np.nan

        # lattice mapping results
        # TODO: Explain parameters
        self.deformation_gradient = np.nan
        self.transformation_matrix_to_super = np.nan
        self.reorientation = np.nan
        self.isometry = np.nan
        self.left_stretch = np.nan

        # atom mapping results
        # TODO: Explain parameters
        self.displacement = np.nan
        self.permutation = np.nan
        self.translation = np.nan

    def is_dummy(self):
        """Returns if MappingResult is a dummy

        Returns
        -------
        TODO

        """
        return np.isnan(self.total_cost)

    def __str__(self):
        """TODO: Docstring for __str__.

        Returns
        -------
        TODO

        """
        string = "Total mapping cost is: "
        string += str(self.total_cost)
        return string

    @classmethod
    def from_casm_mapping_result(
        cls,
        casm_mapping_result: tuple[
            mapperinfo.StructureMappingCost, mapperinfo.StructureMapping
        ],
        parent_path="Not available",
        child_path="Not available",
    ):
        """TODO: Docstring for from_casm_mapping_result.

        Parameters
        ----------
        casm_mapping_result: TODO

        Returns
        -------
        TODO

        """
        result = cls()

        result.parent_path = parent_path
        result.child_path = child_path

        # populate mapping costs
        result.atomic_cost = casm_mapping_result[0].atom_cost()
        result.lattice_cost = casm_mapping_result[0].lattice_cost()
        result.total_cost = casm_mapping_result[0].total_cost()

        # populate lattice mapping attributes
        lattice_mapping = casm_mapping_result[1].lattice_mapping()
        result.deformation_gradient = lattice_mapping.deformation_gradient()
        result.transformation_matrix_to_super = (
            lattice_mapping.transformation_matrix_to_super()
        )
        result.reorientation = lattice_mapping.reorientation()
        result.isometry = lattice_mapping.isometry()
        result.left_stretch = lattice_mapping.left_stretch()

        # populate atom mapping attributes
        atom_mapping = casm_mapping_result[1].atom_mapping()
        result.displacement = atom_mapping.displacement()
        result.permutation = atom_mapping.permutation()
        result.translation = atom_mapping.translation()

        return result


def get_casm_root_dir() -> str:
    """Find casm root directory and return the same.
    Directory where .casm lies

    Returns
    -------
    os.path
        Root of casm project

    Raises
    ------
    RuntimeError
        If not in a casm project

    """
    cwd = os.getcwd()
    while True:
        if os.path.isdir(os.path.join(cwd, ".casm")):
            return cwd
        elif cwd == os.path.dirname(cwd):
            # escape while
            raise RuntimeError("Not in a casm project")
        else:
            # crawl up dir before
            cwd = os.path.dirname(cwd)


def default_parent_crystal_structures_with_paths() -> tuple[
    list[casm.xtal.Prim], list[str]
]:
    """Makes casm ``Prim`` for BCC, FCC, HCP, Omega, SC, DHCP

    Returns
    -------
    List[casm.xtal.Prim]
        A list of casm ``Prim`` which can be used as an input for mapping

    """
    prims_with_names = [
        (
            casm.xtal.Prim.from_poscar(str(file)),
            str(file),
        )
        for file in importlib.resources.files("casmam.xtallib.frequent").iterdir()
        if ".vasp" in str(file)
    ]

    prims = [entry[0] for entry in prims_with_names]
    paths = [entry[1] for entry in prims_with_names]

    return prims, paths


def all_parent_crystal_structures_with_paths() -> tuple[
    list[casm.xtal.Prim], list[str]
]:
    """Makes casm ``Prim`` for all the parent crystals from
    Sanjeev's database

    Returns
    -------
    List[casm.xtal.Prim]

    """
    prims_with_names = default_parent_crystal_structures_with_paths()
    prims_with_names += [
        (
            casm.xtal.Prim.from_poscar(str(file)),
            str(file),
        )
        for file in importlib.resources.files("casmam.xtallib").iterdir()
        if ".vasp" in str(file)
    ]

    prims = [entry[0] for entry in prims_with_names]
    paths = [entry[2] for entry in prims_with_names]

    return prims, paths


def max_vol(parent: casm.xtal.Prim, child: casm.xtal.Structure) -> int | None:
    """Returns if number of atoms in child is divisible by number of atoms
    in parent structure. Mapping algorithm by default finds supercells of parent,
    but not child. If number of atoms in child is not divisible by parent,
    mapping algorithm cannot find supercells of the parent

    Parameters
    ----------
    parent : casm.xtal.Prim
        casm ``Prim`` of a parent crystal structure
    child : casm.xtal.Structure
        casm ``Structure`` of a child crystal structure

    Returns
    -------
    bool
        Returns ``True`` if number of atoms in child is divisible
        by number of atoms in parent. Else returns ``False``. Returns
        ``None`` if it's not divisible

    """
    max_vol = len(child.atom_type()) / len(parent.occ_dof())
    is_divisible = len(child.atom_type()) % len(parent.occ_dof())
    if is_divisible == 0:
        return int(max_vol)

    return None


def mask_child_structure_atom_types(
    child_structures: list[casm.xtal.Structure], masking_atom_type: str = "A"
) -> list[casm.xtal.Structure]:
    """Given a list of child structures, change the atom types at all sites
    in all of the structures to ``masking_atom_type`` which is "A"

    Parameters
    ----------
    child_structures : List[casm.xtal.Structure]
        List of child structures to change the atom type
    masking_atom_type : str, optional
        Replace atom type at each site to ``masking_atom_type``
        By default it will be replace to ``A``

    Returns
    -------
    List[casm.xtal.Structure]
        List of child structures with atom type at each site to ``masking_atom_type``

    """
    return [
        casmamxtal.change_atom_types_in_casm_structure(
            child_structure, [masking_atom_type] * len(child_structure.atom_type())
        )
        for child_structure in child_structures
    ]


def get_child_structures(child_paths: list[str]) -> list[casm.xtal.Structure]:
    """Given a list of child paths, if the file type is json,
    it assumes it's of proeprties.calc.json/structure.json type
    and constructs a casm ``Structure``

    Parameters
    ----------
    child_paths : list[str]
        List of child paths

    Returns
    -------
    list[casm.xtal.Structure]
        List of casm ``Structure`` objects

    """
    child_structures = []
    for child_path in child_paths:
        if ".json" in os.path.basename(child_path):
            with open(child_path, "r") as f:
                properties_dictionary = json.load(f)

            child_structures.append(
                casmamxtal.casm_structure_from_properties_json(properties_dictionary)
            )
        else:
            raise NotImplementedError(
                "Reading child structures from POSCAR not implemented yet"
            )

    return child_structures


def get_properties_json_paths(
    config_names: list[str], calctype="default", relaxed=True
):
    """TODO: Docstring for get_properties_json_paths.

    Parameters
    ----------
    config_names : TODO

    Returns
    -------
    TODO

    """
    casm_root_dir = get_casm_root_dir()
    if relaxed:
        return [
            os.path.join(
                casm_root_dir,
                "training_data",
                config_name,
                "calctype." + calctype,
                "properties.calc.json",
            )
            for config_name in config_names
        ]

    else:
        return [
            os.path.join(casm_root_dir, "training_data", config_name, "structure.json")
            for config_name in config_names
        ]


# TODO: Currently only works if child_paths is a list of .json files
# TODO: If child_paths is poscar types it doesn't work, needs implementing Structure.from_poscar
# TODO: Need to add shorten parent paths argmument
def map_configurations_onto_parent_structures(
    child_paths: list[str],
    parent_paths: str | list[str],
    quiet=False,
    **kwargs,
):
    """Top-level function that constructs child structures,
    parent structures and maps them onto each other

    Parameters
    ----------
    child_paths : TODO
    parent_paths : TODO
    **kwargs : TODO

    Returns
    -------
    TODO

    """
    # sanitize kwargs
    if isinstance(parent_paths, str):
        if parent_paths == "frequent":
            (
                parent_structures,
                parent_paths,
            ) = default_parent_crystal_structures_with_paths()

        elif parent_paths == "all":
            parent_structures, parent_paths = all_parent_crystal_structures_with_paths()

        else:
            raise RuntimeError("Invalid library (" + parent_paths + ") of structures")

    if isinstance(parent_paths, list):
        parent_structures = [
            casm.xtal.Prim.from_poscar(parent_path) for parent_path in parent_paths
        ]

    # parent_paths = [os.path.basename(path) for path in parent_paths]

    child_structures = get_child_structures(child_paths)
    masked_child_structures = mask_child_structure_atom_types(child_structures)

    mapping_results = map_child_structures_onto_parent_structures(
        parent_structures, masked_child_structures, parent_paths, child_paths, quiet
    )
    mapping_results = organize_mapping_results(mapping_results)

    return mapping_results


def default_mapping_options() -> dict:
    """Returns a dictionary of default mapping options
    used

    Returns
    -------
    dict
        Dictionary of default mapping options

    """
    return {
        "use_parent_symmetry": True,
        "use_child_symmetry": True,
        "min_cost": -0.0001,
        "max_cost": 0.1,
        "strain_cost_method": "symmetry_breaking_strain_cost",
        "atom_cost_method": "symmetry_breaking_atom_cost",
    }


def map_child_structures_onto_parent_structures(
    parent_structures: list[casm.xtal.Prim],
    child_structures: list[casm.xtal.Structure],
    parent_paths: list[str] = None,
    child_paths: list[str] = None,
    quiet=True,
    **kwargs,
) -> list[list[list[MappingResult]]]:
    """Cycle through child crystal structures and map each of them
    onto the parent structures. Assumes ``parent_structures`` and
    ``child_structures`` have the desired atom types. If not, use
    helper :func:``casmam.mapping.mask_child_structure_atom_types``
    function. Need ``parent_paths`` and ``child_paths`` to keep track of
    these files in ``MappingResult``.

    Parameters
    ----------
    parent_structures : List[casm.xtal.Prim]
        List of parent crystal structures as casm ``Prim``
    child_structures : List[casm.xtal.Structure]
        List of child crystal structures as casm ``Structure``
    **kwargs : TODO

    Returns
    -------
    list[list[MappingResult]]

    """
    if parent_paths is None:
        parent_paths = ["not available"] * len(parent_structures)

    if child_paths is None:
        child_paths = ["not available"] * len(child_structures)

    # TODO: Sanitize args and kwargs. Think about what to expose to the user
    mapping_results = []

    for child_structure, child_path in zip(child_structures, child_paths):
        child_fg = casm.xtal.make_structure_factor_group(child_structure)

        mapping_results_for_one_child = []
        for parent_structure, parent_path in zip(parent_structures, parent_paths):
            # make child structure factor group
            parent_fg = casm.xtal.make_prim_factor_group(parent_structure)
            # map child onto parent
            max_volume = max_vol(parent_structure, child_structure)
            if max_volume is not None:
                results = mappermethods.map_structures(
                    parent_structure,
                    child_structure,
                    max_vol=max_volume,
                    prim_factor_group=parent_fg,
                    structure_factor_group=child_fg,
                    strain_cost_method="symmetry_breaking_strain_cost",
                    atom_cost_method="symmetry_breaking_atom_cost",
                    min_cost=-0.0001,
                    max_cost=0.1,
                )
            else:
                results = []

            if len(results) == 0:
                empty_mapping_result = MappingResult()
                empty_mapping_result.child_path = child_path
                empty_mapping_result.parent_path = parent_path
                casmam_results = [empty_mapping_result]

            else:
                casmam_results = [
                    MappingResult.from_casm_mapping_result(
                        result, parent_path, child_path
                    )
                    for result in results
                ]

            if not quiet:
                print("Finished mapping " + child_path + " to " + parent_path + "...")

            # convert the results to casmam MappingResult object which is
            # pickle dumpable
            mapping_results_for_one_child.append(casmam_results)

        mapping_results.append(mapping_results_for_one_child)

    return mapping_results


def get_casm_config_name_from_child_path(child_path: str) -> str:
    """This assumes that the child path is in casm directory
    style like "*/training_data/SCEL.../*/structure.json"
    or "*/training_data/SCEL.../*/calctype.*/properties.json"
    and only return SCEL.../* part

    Returns
    -------
    str
        config name

    """
    if os.path.basename(child_path) == "structure.json":
        i = 0
        dirname = child_path
        while i <= 2:
            dirname = os.path.dirname(dirname)
            i += 1

        config_name = child_path.replace(dirname + "/", "").replace(
            "/structure.json", ""
        )

    elif os.path.basename(child_path) == "properties.calc.json":
        i = 0
        dirname = child_path
        while i <= 3:
            dirname = os.path.dirname(dirname)
            if i == 0:
                child_path = child_path.replace(os.path.basename(dirname) + "/", "")
            i += 1

        config_name = child_path.replace(dirname + "/", "").replace(
            "/properties.calc.json", ""
        )

    else:
        raise RuntimeError(
            "Should be casm directory styles with file name either structure.json or properties.calc.json"
        )

    return config_name


def organize_mapping_results(
    mapping_results: list[list[list[MappingResult]]],
    **kwargs,
) -> pd.DataFrame:
    """Organize ``mapping_results`` into a pandas ``DataFrame`` with each
    row containing ``atomic_cost``, ``lattice_cost`` and ``total_cost`` of mapping
    one child structure onto all the parent crystal structures. By default, will only
    include the first best map to each parent crystal structure. If
    ``include_all_valid_mappings`` is ``True``, every valid map to each
    parent crystal structure will be included. If ``include_mapping_objects``
    is ``True``, DataFrame will also include mapping result objects in
    addition to ``atomic_cost`` , ``lattice_cost`` and ``total_cost``

    Parameters
    ----------
    mapping_results : list[list[list[MappingResult]]]
        Mapping results returned by :func:``map_child_structures_onto_given_parent_structures``

    Returns
    -------
    pd.DataFrame
        Mapping results organized into a pandas DataFrame

    """

    # sanitize kwargs and get these
    shorten_child_names = "casm_style"
    shorten_parent_names = "normal"

    # shorten child names
    if shorten_child_names == "casm_style":
        child_structure_names = [
            get_casm_config_name_from_child_path(
                map_results_of_one_child[0][0].child_path
            )
            for map_results_of_one_child in mapping_results
        ]
    elif shorten_child_names == "normal":
        child_structure_names = [
            os.path.basename(map_results_of_one_child[0][0].child_path)
            for map_results_of_one_child in mapping_results
        ]

    else:
        child_structure_names = [
            map_results_of_one_child[0][0].child_path
            for map_results_of_one_child in mapping_results
        ]

    # shorten parent names
    if shorten_parent_names == "normal":
        parent_structure_names = [
            os.path.basename(result[0].parent_path) for result in mapping_results[0]
        ]
    else:
        parent_structure_names = [
            result[0].parent_path for result in mapping_results[0]
        ]

    mapping_table_entries = []
    for map_results_of_one_child in mapping_results:
        table_row_entry = []
        for map_result_of_one_parent in map_results_of_one_child:
            atomic_cost = map_result_of_one_parent[0].atomic_cost
            lattice_cost = map_result_of_one_parent[0].lattice_cost
            total_cost = map_result_of_one_parent[0].total_cost

            table_row_entry.extend([atomic_cost, lattice_cost, total_cost])
            table_row_entry.append(map_result_of_one_parent[0])

        mapping_table_entries.append(table_row_entry)

    mapping_results_table = pd.DataFrame(
        mapping_table_entries,
        index=child_structure_names,
        columns=pd.MultiIndex.from_product(
            [
                parent_structure_names,
                # ["atomic_cost", "lattice_cost", "total_cost"]
                ["atomic_cost", "lattice_cost", "total_cost", "mapping_results"],
            ]
        ),
    )

    return mapping_results_table


def find_best_map_and_flag_conflicts(
    mapping_results: list[MappingResult], tol: float = 1e-4
) -> tuple[MappingResult, list[MappingResult] | None]:
    """Given a list of mapping results, finds
    the best mapping by finding map with least ``total_cost``
    If there are multiple maps with ``total_cost`` close to
    the best map's total cost, they will also be returned

    Parameters
    ----------
    mapping_results : list[MappingResult]
        list of MappingResults

    Returns
    -------
    tuple[MappingResult, list[MappingResult]]
        A tuple of best map with a list of maps within a given ``tol``
        of the best map
    """

    total_costs = np.array(
        [mapping_result.total_cost for mapping_result in mapping_results]
    )

    best_map_index = np.nanargmin(total_costs)
    best_map = mapping_results[best_map_index]
    best_total_cost = best_map.total_cost

    conflicting_indices = [
        index
        for index, total_cost in enumerate(total_costs)
        if np.isclose(total_cost, best_total_cost, tol, tol) and index != best_map_index
    ]

    conflicting_maps = [mapping_results[index] for index in conflicting_indices]
    if len(conflicting_maps) == 0:
        conflicting_maps = None

    return best_map, conflicting_maps


def analyze_mapping_data(
    mapping_data: pd.DataFrame, tol: float = 1e-4, **kwargs
) -> pd.DataFrame:
    """TODO: Docstring for analyze_mapping_data.

    Parameters
    ----------
    mapping_data : TODO

    Returns
    -------
    TODO

    """

    keys_with_mapping_results = [
        key for key in mapping_data if "mapping_results" in key
    ]

    if len(keys_with_mapping_results) == 0:
        raise RuntimeError("Provided DataFrame does not contain MappingResult objects")

    mapping_results_of_all_configs = mapping_data.loc[:, keys_with_mapping_results]

    config_names = []
    table_entries = []
    for config_name, config_data in mapping_results_of_all_configs.iterrows():
        config_mapping_results = [
            config_mapping_result
            for _, config_mapping_result in config_data.iteritems()
        ]
        best_config_map, conflicting_maps = find_best_map_and_flag_conflicts(
            config_mapping_results, tol
        )

        table_entries.append(
            [best_config_map.parent_path, best_config_map, conflicting_maps]
        )
        config_names.append(config_name)

    best_map_table = pd.DataFrame(
        table_entries,
        index=config_names,
        columns=[
            "Best parent map name",
            "Best parent mapping object",
            "Conflicting maps",
        ],
    )

    return best_map_table
