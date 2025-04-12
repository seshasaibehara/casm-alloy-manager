import os
import copy
import json
import casmam
import numpy as np
import pandas as pd
import libcasm.xtal
from tqdm import tqdm
import importlib.resources
import libcasm.mapping.info
import libcasm.mapping.methods


def common_parent_crystal_structures_with_paths() -> (
    tuple[list[libcasm.xtal.Prim], list[str]]
):
    """Construct libcasm.Prim objects for common parent crystal structures
    found in xtallib.
    Also returns the basepaths of the filenames of these crystal structures.

    Common parent crystal structures include BCC, FCC, HCP, Omega, SC, DHCP

    Returns
    -------
    tuple[list[libcasm.xtal.Prim], list[str]]
        A list of casm ``Prim`` which can be used as an input for mapping
        along with a list of basepath names for crystal structures

    """
    prims_with_names = [
        (
            libcasm.xtal.Prim.from_poscar(str(file)),
            str(os.path.basename(file)),
        )
        for file in importlib.resources.files("casmam.xtallib.common").iterdir()
        if ".vasp" in str(file)
    ]

    prims = [entry[0] for entry in prims_with_names]
    paths = [entry[1] for entry in prims_with_names]

    return prims, paths


def all_parent_crystal_structures_with_paths() -> (
    tuple[list[libcasm.xtal.Prim], list[str]]
):
    """Makes casm ``Prim`` for all the parent crystals from Sanjeev's database
    which are included in xtallib.
    Also returns the basepaths of the filenames of these crystal structures.

    Returns
    -------
    tuple[list[libcasm.xtal.Prim], list[str]]
     A list of casm ``Prim`` which can be used as an input for mapping
        along with a list of basepath names for crystal structures

    """
    prims_with_names = common_parent_crystal_structures_with_paths()
    prims_with_names += [
        (
            libcasm.xtal.Prim.from_poscar(str(file)),
            str(file),
        )
        for file in importlib.resources.files("casmam.xtallib").iterdir()
        if ".vasp" in str(file)
    ]

    prims = [entry[0] for entry in prims_with_names]
    paths = [entry[1] for entry in prims_with_names]

    return prims, paths


def max_vol(parent: libcasm.xtal.Prim, child: libcasm.xtal.Structure) -> int | None:
    """Returns if number of atoms in child is divisible by number of atoms
    in parent structure. Mapping algorithm by default finds supercells of parent,
    but not child. If number of atoms in child is not divisible by parent,
    mapping algorithm cannot find supercells of the parent

    # TODO: Add support for vacancies

    Parameters
    ----------
    parent : libcasm.xtal.Prim
        casm ``Prim`` of a parent crystal structure
    child : libcasm.xtal.Structure
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
    child_structures: list[libcasm.xtal.Structure], masking_atom_type: str = "A"
) -> list[libcasm.xtal.Structure]:
    """Given a list of child structures, change the atom types at all sites
    in all of the structures to ``masking_atom_type`` which is "A"

    Parameters
    ----------
    child_structures : List[libcasm.xtal.Structure]
        List of child structures to change the atom type
    masking_atom_type : str, optional
        Replace atom type at each site to ``masking_atom_type``
        By default it will be replace to ``A``

    Returns
    -------
    List[libcasm.xtal.Structure]
        List of child structures with atom type at each site to ``masking_atom_type``

    """
    return [
        casmam.mask.change_atom_types_in_casm_structure(
            child_structure, [masking_atom_type] * len(child_structure.atom_type())
        )
        for child_structure in child_structures
    ]


def mask_parent_structure_atom_dofs(
    parent_structures: list[libcasm.xtal.Prim], masking_atom_dof_type: str = "A"
) -> list[libcasm.xtal.Structure]:
    """Given a list of child structures, change the atom types at all sites
    in all of the structures to ``masking_atom_type`` which is "A"

    Parameters
    ----------
    child_structures : List[libcasm.xtal.Structure]
        List of child structures to change the atom type
    masking_atom_type : str, optional
        Replace atom type at each site to ``masking_atom_type``
        By default it will be replace to ``A``

    Returns
    -------
    List[libcasm.xtal.Structure]
        List of child structures with atom type at each site to ``masking_atom_type``

    """
    return [
        casmam.mask.change_dofs_in_casm_prim(
            parent_structure,
            [[masking_atom_dof_type]] * len(parent_structure.occ_dof()),
        )
        for parent_structure in parent_structures
    ]


def make_child_structures(child_paths: list[str]) -> list[libcasm.xtal.Structure]:
    """Given a list of child paths, if the file type is json,
    it assumes it's of proeprties.calc.json/structure.json type
    and constructs a casm ``Structure``

    Parameters
    ----------
    child_paths : list[str]
        List of child paths

    Returns
    -------
    list[libcasm.xtal.Structure]
        List of casm ``Structure`` objects

    """
    child_structures = []
    for child_path in child_paths:
        if ".json" in os.path.basename(child_path):
            with open(child_path, "r") as f:
                properties_dictionary = json.load(f)
            child_structures.append(
                libcasm.xtal.Structure.from_dict(properties_dictionary)
            )
        else:
            child_structures.append(libcasm.xtal.Structure.from_poscar(child_path))

    return child_structures


def make_parent_crystal_structures_and_paths(
    parent_paths: str | list[str],
) -> list[libcasm.xtal.Prim, str]:
    """Make libcasm.xtal.Prim from parent crystal structure library

    Parameters
    ----------
    parent_paths : str | list[str]
        Valid options are "common", "all" or a list of paths to POSCARs

    Returns
    -------
    list[tuple[libcasm.xtal.Prim, str]]

    """
    if isinstance(parent_paths, str):
        if parent_paths == "common":
            (
                parent_structures,
                parent_names,
            ) = common_parent_crystal_structures_with_paths()

        elif parent_paths == "all":
            parent_structures, parent_names = all_parent_crystal_structures_with_paths()

        else:
            raise RuntimeError("Invalid library (" + parent_paths + ") of structures")

    if isinstance(parent_paths, list):
        parent_structures = [
            libcasm.xtal.Prim.from_poscar(parent_path) for parent_path in parent_paths
        ]
        parent_names = copy.deepcopy(parent_paths)

    return parent_structures, parent_names


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
        "max_cost": 0.1,
        "lattice_cost_method": "symmetry_breaking_strain_cost",
        "atom_cost_method": "symmetry_breaking_disp_cost",
    }


def default_child_structure_options() -> dict:
    """Returns a dictionary of default child structure
    options

    Returns
    -------
    dict

    """
    return {"mask_occupants": True}


def default_parent_structure_options() -> dict:
    """Returns a dictionary of default child structure
    options

    Returns
    -------
    dict

    """
    return {"mask_occupants": True}


def map(
    child_paths: list[str], parent_paths: str | list[str], **kwargs
) -> pd.DataFrame:
    """Top-level function that constructs child structures,
    parent structures and maps them onto each other

    Parameters
    ----------
    child_paths : list[str]
        List of paths to child structures
    parent_paths : str | list[str]
        "common", "all" from xtallib or list of paths to parent structures
    **kwargs : TODO

    Returns
    -------
    pd.DataFrame

    """
    # sanitize kwargs
    mapping_options = default_mapping_options()
    child_structure_options = default_child_structure_options()
    parent_structure_options = default_parent_structure_options()

    for key, value in kwargs.items():
        if key == "mapping_options":
            for mapping_key, mapping_option in value.items():
                mapping_options[mapping_key] = mapping_option

        if key == "child_structure_options":
            for child_struc_key, child_option in value.items():
                child_structure_options[child_struc_key] = child_option

        if key == "parent_structure_options":
            for parent_struc_key, parent_option in value.items():
                parent_structure_options[parent_struc_key] = parent_option

    # make prims from parent crystal paths
    parent_structures, parent_names = make_parent_crystal_structures_and_paths(
        parent_paths
    )
    child_structures = make_child_structures(child_paths)

    if child_structure_options["mask_occupants"] is True:
        child_structures = mask_child_structure_atom_types(child_structures)

    if parent_structure_options["mask_occupants"] is True:
        parent_structures = mask_parent_structure_atom_dofs(parent_structures)

    mapping_results = map_child_structures_onto_parent_structures(
        parent_structures,
        child_structures,
        mapping_options,
    )
    mapping_results = organize_mapping_results(
        mapping_results, child_paths, parent_names
    )

    return mapping_results


def map_child_structures_onto_parent_structures(
    parent_structures: list[libcasm.xtal.Prim],
    child_structures: list[libcasm.xtal.Structure],
    mapping_options: dict,
) -> list[list[libcasm.mapping.info.StructureMappingResults]]:
    """Cycle through child crystal structures and map each of them
    onto the parent structures. Assumes ``parent_structures`` and
    ``child_structures`` have the desired atom types. If not, use
    helper :func:``casmam.mapping.mask_child_structure_atom_types``
    function. Need ``parent_paths`` and ``child_paths`` to keep track of
    these files in ``MappingResult``.

    Parameters
    ----------
    parent_structures : List[libcasm.xtal.Prim]
        List of parent crystal structures as casm ``Prim``
    child_structures : List[libcasm.xtal.Structure]
        List of child crystal structures as casm ``Structure``

    Returns
    -------
    list[list[libcasm.mapping.info.StructureMappingResults]]

    """
    # sanitize mapping_options
    use_child_symmetry = False
    if mapping_options["use_child_symmetry"] is True:
        use_child_symmetry = True

    use_parent_symmetry = False
    if mapping_options["use_parent_symmetry"] is True:
        use_parent_symmetry = True

    mapping_options.pop("use_child_symmetry")
    mapping_options.pop("use_parent_symmetry")

    mapping_results = []
    for child_structure in tqdm(child_structures, desc="Mapping child structure: "):
        # Make child factor group
        if use_child_symmetry is True:
            child_fg = libcasm.xtal.make_factor_group(child_structure)
        else:
            child_fg = []

        mapping_results_for_one_child = []
        for parent_structure in tqdm(
            parent_structures, desc="Mapping one child onto parent: "
        ):

            # make parent factor group
            if use_parent_symmetry is True:
                parent_fg = libcasm.xtal.make_factor_group(parent_structure)
            else:
                parent_fg = []

            # map child onto parent
            max_volume = max_vol(parent_structure, child_structure)

            if max_volume is not None:
                results = libcasm.mapping.methods.map_structures(
                    prim=parent_structure,
                    structure=child_structure,
                    max_vol=max_volume,
                    prim_factor_group=parent_fg,
                    structure_factor_group=child_fg,
                    **mapping_options,
                )
            else:
                results = []

            mapping_results_for_one_child.append(results)

        mapping_results.append(mapping_results_for_one_child)

    return mapping_results


def organize_mapping_results(
    mapping_results, child_structure_names, parent_structure_names
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
    mapping_table_entries = []
    for map_results_of_one_child in mapping_results:
        table_row_entry = []
        for map_result_of_one_parent in map_results_of_one_child:
            if len(map_result_of_one_parent) == 0:
                atomic_cost = np.nan
                lattice_cost = np.nan
                total_cost = np.nan
                mapping_result_dict = {}

            else:
                atomic_cost = map_result_of_one_parent[0].atom_cost()
                lattice_cost = map_result_of_one_parent[0].lattice_cost()
                total_cost = map_result_of_one_parent[0].total_cost()
                mapping_result_dict = map_result_of_one_parent[0].to_dict()

            table_row_entry.extend(
                [atomic_cost, lattice_cost, total_cost, mapping_result_dict]
            )

        mapping_table_entries.append(table_row_entry)

    mapping_results_table = pd.DataFrame(
        mapping_table_entries,
        index=child_structure_names,
        columns=pd.MultiIndex.from_product(
            [
                parent_structure_names,
                ["atom_cost", "lattice_cost", "total_cost", "mapping_results"],
            ]
        ),
    )

    return mapping_results_table


# def find_best_map_and_flag_conflicts(
#    mapping_results: list[MappingResult], tol: float = 1e-4
# ) -> tuple[MappingResult, list[MappingResult] | None]:
#    """Given a list of mapping results, finds
#    the best mapping by finding map with least ``total_cost``
#    If there are multiple maps with ``total_cost`` close to
#    the best map's total cost, they will also be returned
#
#    Parameters
#    ----------
#    mapping_results : list[MappingResult]
#        list of MappingResults
#
#    Returns
#    -------
#    tuple[MappingResult, list[MappingResult]]
#        A tuple of best map with a list of maps within a given ``tol``
#        of the best map
#    """
#
#    total_costs = np.array(
#        [mapping_result.total_cost for mapping_result in mapping_results]
#    )
#
#    if all([mapping_result.is_dummy() for mapping_result in mapping_results]):
#        empty_mapping_result = MappingResult()
#        empty_mapping_result.child_path = mapping_results[0].child_path
#        return empty_mapping_result, None
#
#    best_map_index = np.nanargmin(total_costs)
#    best_map = mapping_results[best_map_index]
#    best_total_cost = best_map.total_cost
#
#    conflicting_indices = [
#        index
#        for index, total_cost in enumerate(total_costs)
#        if np.isclose(total_cost, best_total_cost, tol, tol) and index != best_map_index
#    ]
#
#    conflicting_maps = [mapping_results[index] for index in conflicting_indices]
#    if len(conflicting_maps) == 0:
#        conflicting_maps = None
#
#    return best_map, conflicting_maps
#
#
# def analyze_mapping_data(
#    mapping_data: pd.DataFrame, tol: float = 1e-4, **kwargs
# ) -> pd.DataFrame:
#    """TODO: Docstring for analyze_mapping_data.
#
#    Parameters
#    ----------
#    mapping_data : TODO
#
#    Returns
#    -------
#    TODO
#
#    """
#
#    keys_with_mapping_results = [
#        key for key in mapping_data if "mapping_results" in key
#    ]
#
#    if len(keys_with_mapping_results) == 0:
#        raise RuntimeError("Provided DataFrame does not contain MappingResult objects")
#
#    mapping_results_of_all_configs = mapping_data.loc[:, keys_with_mapping_results]
#
#    config_names = []
#    table_entries = []
#    for config_name, config_data in mapping_results_of_all_configs.iterrows():
#        config_mapping_results = [
#            config_mapping_result
#            for _, config_mapping_result in config_data.iteritems()
#        ]
#        best_config_map, conflicting_maps = find_best_map_and_flag_conflicts(
#            config_mapping_results, tol
#        )
#
#        table_entries.append(
#            [best_config_map.parent_path, best_config_map, conflicting_maps]
#        )
#        config_names.append(config_name)
#
#    best_map_table = pd.DataFrame(
#        table_entries,
#        index=config_names,
#        columns=[
#            "Best parent map name",
#            "Best parent mapping object",
#            "Conflicting maps",
#        ],
#    )
#
#    return best_map_table
