import casm.xtal
import numpy as np


def get_structure_info_from_properties_json(
    properties_json: dict,
) -> tuple[casm.xtal.Lattice, np.ndarray, list[str]]:
    """Read casm properties.calc.json/structure.json dictionary and return
    relavent information to construct a structure object

    Parameters
    ----------
    properties_json : dict
        properties.calc.json dictionary

    Returns
    -------
    Tuple[casm.xtal.Lattice, np.ndarray, List[str]]
        Returns casm Lattice, fractional coordinates as a :math:`3 \\times \\mathbf{N}` matrix,
        atom types at each site

    """
    # read lattice
    lattice_row_vector_matrix = np.array(properties_json["lattice_vectors"])
    casm_lattice = casm.xtal.Lattice(np.transpose(lattice_row_vector_matrix))

    # read coords
    frac_coords = np.array(properties_json["atom_coords"]).transpose()
    # convert coord to frac if they are cartesian
    if properties_json["coordinate_mode"] == "Cartesian":
        frac_coords = casm.xtal.cartesian_to_fractional(casm_lattice, frac_coords)

    atom_types = properties_json["atom_type"]

    return casm_lattice, frac_coords, atom_types


def casm_structure_from_properties_json(properties_json: dict) -> casm.xtal.Structure:
    """Construct a casm ``Structure`` from properties.calc.json/structure.json dictionary

    Parameters
    ----------
    properties_json : dict
        properties.calc.json dictionary with lattice, atom_coords, atom_types
        information

    Returns
    -------
    casm.xtal.Structure
        A casm ``Structure`` with lattice, atom_coords, atom_types. All the
        molecule properties (coords, types) are default values which are ``None``.

    """
    casm_lattice, frac_coords, atom_types = get_structure_info_from_properties_json(
        properties_json
    )

    return casm.xtal.Structure(casm_lattice, frac_coords, atom_types)


def casm_prim_from_properties_json(properties_json: dict) -> casm.xtal.Prim:
    """Construct a casm ``Prim`` from properties.calc.json/structure.json dictionary

    Parameters
    ----------
    properties_json : dict
        properties.calc.json dictionary with lattice, atom_coords, atom_types
        information

    Returns
    -------
    casm.xtal.Prim
        A casm ``Prim`` with lattice, atom_coords, atom_dofs. All the molecule
        properties (coords, types) are default values which are ``None``

    """
    casm_lattice, frac_coords, atom_types = get_structure_info_from_properties_json(
        properties_json
    )
    atom_dofs = [[atom_type] for atom_type in atom_types]

    return casm.xtal.Prim(casm_lattice, frac_coords, atom_dofs)


def casm_structure_from_poscar(arg1):
    """TODO: Docstring for casm_structure_from_poscar.
    Need this function is casm.xtal.Structure

    Parameters
    ----------
    arg1 : TODO

    Returns
    -------
    TODO

    """
    pass


def change_atom_types_in_casm_structure(
    casm_structure: casm.xtal.Structure, atom_types: list[str]
) -> casm.xtal.Structure:
    """For a given ``casm_structure`` change the ``atom_types`` for each site.
    Useful for masking any of the sites while mapping.

    Parameters
    ----------
    casm_structure : casm.xtal.Structure
        A casm ``Structure`` for which atom types need to be changed
    atom_types : List[str]
        List of new atom types. For example, if you have two sites in
        your crystal and want elements "A" and "B" as your atom types
        respectively at each site, then``atom_types`` will be ["A", "B"]

    Returns
    -------
    casm.xtal.Structure
        A new casm `Structure` with atom types changed

    """
    if len(atom_types) != len(casm_structure.atom_type()):
        raise RuntimeError(
            "Given number of atom types ("
            + str(len(atom_types))
            + ") is not the same as number of atoms ("
            + str(len(casm_structure.atom_type()))
            + ") in the structure"
        )
    return casm.xtal.Structure(
        casm_structure.lattice(), casm_structure.atom_coordinate_frac(), atom_types
    )


def change_dofs_in_casm_prim(
    casm_prim: casm.xtal.Prim, atom_dofs: list[list[str]]
) -> casm.xtal.Prim:
    """For a given ``casm_prim`` change the ``atom_dofs`` for each site.
    Useful for masking any of the sites while mapping.

    Parameters
    ----------
    casm_prim : casm.xtal.Prim
        A casm ``Prim`` for which occupant dofs need to be changed
    atom_dofs : List[List[str]]
        Atom dofs that need to be applied at each site. For example, if you
        have two sites in your crystal and want elements "A" and "B" at both
        sites, atom_dofs will be [["A", "B"],["A", "B"]]

    Returns
    -------
    casm.xtal.Prim
        A new casm ``Prim`` with the given ``atom_dofs``

    Raises
    ------
    RuntimeError
        If length of sites in ``casm_prim`` don't match with the length
        of provided ``atom_dofs``

    """
    if len(atom_dofs) != len(casm_prim.occ_dof()):
        raise RuntimeError(
            "Given number of atom dofs ("
            + str(len(atom_dofs))
            + ") is not the same as number of atoms ("
            + str(len(casm_prim.occ_dof()))
            + ") in the structure"
        )

    return casm.xtal.Prim(casm_prim.lattice(), casm_prim.coordinate_frac(), atom_dofs)
