import libcasm.xtal


def change_atom_types_in_casm_structure(
    casm_structure: libcasm.xtal.Structure, atom_types: list[str]
) -> libcasm.xtal.Structure:
    """For a given ``casm_structure`` change the ``atom_types`` for each site.
    Useful for masking any of the sites while mapping.

    Parameters
    ----------
    casm_structure : libcasm.xtal.Structure
        A casm ``Structure`` for which atom types need to be changed
    atom_types : List[str]
        List of new atom types. For example, if you have two sites in
        your crystal and want elements "A" and "B" as your atom types
        respectively at each site, then``atom_types`` will be ["A", "B"]

    Returns
    -------
    libcasm.xtal.Structure
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
    return libcasm.xtal.Structure(
        casm_structure.lattice(), casm_structure.atom_coordinate_frac(), atom_types
    )


def change_dofs_in_casm_prim(
    casm_prim: libcasm.xtal.Prim, atom_dofs: list[list[str]]
) -> libcasm.xtal.Prim:
    """For a given ``casm_prim`` change the ``atom_dofs`` for each site.
    Useful for masking any of the sites while mapping.

    Parameters
    ----------
    casm_prim : libcasm.xtal.Prim
        A casm ``Prim`` for which occupant dofs need to be changed
    atom_dofs : List[List[str]]
        Atom dofs that need to be applied at each site. For example, if you
        have two sites in your crystal and want elements "A" and "B" at both
        sites, atom_dofs will be [["A", "B"],["A", "B"]]

    Returns
    -------
    libcasm.xtal.Prim
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

    return libcasm.xtal.Prim(
        casm_prim.lattice(), casm_prim.coordinate_frac(), atom_dofs
    )
