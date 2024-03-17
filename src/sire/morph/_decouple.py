__all__ = ["annihilate", "decouple"]


def annihilate(mol, as_new_molecule: bool = True, map=None):
    """
    Return a merged molecule that represents the perturbation that
    completely annihilates the molecule. The returned merged molecule
    will be suitable for using in a double-annihilation free energy
    simulation, e.g. to calculate absolute binding free energies.

    Parameters
    ----------
    mol : Molecule view
        The molecule (or part of molecule) to annihilate.
        This will only annihilate the atoms in this molecule view.
        Normally, you would want to pass in the entire molecule.
    as_new_molecule : bool, optional
        Whether to return the merged molecule as a new molecule,
        or to assign a new molecule number to the result. Default is True.
    map : dict, optional
        Property map to assign properties in the returned,
        merged molecule, plus to find the properties that will be
        annihilated.

    Returns
    -------
    Molecule
        The merged molecule representing the annihilation perturbation
    """
    try:
        # make sure we have only the reference state
        mol = mol.perturbation().extract_reference(remove_ghosts=True)
    except Exception:
        pass

    from ..base import create_map
    from ..mm import LJParameter
    from ..mol import Element
    from ..units import kcal_per_mol, mod_electron, g_per_mol

    map = create_map(map)

    c = mol.cursor()
    c_mol = c.molecule()

    c["is_perturbable"] = True

    has_key = {}

    for key in [
        "charge",
        "LJ",
        "bond",
        "angle",
        "dihedral",
        "improper",
        "forcefield",
        "intrascale",
        "mass",
        "element",
        "atomtype",
        "ambertype",
        "connectivity",
    ]:
        key = map[key].source()

        if key in c:
            c_mol[f"{key}0"] = c_mol[key]
            c_mol[f"{key}1"] = c_mol[key]

            has_key[key] = True

            if key != "connectivity":
                del c_mol[key]
        else:
            has_key[key] = False

    lj_prop = map["LJ"].source()
    chg_prop = map["charge"].source()
    elem_prop = map["element"].source()
    ambtype_prop = map["ambertype"].source()
    atomtype_prop = map["atomtype"].source()
    mass_prop = map["mass"].source()

    # destroy all of the atoms
    for atom in c.atoms():
        lj = atom[f"{lj_prop}0"]

        atom[f"{lj_prop}1"] = LJParameter(lj.sigma(), 0.0 * kcal_per_mol)
        atom[f"{chg_prop}1"] = 0 * mod_electron

        if has_key[elem_prop]:
            atom[f"{elem_prop}1"] = Element(0)

        if has_key[ambtype_prop]:
            atom[f"{ambtype_prop}1"] = "Xx"

        if has_key[atomtype_prop]:
            atom[f"{atomtype_prop}1"] = "Xx"

        if has_key[mass_prop]:
            atom[f"{mass_prop}1"] = 0.0 * g_per_mol

    # now remove all of the bonds, angles, dihedrals, impropers
    for key in ["bond", "angle", "dihedral", "improper"]:
        if has_key[key]:
            p = c[f"{key}1"]
            p.clear()
            c[f"{key}1"] = p

    # we will leave the intrascale property as is, as this accounts
    # for the connectivity of this molecule, and would likely break
    # things if we scaled it with lambda (the charge and LJ are already
    # being scaled down)

    mol = c_mol.commit()

    c_mol["molecule0"] = mol.perturbation().extract_reference(remove_ghosts=True)
    c_mol["molecule1"] = mol.perturbation().extract_perturbed(remove_ghosts=True)

    if "parameters" in c_mol:
        del c_mol["parameters"]

    if "amberparams" in c_mol:
        del c_mol["amberparams"]

    if as_new_molecule:
        c_mol.renumber()

    # need to add a LambdaSchedule that could be used to decouple
    # the molecule
    from ..cas import LambdaSchedule

    # we decouple via a standard morph which does not scale the
    # intramolecular terms
    c_mol["schedule"] = LambdaSchedule.standard_annihilate(
        perturbed_is_annihilated=True
    )

    mol = c_mol.commit().perturbation().link_to_reference()

    return mol


def decouple(mol, as_new_molecule: bool = True, map=None):
    """
    Return a merged molecule that represents the perturbation that
    completely decouples the molecule. The returned merged molecule
    will be suitable for using in a double-decoupling free energy
    simulation, e.g. to calculate absolute binding free energies.

    Parameters
    ----------
    mol : Molecule view
        The molecule (or part of molecule) to decouple.
        This will only decouple the atoms in this molecule view.
        Normally, you would want to pass in the entire molecule.
    as_new_molecule : bool, optional
        Whether to return the merged molecule as a new molecule,
        or to assign a new molecule number to the result. Default is True.
    map : dict, optional
        Property map to assign properties in the returned,
        merged molecule, plus to find the properties that will be
        decoupled.

    Returns
    -------
    Molecule
        The merged molecule representing the decoupling perturbation
    """
    try:
        # make sure we have only the reference state
        mol = mol.perturbation().extract_reference(remove_ghosts=True)
    except Exception:
        pass

    from ..base import create_map
    from ..mm import LJParameter
    from ..units import kcal_per_mol, mod_electron

    map = create_map(map)

    c = mol.cursor()

    c_mol = c.molecule()
    c_mol["is_perturbable"] = True

    for key in [
        "charge",
        "LJ",
        "bond",
        "angle",
        "dihedral",
        "improper",
        "forcefield",
        "intrascale",
        "mass",
        "element",
        "atomtype",
        "ambertype",
        "connectivity",
    ]:
        key = map[key].source()

        if key in c:
            c_mol[f"{key}0"] = c_mol[key]
            c_mol[f"{key}1"] = c_mol[key]

            if key != "connectivity":
                del c_mol[key]

    lj_prop = map["LJ"].source()
    chg_prop = map["charge"].source()

    for atom in c.atoms():
        lj = atom[f"{lj_prop}0"]

        atom[f"{lj_prop}1"] = LJParameter(lj.sigma(), 0.0 * kcal_per_mol)
        atom[f"{chg_prop}1"] = 0 * mod_electron

    mol = c_mol.commit()

    c_mol["molecule0"] = mol.perturbation().extract_reference(remove_ghosts=True)
    c_mol["molecule1"] = mol.perturbation().extract_perturbed(remove_ghosts=True)

    if "parameters" in c_mol:
        del c_mol["parameters"]

    if "amberparams" in c_mol:
        del c_mol["amberparams"]

    if as_new_molecule:
        c_mol.renumber()

    # need to add a LambdaSchedule that could be used to decouple
    # the molecule
    from ..cas import LambdaSchedule

    # we decouple via a standard morph which does not scale the
    # intramolecular terms
    c_mol["schedule"] = LambdaSchedule.standard_decouple(perturbed_is_decoupled=True)

    mol = c_mol.commit().perturbation().link_to_reference()

    return mol
