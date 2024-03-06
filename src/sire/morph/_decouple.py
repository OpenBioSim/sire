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
    pass


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

    if as_new_molecule:
        c_mol.renumber()

    return c.commit()
