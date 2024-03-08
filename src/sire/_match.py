__all__ = ["match_atoms"]


def match_atoms(
    mol0, mol1, match=None, prematch=None, match_light_atoms=False, map0=None, map1=None
):
    """
    Perform a simple match that tries to identify the mapping from
    atoms in 'mol0' to the atoms in 'mol1'. This uses the `AtomMCSMatcher`
    to match the atoms, using the passed `prematch` argument.

    However, if the `match` argument is provided, this will be used
    as the atom mapping directly (it can either be a dictionary mapping
    atom identifiers, or an `AtomMatcher` object).

    Parameters
    ----------
    mol0 : Molecule view
        The reference state molecule (or part of molecule)
    mol1 : Molecule view
        The perturbed state molecule (or part of molecule)
    match : dict, AtomMatcher, optional
        The atom matcher to use to match atoms. If this is a dictionary
        of atom identifiers, then this will be passed to a
        `AtomIDMatcher` object. If this is an `AtomMatcher` object, then
        this will be used directly.
    prematch : dict, AtomMatcher, optional
        The atom matcher to use to prematch atoms. If `match` is not
        supplied, then this will be used as the `prematch` argument
        to the `AtomMCSMatcher` used to find the maximum common subgraph
        match.
    match_light_atoms : bool, optional
        Whether to match light atoms (i.e. hydrogen atoms) if using the
        default `AtomMCSMatcher`. Default is False.
    map0 : dict, optional
        Property map to find properties in `mol0`
    map1 : dict, optional
        Property map to find properties in `mol1`

    Returns
    -------
    AtomMapping
        The atom mapping between the two molecules (or parts of molecules)
    """
    from .mol import AtomMapping
    from .legacy.Mol import AtomMCSMatcher, AtomMatcher
    from .base import create_map

    if mol0.num_molecules() > 1 or mol1.num_molecules() > 1:
        raise ValueError("You cannot match multiple molecules at once")

    from .system import System

    if System.is_system(mol0):
        mol0 = mol0[0]

    if System.is_system(mol1):
        mol1 = mol1[0]

    map0 = create_map(map0)
    map1 = create_map(map1)

    if match is not None:
        if prematch is not None:
            raise ValueError("You cannot provide both a `match` and a `prematch`")

        if not isinstance(match, AtomMatcher):
            from .legacy.Mol import AtomIDMatcher

            # create a dictionary of atom identifiers
            if isinstance(match, dict):
                from . import atomid

                matches = {}

                for atom0, atom1 in match.items():
                    if isinstance(atom0, int):
                        atom0 = atomid(idx=atom0)
                    elif isinstance(atom0, str):
                        atom0 = atomid(name=atom0)

                    if isinstance(atom1, int):
                        atom1 = atomid(idx=atom1)
                    elif isinstance(atom1, str):
                        atom1 = atomid(name=atom1)

                    matches[atom0] = atom1

                match = matches

            elif "KartografAtomMapper" in str(match.__class__):
                # use Kartograf to get the mapping - convert to RDKit then Kartograf
                from kartograf.atom_aligner import align_mol_shape
                from kartograf import KartografAtomMapper, SmallMoleculeComponent

                if not isinstance(match, KartografAtomMapper):
                    raise TypeError("match must be a KartografAtomMapper")

                from .convert import to
                from . import atomid

                rd_mol0 = to(mol0, "rdkit")
                rd_mol1 = to(mol1, "rdkit")

                k_mol0, k_mol1 = [
                    SmallMoleculeComponent.from_rdkit(m) for m in [rd_mol0, rd_mol1]
                ]

                k_0mol1 = align_mol_shape(k_mol1, ref_mol=k_mol0)

                mapping = next(match.suggest_mappings(k_mol0, k_0mol1))

                match = {}

                for k, v in mapping.componentA_to_componentB.items():
                    match[atomid(idx=k)] = atomid(idx=v)

            matcher = AtomIDMatcher(match)
        else:
            matcher = match

    elif prematch is not None:
        if not isinstance(prematch, AtomMatcher):
            from .legacy.Mol import AtomIDMatcher

            prematch = AtomIDMatcher(prematch)

        matcher = AtomMCSMatcher(
            prematcher=prematch, match_light_atoms=match_light_atoms, verbose=False
        )
    else:
        matcher = AtomMCSMatcher(match_light_atoms=match_light_atoms, verbose=False)

    m = matcher.match(mol0, map0, mol1, map1)

    atoms0 = []
    atoms1 = []

    for atom0, atom1 in m.items():
        atoms0.append(atom0.value())
        atoms1.append(atom1.value())

    return AtomMapping(
        mol0.atoms(),
        mol1.atoms(),
        mol0.molecule().atoms()[atoms0],
        mol1.molecule().atoms()[atoms1],
        map0,
        map1,
    )
