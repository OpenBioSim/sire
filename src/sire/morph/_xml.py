__all__ = ["evaluate_xml_force"]


def evaluate_xml_force(mols, xml, force):
    """
    Evaluate the custom force defined in the passed XML file.
    The passed molecules must be the ones used to create the
    OpenMM context associated with the XML file.


    Parameters
    ----------

    mols : sire.system.System, sire.mol.Molecule
        The perturbable molecular system or molecule to evaluate the force on.
        This should have already been linked to the appropriate end state.

    xml : str
        The path to the XML file containing the custom force.

    force : str
        The name of the custom force to evaluate. Options are:
        "ghost-ghost-lj", "ghost-ghost-coulomb",
        "ghost-nonghost-lj", "ghost-nonghost-coulomb", "ghost-14".

    Returns
    -------

    pairs : [(sire.mol.Atom, sire.mol.Atom)]
        The atom pairs that interacted.

    nrg_coul : [sr.units.GeneralUnit]
        The Coulomb energies for each atom pair.

    nrg_lj : [sr.units.GeneralUnit]
        The Lennard-Jones energies for each atom pair.
    """

    import xml.etree.ElementTree as ET
    import sys

    from .._measure import measure
    from ..legacy.Mol import Molecule
    from ..system import System
    from ..units import nanometer, kJ_per_mol

    # Store the name of the current module.
    module = sys.modules[__name__]

    # Validate the molecules.
    if not isinstance(mols, (System, Molecule)):
        raise TypeError(
            "'mols' must be of type 'sire.system.System' or 'sire.mol.Molecule'."
        )

    # Validate the XML file.
    if not isinstance(xml, str):
        raise TypeError("'xml' must be of type 'str'.")

    # Try to parse the XML file.
    try:
        tree = ET.parse(xml)
    except:
        raise ValueError(f"Could not parse the XML file: {xml}")

    # Validate the force type.
    if not isinstance(force, str):
        raise TypeError("'force' must be of type 'str'.")

    # Sanitize the force name.
    force = (
        force.lower()
        .replace(" ", "")
        .replace("-", "")
        .replace("_", "")
        .replace("/", "")
    )

    # Validate the force name.
    valid = [
        "ghostghostlj",
        "ghostghostcoulomb",
        "ghostnonghostlj",
        "ghostnonghostcoulomb",
        "ghost14",
    ]
    if force not in valid:
        raise ValueError(
            "'force' must be one of 'ghost-ghost-lj', 'ghost-ghost-coulomb', "
            "'ghost-nonghost-lj', 'ghost-nonghost-coulomb', or 'ghost-14'."
        )

    # Map sanitised name to the OpenMM force name in the XML.
    _force_name_map = {
        "ghostghostlj": "GhostGhostLJForce",
        "ghostghostcoulomb": "GhostGhostCoulombForce",
        "ghostnonghostlj": "GhostNonGhostLJForce",
        "ghostnonghostcoulomb": "GhostNonGhostCoulombForce",
        "ghost14": "Ghost14BondForce",
    }
    name = _force_name_map[force]

    # Loop over the forces until we find the named CustomNonbondedForce.
    is_found = False
    for force in tree.find("Forces"):
        if force.get("name") == name:
            is_found = True
            break

    # Raise an error if the force was not found.
    if not is_found:
        raise ValueError(f"Could not find the force: {name}")

    # Get the energy terms.
    terms = list(reversed(force.get("energy").split(";")[1:-1]))

    # Create a list to store the results.
    pairs = []
    nrg_coul_list = []
    nrg_lj_list = []

    # CustomNonbondedForce: ghost-ghost or ghost-nonghost.
    if name != "Ghost14BondForce":
        # Get the parameters for this force.
        parameters = [p.get("name") for p in force.find("PerParticleParameters")]

        # Get all the particle parameters.
        particles = force.find("Particles")

        # Get the two sets of particles that interact.
        set1 = [
            int(p.get("index"))
            for p in force.find("InteractionGroups")
            .find("InteractionGroup")
            .find("Set1")
        ]
        set2 = [
            int(p.get("index"))
            for p in force.find("InteractionGroups")
            .find("InteractionGroup")
            .find("Set2")
        ]

        # Get the exclusions.
        exclusions = [
            (int(e.get("p1")), int(e.get("p2")))
            for e in force.find("Exclusions").findall("Exclusion")
        ]
        for x, (i, j) in enumerate(exclusions):
            if i > j:
                exclusions[x] = (j, i)
        exclusions = set(exclusions)

        # Get the cutoff distance.
        cutoff = float(force.get("cutoff"))

        # Get the list of atoms.
        atoms = mols.atoms()

        # Loop over all particles in set1.
        for x in range(len(set1)):
            # Get the index from set1.
            i = set1[x]

            # Get the parameters for this particle.
            particle_i = particles[i]

            # Get the atom.
            atom_i = atoms[i]

            # Set the parameters for this particle.
            for k, param in enumerate(parameters):
                setattr(module, param + "1", float(particle_i.get(f"param{k + 1}")))

            # Loop over particles in set2.
            for y in range(len(set2)):
                # Get the index from set2.
                j = set2[y]

                # Ignore self-interaction.
                if i == j:
                    continue

                # Check if this pair is excluded.
                pair = (i, j) if i < j else (j, i)
                if pair in exclusions:
                    continue

                # Get the parameters for this particle.
                particle_j = particles[j]

                # Get the atom.
                atom_j = atoms[j]

                # Set the parameters for this particle.
                for k, param in enumerate(parameters):
                    setattr(module, param + "2", float(particle_j.get(f"param{k + 1}")))

                # Get the distance between the particles.
                r = measure(atom_i, atom_j).to(nanometer)

                # Atoms are within the cutoff.
                if r < cutoff:
                    # Evaluate the energy term by term.
                    for term in terms:
                        # Replace any instances of ^ with **.
                        term = term.replace("^", "**")

                        # Split the term into the result and the expression.
                        result, expression = term.split("=")

                        # Evaluate the expression.
                        setattr(module, result, eval(expression))

                    # Give energies units.
                    coul_nrg = module.coul_nrg * kJ_per_mol
                    lj_nrg = module.lj_nrg * kJ_per_mol

                    # Append the results for this pair.
                    pairs.append((atom_i, atom_j))
                    nrg_coul_list.append(coul_nrg)
                    nrg_lj_list.append(lj_nrg)

    # CustomBondForce: ghost-14.
    else:
        # Get the parameters for this force.
        parameters = [p.get("name") for p in force.find("PerBondParameters")]

        # Get all the bond parameters.
        bonds = force.find("Bonds").findall("Bond")

        # Get the list of atoms.
        atoms = mols.atoms()

        # Loop over all bonds.
        for bond in bonds:
            # Get the atoms involved in the bond.
            atom_i = atoms[int(bond.get("p1"))]
            atom_j = atoms[int(bond.get("p2"))]

            # Set the parameters for this bond.
            for k, param in enumerate(parameters):
                setattr(module, param, float(bond.get(f"param{k + 1}")))

            # Get the distance between the particles.
            r = measure(atom_i, atom_j).to(nanometer)

            # Evaluate the energy term by term.
            for term in terms:
                # Replace any instances of ^ with **.
                term = term.replace("^", "**")

                # Split the term into the result and the expression.
                result, expression = term.split("=")

                # Evaluate the expression.
                setattr(module, result, eval(expression))

            # Give energies units.
            coul_nrg = module.coul_nrg * kJ_per_mol
            lj_nrg = module.lj_nrg * kJ_per_mol

            # Append the results for this bond.
            pairs.append((atom_i, atom_j))
            nrg_coul_list.append(coul_nrg)
            nrg_lj_list.append(lj_nrg)

    # Return the results.
    return pairs, nrg_coul_list, nrg_lj_list
