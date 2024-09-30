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
        "ghost-ghost", "ghost-nonghost", "ghost-14".

    Returns
    -------

    [((sire.mol.Atom, sr.mol.Atom), (sr.units.GeneralUnit, sr.units.GeneralUnit))]
        The atom pairs and pairwise Coulomb and Lennard-Jones energies.
    """

    from math import sqrt

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
    if not force in ["ghostghost", "ghostnonghost", "ghost14"]:
        raise ValueError(
            "'force' must be one of 'ghost-ghost', 'ghost-nonghost', or 'ghost-14'."
        )

    # Create the name and index based on the force type.
    if force == "ghostghost":
        name = "CustomNonbondedForce"
        index = 0
    elif force == "ghostnonghost":
        name = "CustomNonbondedForce"
        index = 1
    elif force == "ghost14":
        name = "CustomBondForce"
        index = 0

    # Get the root of the XML tree.
    root = tree.getroot()

    # Loop over the forces until we find the first CustomNonbondedForce.
    force_index = 0
    for force in tree.find("Forces"):
        if force.get("name") == name:
            if force_index == index:
                break
            force_index += 1

    # Get the energy terms.
    terms = list(reversed(force.get("energy").split(";")[1:-1]))

    # Create a list to store the results.
    results = []

    # CustomNonbondedForce: ghost-ghost or ghost-nonghost.
    if name == "CustomNonbondedForce":
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
            setattr(module, parameters[0] + "1", float(particle_i.get("param1")))
            setattr(module, parameters[1] + "1", float(particle_i.get("param2")))
            setattr(module, parameters[2] + "1", float(particle_i.get("param3")))
            setattr(module, parameters[3] + "1", float(particle_i.get("param4")))
            setattr(module, parameters[4] + "1", float(particle_i.get("param5")))

            # Loop over all other particles in set1.
            for y in range(x + 1, len(set1)):
                # Get the index from set2.
                j = set1[y]

                # Check if this pair is excluded.
                pair = (i, j) if i < j else (j, i)
                if pair in exclusions:
                    continue

                # Get the parameters for this particle.
                particle_j = particles[j]

                # Get the atom.
                atom_j = atoms[j]

                # Set the parameters for this particle.
                setattr(module, parameters[0] + "2", float(particle_j.get("param1")))
                setattr(module, parameters[1] + "2", float(particle_j.get("param2")))
                setattr(module, parameters[2] + "2", float(particle_j.get("param3")))
                setattr(module, parameters[3] + "2", float(particle_j.get("param4")))
                setattr(module, parameters[4] + "2", float(particle_j.get("param5")))

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
                    results.append(((atom_i, atom_j), (coul_nrg, lj_nrg)))

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
            setattr(module, parameters[0], float(bond.get("param1")))
            setattr(module, parameters[1], float(bond.get("param2")))
            setattr(module, parameters[2], float(bond.get("param3")))
            setattr(module, parameters[3], float(bond.get("param4")))
            setattr(module, parameters[4], float(bond.get("param5")))

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
            results.append(((atom_i, atom_j), (coul_nrg, lj_nrg)))

    # Return the results.
    return results
