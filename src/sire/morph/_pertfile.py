__all__ = ["create_from_pertfile"]


def create_from_pertfile(mol, pertfile, map=None):
    """
    Create a merged molecule from the passed molecule and pertfile.
    This will create and return a merged molecule with an initial
    and reference state that follows the instructions in the
    pertfile.

    Parameters
    ----------

    mol : sire.mol.Molecule
        The molecule to be merged

    pertfile : str
        The path to the pertfile

    Returns
    -------

    sire.mol.Molecule
        The merged molecule
    """
    from ..legacy.IO import PerturbationsLibrary

    from ..base import create_map

    map = create_map(map)

    # we operate on the whole molecule
    mol = mol.molecule()

    # find the name of the molecule in the pertfile...
    molname = None

    for line in open(pertfile).readlines():
        if line.startswith("molecule"):
            try:
                molname = line.split()[1]
                break
            except IndexError:
                pass

    if molname is None:
        raise ValueError(f"Could not find molecule name in pertfile: {pertfile}")

    perturbations = PerturbationsLibrary(pertfile)

    # get the template for the first molecule in this file
    # (there only ever seems to be a single template per file)
    template = perturbations.get_template(molname)

    # The pert file identifies using atom names - create a map to
    # identify atoms by index
    name_to_idx = {}

    for atom in mol.atoms():
        name_to_idx[atom.name().value()] = atom.index().value()

    # update all of the atoms from the data in the pertfile
    c = mol.cursor()

    chg_prop = map["charge"].source()
    lj_prop = map["LJ"].source()
    typ_prop = map["ambertype"].source()

    c["charge0"] = c[chg_prop]
    c["charge1"] = c[chg_prop]

    c["LJ0"] = c[lj_prop]
    c["LJ1"] = c[lj_prop]

    c["ambertype0"] = c[typ_prop]
    c["ambertype1"] = c[typ_prop]

    for atom in c.atoms():
        atomname = atom.name

        try:
            q0 = template.get_init_charge(atomname)
            q1 = template.get_final_charge(atomname)
            lj0 = template.get_init_lj(atomname)
            lj1 = template.get_final_lj(atomname)
            typ0 = template.get_init_type(atomname)
            typ1 = template.get_final_type(atomname)
        except Exception as e:
            print(e)
            continue

        atom["charge0"] = q0
        atom["charge1"] = q1

        atom["LJ0"] = lj0
        atom["LJ1"] = lj1

        atom["ambertype0"] = typ0
        atom["ambertype1"] = typ1

    # now update all of the internals
    bond_prop = map["bond"].source()
    ang_prop = map["angle"].source()
    dih_prop = map["dihedral"].source()
    imp_prop = map["improper"].source()

    bonds0 = mol.property(bond_prop)
    bonds1 = bonds0.clone()

    angles0 = mol.property(ang_prop)
    angles1 = angles0.clone()

    dihedrals0 = mol.property(dih_prop)
    dihedrals1 = dihedrals0.clone()

    impropers0 = mol.property(imp_prop)
    impropers1 = impropers0.clone()

    from ..legacy.MM import AmberBond, AmberAngle, AmberDihedral, AmberDihPart
    from ..cas import Symbol
    from ..mol import AtomIdx, BondID, AngleID, DihedralID, ImproperID

    r = Symbol("r")
    theta = Symbol("theta")
    phi = Symbol("phi")

    for bond in template.get_bonds():
        atom1 = bond.atom0()
        atom2 = bond.atom1()

        atom1_idx = AtomIdx(name_to_idx[atom1.value()])
        atom2_idx = AtomIdx(name_to_idx[atom2.value()])

        k0 = template.get_init_bond_k(bond)
        k1 = template.get_final_bond_k(bond)

        r0 = template.get_init_bond_r(bond)
        r1 = template.get_final_bond_r(bond)

        bonds0.set(
            BondID(atom1_idx, atom2_idx),
            AmberBond(k0, r0).to_expression(r),
        )

        bonds1.set(
            BondID(atom1_idx, atom2_idx),
            AmberBond(k1, r1).to_expression(r),
        )

    c["bond0"] = bonds0
    c["bond1"] = bonds1

    for angle in template.get_angles():
        atom1 = angle.atom0()
        atom2 = angle.atom1()
        atom3 = angle.atom2()

        atom1_idx = AtomIdx(name_to_idx[atom1.value()])
        atom2_idx = AtomIdx(name_to_idx[atom2.value()])
        atom3_idx = AtomIdx(name_to_idx[atom3.value()])

        k0 = template.get_init_angle_k(angle)
        k1 = template.get_final_angle_k(angle)

        theta0 = template.get_init_angle_t(angle)
        theta1 = template.get_final_angle_t(angle)

        angles0.set(
            AngleID(atom1_idx, atom2_idx, atom3_idx),
            AmberAngle(k0, theta0).to_expression(theta),
        )

        angles1.set(
            AngleID(atom1_idx, atom2_idx, atom3_idx),
            AmberAngle(k1, theta1).to_expression(theta),
        )

    c["angle0"] = angles0
    c["angle1"] = angles1

    for dihedral in template.get_dihedrals():
        atom1 = dihedral.atom0()
        atom2 = dihedral.atom1()
        atom3 = dihedral.atom2()
        atom4 = dihedral.atom3()

        atom1_idx = AtomIdx(name_to_idx[atom1.value()])
        atom2_idx = AtomIdx(name_to_idx[atom2.value()])
        atom3_idx = AtomIdx(name_to_idx[atom3.value()])
        atom4_idx = AtomIdx(name_to_idx[atom4.value()])

        params0 = template.get_init_dih_params(dihedral)
        params1 = template.get_final_dih_params(dihedral)

        func0 = 0
        func1 = 0

        for i in range(0, len(params0), 3):
            k0 = params0[i]
            n0 = params0[i + 1]
            phi0 = params0[i + 2]

            func0 += AmberDihedral(AmberDihPart(k0, n0, phi0)).to_expression(phi)

        for i in range(0, len(params1), 3):
            k1 = params1[i]
            n1 = params1[i + 1]
            phi1 = params1[i + 2]

            func1 += AmberDihedral(AmberDihPart(k1, n1, phi1)).to_expression(phi)

        dihedrals0.set(DihedralID(atom1_idx, atom2_idx, atom3_idx, atom4_idx), func0)
        dihedrals1.set(DihedralID(atom1_idx, atom2_idx, atom3_idx, atom4_idx), func1)

    c["dihedral0"] = dihedrals0
    c["dihedral1"] = dihedrals1

    for improper in template.get_impropers():
        atom1 = improper.atom0()
        atom2 = improper.atom1()
        atom3 = improper.atom2()
        atom4 = improper.atom3()

        atom1_idx = AtomIdx(name_to_idx[atom1.value()])
        atom2_idx = AtomIdx(name_to_idx[atom2.value()])
        atom3_idx = AtomIdx(name_to_idx[atom3.value()])
        atom4_idx = AtomIdx(name_to_idx[atom4.value()])

        params0 = template.get_init_imp_params(improper)
        params1 = template.get_final_imp_params(improper)

        func0 = 0
        func1 = 0

        for i in range(0, len(params0), 3):
            k0 = params0[i]
            n0 = params0[i + 1]
            phi0 = params0[i + 2]

            func0 += AmberDihedral(AmberDihPart(k0, n0, phi0)).to_expression(phi)

        for i in range(0, len(params1), 3):
            k1 = params1[i]
            n1 = params1[i + 1]
            phi1 = params1[i + 2]

            func1 += AmberDihedral(AmberDihPart(k1, n1, phi1)).to_expression(phi)

        impropers0.set(ImproperID(atom1_idx, atom2_idx, atom3_idx, atom4_idx), func0)
        impropers1.set(ImproperID(atom1_idx, atom2_idx, atom3_idx, atom4_idx), func1)

    c["improper0"] = impropers0
    c["improper1"] = impropers1

    # duplicate the coordinates, mass, and element properties
    for prop in ["coordinates", "mass", "element", "forcefield", "intrascale"]:
        orig_prop = map[prop].source()
        c[prop + "0"] = c[orig_prop]
        c[prop + "1"] = c[orig_prop]
        del c[orig_prop]

    # now remove all of the "default" properties
    del c[chg_prop]
    del c[lj_prop]
    del c[typ_prop]
    del c[bond_prop]
    del c[ang_prop]
    del c[dih_prop]
    del c[imp_prop]

    mol = c.commit()

    # we now need to generat the AmberParameters for the two end states
    from ..legacy.MM import AmberParams

    map0 = {
        "charge": "charge0",
        "LJ": "LJ0",
        "ambertype": "ambertype0",
        "bond": "bond0",
        "angle": "angle0",
        "dihedral": "dihedral0",
        "improper": "improper0",
    }

    params0 = AmberParams(mol, map0)

    map1 = {
        "charge": "charge1",
        "LJ": "LJ1",
        "ambertype": "ambertype1",
        "bond": "bond1",
        "angle": "angle1",
        "dihedral": "dihedral1",
        "improper": "improper1",
    }

    params1 = AmberParams(mol, map1)

    c["parameters0"] = params0
    c["parameters1"] = params1

    c["is_perturbable"] = True

    # make sure that we link to the reference state
    return c.commit().perturbation().link_to_reference(auto_commit=True)
