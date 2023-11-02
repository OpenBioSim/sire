
#include "openmmmolecule.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomproperty.hpp"
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"
#include "SireMol/bondorder.h"
#include "SireMol/atomvelocities.h"

#include "SireMM/atomljs.h"
#include "SireMM/selectorbond.h"
#include "SireMM/amberparams.h"
#include "SireMM/twoatomfunctions.h"

#include "SireMaths/vector.h"

#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "tostring.h"

#include <QSet>

#include <QReadWriteLock>

#include <QDebug>

using namespace SireOpenMM;
using namespace SireMM;
using namespace SireMol;
using namespace SireBase;

////////
//////// Implementation of OpenMMMolecule
////////

OpenMMMolecule::OpenMMMolecule()
{
}

OpenMMMolecule::OpenMMMolecule(const Molecule &mol,
                               const PropertyMap &map)
{
    molinfo = mol.info();
    number = mol.number();

    if (molinfo.nAtoms() == 0)
    {
        // nothing to extract
        molinfo = MoleculeInfo();
        return;
    }

    bool is_perturbable = false;
    bool swap_end_states = false;

    if (mol.hasProperty(map["is_perturbable"]))
    {
        is_perturbable = mol.property(map["is_perturbable"]).asABoolean();

        if (map.specified("swap_end_states"))
        {
            swap_end_states = map["swap_end_states"].value().asABoolean();
        }
    }

    if (is_perturbable)
    {
        ffinfo = mol.property(map["forcefield0"]).asA<MMDetail>();
    }
    else
    {
        ffinfo = mol.property(map["forcefield"]).asA<MMDetail>();
    }

    if (map.specified("constraint"))
    {
        const auto c = map["constraint"].source().toLower().simplified();

        if (c == "none")
        {
            constraint_type = CONSTRAIN_NONE;
        }
        else if (c == "h-bonds" or c == "h_bonds")
        {
            constraint_type = CONSTRAIN_HBONDS;
        }
        else if (c == "bonds")
        {
            constraint_type = CONSTRAIN_BONDS;
        }
        else if (c == "h-bonds-h-angles" or c == "h_bonds_h_angles")
        {
            constraint_type = CONSTRAIN_HBONDS | CONSTRAIN_HANGLES;
        }
        else if (c == "bonds-h-angles" or c == "bonds_h_angles")
        {
            constraint_type = CONSTRAIN_BONDS | CONSTRAIN_HANGLES;
        }
        else
        {
            throw SireError::invalid_key(QObject::tr(
                                             "Unrecognised constraint type '%1'. Valid values are "
                                             "'none', 'h-bonds', 'bonds', 'h-bonds-h-angles' or "
                                             "'bonds-h-angles',")
                                             .arg(c),
                                         CODELOC);
        }
    }
    else
    {
        constraint_type = CONSTRAIN_NONE;
    }

    if (map.specified("perturbable_constraint"))
    {
        const auto c = map["perturbable_constraint"].source().toLower().simplified();

        if (c == "none")
        {
            perturbable_constraint_type = CONSTRAIN_NONE;
        }
        else if (c == "h-bonds")
        {
            perturbable_constraint_type = CONSTRAIN_HBONDS;
        }
        else if (c == "bonds")
        {
            perturbable_constraint_type = CONSTRAIN_BONDS;
        }
        else if (c == "h-bonds-h-angles")
        {
            perturbable_constraint_type = CONSTRAIN_HBONDS | CONSTRAIN_HANGLES;
        }
        else if (c == "bonds-h-angles")
        {
            perturbable_constraint_type = CONSTRAIN_BONDS | CONSTRAIN_HANGLES;
        }
        else
        {
            throw SireError::invalid_key(QObject::tr(
                                             "Unrecognised perturbable constraint type '%1'. Valid values are "
                                             "'none', 'h-bonds', 'bonds', 'h-bonds-h-angles' or "
                                             "'bonds-h-angles',")
                                             .arg(c),
                                         CODELOC);
        }
    }
    else
    {
        perturbable_constraint_type = constraint_type;
    }

    // extract out the indicies of any field atoms, if there are any.
    // Field atoms are those that only exist via the interaction field
    // they produce (i.e. their potentials are combined into a grid,
    // which affects the energy/forces of other particles that are
    // on that grid)
    const auto field_atom_prop = map["is_field_atom"];

    if (mol.hasProperty(field_atom_prop))
    {
        const auto atoms = mol.atoms("atom property is_field_atom == True");

        for (int i = 0; i < atoms.count(); ++i)
        {
            field_atoms.append(atoms(i).index().value());
        }

        // should not need to sort it, but just in case
        qSort(field_atoms);
    }

    if (ffinfo.isAmberStyle())
    {
        if (is_perturbable)
        {
            // update the map to find the lambda=0 properties
            // Note that we don't use the coordinates0 or coordinates1
            // properties, because we need to build the molecule from
            // its current coordinates (which should represent the
            // current lambda state)
            QStringList props = {"LJ", "ambertype", "angle", "atomtype",
                                 "bond", "charge",
                                 "dihedral", "element", "forcefield",
                                 "gb_radii", "gb_screening", "improper",
                                 "intrascale", "mass", "name",
                                 "parameters", "treechain"};

            // we can't specialise these globally in case other molecules
            // are not of amber type
            auto map0 = map.addSuffix("0", props);
            auto map1 = map.addSuffix("1", props);

            if (swap_end_states)
            {
                std::swap(map0, map1);
            }

            // save this perturbable map - this will help us set
            // new properties from the results of dynamics, e.g.
            // updating coordinates after minimisation
            perturtable_map = map0;

            // extract the parameters in amber format - this should work,
            // as the 'forcefield' property has promised that this is
            // an amber-style molecule
            const auto params = SireMM::AmberParams(mol, map0);
            const auto params1 = SireMM::AmberParams(mol, map1);

            perturbed.reset(new OpenMMMolecule(*this));
            perturbed->constructFromAmber(mol, params1, params, map1, true);

            this->constructFromAmber(mol, params, params1, map0, true);

            this->alignInternals(map);
        }
        else
        {
            const auto params = SireMM::AmberParams(mol, map);
            this->constructFromAmber(mol, params, params, map, false);
        }
    }
    else
    {
        throw SireError::unsupported(QObject::tr(
                                         "We currently only support creating OpenMM molecules from "
                                         "molecules that have been parameterised using Amber-style "
                                         "forcefields. The molecule %1 has been parameterised using "
                                         "the forcefield %2.")
                                         .arg(mol.toString())
                                         .arg(ffinfo.toString()),
                                     CODELOC);
    }

    if (this->hasFieldAtoms())
    {
        this->processFieldAtoms();
    }
}

OpenMMMolecule::~OpenMMMolecule()
{
}

bool OpenMMMolecule::operator==(const OpenMMMolecule &other) const
{
    return this == &other;
}

bool OpenMMMolecule::operator!=(const OpenMMMolecule &other) const
{
    return not this->operator==(other);
}

bool OpenMMMolecule::isPerturbable() const
{
    return perturbed.get() != 0;
}

bool OpenMMMolecule::isGhostAtom(int atom) const
{
    return from_ghost_idxs.contains(atom) or to_ghost_idxs.contains(atom);
}

bool OpenMMMolecule::isFieldAtom(int atom) const
{
    return field_atoms.contains(atom);
}

bool OpenMMMolecule::isFieldMolecule() const
{
    return field_atoms.count() == atoms.count();
}

bool OpenMMMolecule::hasFieldAtoms() const
{
    return not field_atoms.isEmpty();
}

OpenMM::Vec3 to_vec3(const SireMaths::Vector &coords)
{
    const double internal_to_nm = (1 * SireUnits::angstrom).to(SireUnits::nanometer);

    return OpenMM::Vec3(internal_to_nm * coords.x(),
                        internal_to_nm * coords.y(),
                        internal_to_nm * coords.z());
}

OpenMM::Vec3 to_vec3(const SireMol::Velocity3D &vel)
{
    return OpenMM::Vec3(vel.x().to(SireUnits::nanometers_per_ps),
                        vel.y().to(SireUnits::nanometers_per_ps),
                        vel.z().to(SireUnits::nanometers_per_ps));
}

inline qint64 to_pair(qint64 x, qint64 y)
{
    if (y < x)
        return to_pair(y, x);
    else
        return x << 32 | y & 0x00000000FFFFFFFF;
}

std::tuple<int, int, double, double, double> OpenMMMolecule::getException(
    int atom0, int atom1, int start_index, double coul_14_scl, double lj_14_scl) const
{
    double charge = 0.0;
    double sigma = 0.0;
    double epsilon = 0.0;

    if (coul_14_scl != 0 or lj_14_scl != 0)
    {
        if (atom0 < 0 or atom0 >= cljs.count() or atom1 < 0 or atom1 >= cljs.count())
            throw SireError::invalid_index(QObject::tr(
                                               "Cannot get CLJ parameters for atom %1 or atom %2.")
                                               .arg(atom0)
                                               .arg(atom1),
                                           CODELOC);

        const auto &clj0 = cljs.constData()[atom0];
        const auto &clj1 = cljs.constData()[atom1];

        charge = coul_14_scl * std::get<0>(clj0) * std::get<0>(clj1);
        sigma = 0.5 * (std::get<1>(clj0) + std::get<1>(clj1));
        epsilon = lj_14_scl * std::sqrt(std::get<2>(clj0) * std::get<2>(clj1));
    }

    if (this->isPerturbable() and charge == 0 and epsilon == 0)
    {
        // openmm tries to optimise away zero parameters - this is an issue
        // as perturbation requires that we don't remove them!
        // If we don't do this, then we get a
        // "updateParametersInContext: The set of non-excluded exceptions has changed"
        /// exception when we update parameters in context
        sigma = 1e-9;
        epsilon = 1e-9;
    }

    return std::make_tuple(atom0 + start_index,
                           atom1 + start_index,
                           charge, sigma, epsilon);
}

/** Return closest constraint length to 'length' based on what
 *  we have seen before and the constraint_length_tolerance
 */
double getSharedConstraintLength(double length)
{
    // distance below which constraints are considered to be equal (e.g. to
    // r0 or to another angle - units in nanometers)
    const double constraint_length_tolerance = 0.005;

    // is this close to a O-H bond length of water?
    if (std::abs(length - 0.09572) < constraint_length_tolerance)
    {
        return 0.09572;
    }

    static QVector<double> angle_constraint_lengths;
    static QReadWriteLock l;

    // most of the time we expect a hit, so a read lock is ok
    QReadLocker locker(&l);

    // is this close to any of the existing angle constraints?
    for (const auto &cl : angle_constraint_lengths)
    {
        if (std::abs(cl - length) < constraint_length_tolerance)
        {
            // this is close enough to an existing constraint
            // so we will use that instead
            return cl;
        }
    }

    locker.unlock();

    // ok, we didn't find it - we should append it to the list
    // unless another thread has got here first
    QWriteLocker locker2(&l);

    for (const auto &cl : angle_constraint_lengths)
    {
        if (std::abs(cl - length) < constraint_length_tolerance)
        {
            // this is close enough to an existing constraint
            // so we will use that instead
            return cl;
        }
    }

    // this is a new constraint, so add it to the list
    angle_constraint_lengths.append(length);

    return length;
}

void OpenMMMolecule::constructFromAmber(const Molecule &mol,
                                        const AmberParams &params,
                                        const AmberParams &params1,
                                        const PropertyMap &map,
                                        bool is_perturbable)
{
    const auto &moldata = mol.data();
    atoms = mol.atoms();
    const int nats = atoms.count();

    if (nats <= 0)
    {
        return;
    }

    // look up the CGAtomIdx of each atom - this is because we
    // will use AtomIdx for the ordering and atom identifiers
    auto idx_to_cgatomidx = QVector<SireMol::CGAtomIdx>(nats);
    auto idx_to_cgatomidx_data = idx_to_cgatomidx.data();

    for (int i = 0; i < nats; ++i)
    {
        idx_to_cgatomidx_data[i] = molinfo.cgAtomIdx(SireMol::AtomIdx(i));
    }

    // extract the coordinates and convert to OpenMM units
    const auto &c = moldata.property(map["coordinates"]).asA<SireMol::AtomCoords>();

    this->coords = QVector<OpenMM::Vec3>(nats, OpenMM::Vec3(0, 0, 0));
    auto coords_data = coords.data();

    for (int i = 0; i < nats; ++i)
    {
        coords_data[i] = to_vec3(c.at(idx_to_cgatomidx_data[i]));
    }

    // extract the velocities and convert to OpenMM units
    const auto vels_prop = map["velocities"];
    vels = QVector<OpenMM::Vec3>(nats, OpenMM::Vec3(0, 0, 0));

    if (moldata.hasProperty(vels_prop))
    {
        const auto &v = moldata.property(vels_prop).asA<SireMol::AtomVelocities>();
        auto vels_data = vels.data();

        for (int i = 0; i < nats; ++i)
        {
            vels_data[i] = to_vec3(v.at(idx_to_cgatomidx_data[i]));
        }
    }

    // extract the masses and convert to OpenMM units
    const auto params_masses = params.masses();

    this->masses = QVector<double>(nats, 0.0);

    auto masses_data = masses.data();

    if (is_perturbable)
    {
        const auto params1_masses = params1.masses();

        for (int i = 0; i < nats; ++i)
        {
            const auto cgatomidx = idx_to_cgatomidx_data[i];

            // use the largest mass of the perturbing atoms
            double mass = std::max(params_masses.at(cgatomidx).to(SireUnits::g_per_mol),
                                   params1_masses.at(cgatomidx).to(SireUnits::g_per_mol));

            if (mass < 0.05)
            {
                mass = 0.0;
            }

            if (mass < 0.5)
            {
                virtual_sites.append(i);
            }
            else if (mass < 2.5)
            {
                light_atoms.append(i);
            }

            masses_data[i] = mass;
        }
    }
    else
    {
        for (int i = 0; i < nats; ++i)
        {
            double mass = params_masses.at(idx_to_cgatomidx_data[i]).to(SireUnits::g_per_mol);

            if (mass < 0.05)
            {
                mass = 0.0;
            }

            if (mass < 0.5)
            {
                virtual_sites.append(i);
            }
            else if (mass < 2.5)
            {
                light_atoms.append(i);
            }

            masses_data[i] = mass;
        }
    }

    // extract the charges and LJ parameters and convert to OpenMM units
    const auto params_charges = params.charges();
    const auto params_ljs = params.ljs();

    this->cljs = QVector<std::tuple<double, double, double>>(nats, std::make_tuple(0.0, 0.0, 0.0));
    auto cljs_data = cljs.data();

    for (int i = 0; i < nats; ++i)
    {
        const auto &cgatomidx = idx_to_cgatomidx_data[i];

        const double chg = params_charges.at(idx_to_cgatomidx_data[i]).to(SireUnits::mod_electron);

        const auto &lj = params_ljs.at(idx_to_cgatomidx_data[i]);
        const double sig = lj.sigma().to(SireUnits::nanometer);
        const double eps = lj.epsilon().to(SireUnits::kJ_per_mol);

        cljs_data[i] = std::make_tuple(chg, sig, eps);
    }

    this->bond_params.clear();
    this->constraints.clear();

    // initialise all atoms as being unbonded
    this->unbonded_atoms.reserve(nats);

    for (int i = 0; i < nats; ++i)
    {
        this->unbonded_atoms.insert(i);
    }

    // now the bonds
    const double bond_k_to_openmm = 2.0 * (SireUnits::kcal_per_mol / (SireUnits::angstrom * SireUnits::angstrom)).to(SireUnits::kJ_per_mol / (SireUnits::nanometer * SireUnits::nanometer));
    const double bond_r0_to_openmm = SireUnits::angstrom.to(SireUnits::nanometer);

    QSet<qint64> constrained_pairs;

    bool include_constrained_energies = true;

    if (map.specified("include_constrained_energies"))
    {
        include_constrained_energies = map["include_constrained_energies"].value().asABoolean();
    }

    for (auto it = params.bonds().constBegin();
         it != params.bonds().constEnd();
         ++it)
    {
        const auto bondid = it.key().map(molinfo);
        const auto &bondparam = it.value().first;

        int atom0 = bondid.get<0>().value();
        int atom1 = bondid.get<1>().value();

        if (atom0 > atom1)
            std::swap(atom0, atom1);

        const double k = bondparam.k() * bond_k_to_openmm;
        const double r0 = bondparam.r0() * bond_r0_to_openmm;

        if (k != 0)
        {
            this->unbonded_atoms.remove(atom0);
            this->unbonded_atoms.remove(atom1);
        }

        const bool has_light_atom = (masses_data[atom0] < 2.5 or masses_data[atom1] < 2.5);
        const bool has_massless_atom = masses_data[atom0] < 0.5 or masses_data[atom1] < 0.5;

        auto this_constraint_type = constraint_type;

        if (is_perturbable)
        {
            this_constraint_type = perturbable_constraint_type;
        }

        bool bond_is_not_constrained = true;

        if ((not has_massless_atom) and ((this_constraint_type & CONSTRAIN_BONDS) or (has_light_atom and (this_constraint_type & CONSTRAIN_HBONDS))))
        {
            // add the constraint - this constrains the bond to whatever length it has now
            const auto delta = coords[atom1] - coords[atom0];
            auto constraint_length = std::sqrt((delta[0] * delta[0]) +
                                               (delta[1] * delta[1]) +
                                               (delta[2] * delta[2]));

            // use the r0 for the bond if this is close to the measured length and this
            // is not a perturbable molecule
            if (not is_perturbable and std::abs(constraint_length - r0) < 0.01)
            {
                constraint_length = r0;
            }

            this->constraints.append(std::make_tuple(atom0, atom1, constraint_length));
            constrained_pairs.insert(to_pair(atom0, atom1));
            bond_is_not_constrained = false;
        }

        if (include_constrained_energies or bond_is_not_constrained)
            this->bond_params.append(std::make_tuple(atom0, atom1, r0, k));
    }

    // now the angles
    const double angle_k_to_openmm = 2.0 * (SireUnits::kcal_per_mol).to(SireUnits::kJ_per_mol);

    ang_params.clear();

    for (auto it = params.angles().constBegin();
         it != params.angles().constEnd();
         ++it)
    {
        const auto angid = it.key().map(molinfo);
        const auto &angparam = it.value().first;

        int atom0 = angid.get<0>().value();
        int atom1 = angid.get<1>().value();
        int atom2 = angid.get<2>().value();

        if (atom0 > atom2)
            std::swap(atom0, atom2);

        const double k = angparam.k() * angle_k_to_openmm;
        const double theta0 = angparam.theta0(); // already in radians

        const bool is_h_x_h = masses_data[atom0] < 2.5 and masses_data[atom2] < 2.5;

        const auto key = to_pair(atom0, atom2);

        bool angle_is_not_constrained = true;

        if (not constrained_pairs.contains(key))
        {
            auto this_constraint_type = constraint_type;

            if (is_perturbable)
            {
                this_constraint_type = perturbable_constraint_type;
            }

            // only include the angle X-y-Z if X-Z are not already constrained
            if ((this_constraint_type & CONSTRAIN_HANGLES) and is_h_x_h)
            {
                const auto delta = coords[atom2] - coords[atom0];
                auto constraint_length = std::sqrt((delta[0] * delta[0]) +
                                                   (delta[1] * delta[1]) +
                                                   (delta[2] * delta[2]));

                // we can speed up OpenMM by making sure that constraints are
                // equal if they operate on similar molecules (e.g. all water
                // constraints are the same. We will check for this if this is
                // a non-perturbable small molecule
                if (mol.nAtoms() < 10 and not is_perturbable)
                {
                    constraint_length = getSharedConstraintLength(constraint_length);
                }

                constraints.append(std::make_tuple(atom0, atom2,
                                                   constraint_length));
                constrained_pairs.insert(key);
                angle_is_not_constrained = false;
            }
        }
        else
            angle_is_not_constrained = false;

        if (include_constrained_energies or angle_is_not_constrained)
            ang_params.append(std::make_tuple(atom0, atom1, atom2,
                                              theta0, k));
    }

    // now the dihedrals
    const double dihedral_k_to_openmm = (SireUnits::kcal_per_mol).to(SireUnits::kJ_per_mol);

    dih_params.clear();

    for (auto it = params.dihedrals().constBegin();
         it != params.dihedrals().constEnd();
         ++it)
    {
        const auto dihid = it.key().map(molinfo);
        const auto &dihparam = it.value().first;

        int atom0 = dihid.get<0>().value();
        int atom1 = dihid.get<1>().value();
        int atom2 = dihid.get<2>().value();
        int atom3 = dihid.get<3>().value();

        if (atom0 > atom3)
        {
            std::swap(atom0, atom3);
            std::swap(atom1, atom2);
        }

        for (const auto &term : dihparam.terms())
        {
            const double v = term.k() * dihedral_k_to_openmm;
            const double phase = term.phase();                      // already in radians
            const int periodicity = std::round(term.periodicity()); // this is just an integer

            if (periodicity > 0)
            {
                dih_params.append(std::make_tuple(atom0, atom1,
                                                  atom2, atom3,
                                                  periodicity, phase, v));
            }
            else if (periodicity == 0 and v == 0)
            {
                // this is a null dihedral, e.g. for perturbation. Make sure
                // that the periodicity is 1, else otherwise openmm will complain
                dih_params.append(std::make_tuple(atom0, atom1,
                                                  atom2, atom3,
                                                  1, phase, v));
            }
            else
            {
                qWarning() << "IGNORING DIHEDRAL WITH WEIRD PERIODICITY" << dihparam.toString();
            }
        }
    }

    // now the impropers
    for (auto it = params.impropers().constBegin();
         it != params.impropers().constEnd();
         ++it)
    {
        const auto impid = it.key().map(molinfo);
        const auto &impparam = it.value().first;

        const int atom0 = impid.get<0>().value();
        const int atom1 = impid.get<1>().value();
        const int atom2 = impid.get<2>().value();
        const int atom3 = impid.get<3>().value();

        for (const auto &term : impparam.terms())
        {
            const double v = term.k() * dihedral_k_to_openmm;
            const double phase = term.phase();                      // already in radians
            const int periodicity = std::round(term.periodicity()); // this is just an integer

            if (periodicity > 0)
            {
                dih_params.append(std::make_tuple(atom0, atom1,
                                                  atom2, atom3,
                                                  periodicity, phase, v));
            }
            else if (periodicity == 0 and v == 0)
            {
                // this is a null dihedral, e.g. for perturbation. Make sure
                // that the periodicity is 1, else otherwise openmm will complain
                dih_params.append(std::make_tuple(atom0, atom1,
                                                  atom2, atom3,
                                                  1, phase, v));
            }
            else
            {
                qWarning() << "IGNORING IMPROPER WITH WEIRD PERIODICITY" << impparam.toString();
            }
        }
    }

    this->buildExceptions(mol, constrained_pairs, map);
}

bool is_ghost(const std::tuple<double, double, double> &clj)
{
    return std::get<0>(clj) == 0 and (std::get<1>(clj) == 0 or std::get<2>(clj) == 0);
}

bool is_ghost_charge(const std::tuple<double, double, double> &clj)
{
    return std::get<0>(clj) == 0;
}

bool is_ghost_lj(const std::tuple<double, double, double> &clj)
{
    return std::get<1>(clj) == 0 or std::get<2>(clj) == 0;
}

/** Go through all of the internals and compare them to the perturbed
 *  state. Ensure that there is a one-to-one mapping, with them all
 *  in the same order. Any that are missing are added as nulls in
 *  the correct end state.
 */
void OpenMMMolecule::alignInternals(const PropertyMap &map)
{
    // first go through an see which atoms are ghosts
    // While we do this, set the alpha values for the
    // reference and perturbed molecules
    if (cljs.count() != perturbed->cljs.count())
        throw SireError::incompatible_error(QObject::tr(
                                                "Different number of CLJ parameters between the reference "
                                                "(%1) and perturbed (%2) states.")
                                                .arg(cljs.count())
                                                .arg(perturbed->cljs.count()),
                                            CODELOC);

    this->alphas = QVector<double>(cljs.count(), 0.0);
    this->perturbed->alphas = this->alphas;

    for (int i = 0; i < cljs.count(); ++i)
    {
        const auto &clj0 = cljs.at(i);
        const auto &clj1 = perturbed->cljs.at(i);

        if (clj0 != clj1)
        {
            if (is_ghost(clj0))
            {
                from_ghost_idxs.insert(i);
                this->alphas[i] = 1.0;
            }
            else if (is_ghost(clj1))
            {
                to_ghost_idxs.insert(i);
                this->perturbed->alphas[i] = 1.0;
            }
        }
    }

    QVector<std::tuple<int, int, double, double>> bond_params_1;
    bond_params_1.reserve(bond_params.count());

    QVector<bool> found_index_0(bond_params.count(), false);
    QVector<bool> found_index_1(perturbed->bond_params.count(), false);

    for (int i = 0; i < bond_params.count(); ++i)
    {
        const auto &bond0 = bond_params.at(i);

        int atom0 = std::get<0>(bond0);
        int atom1 = std::get<1>(bond0);

        bool found = false;

        for (int j = 0; j < perturbed->bond_params.count(); ++j)
        {
            if (not found_index_1[j])
            {
                const auto &bond1 = perturbed->bond_params.at(j);

                if (std::get<0>(bond1) == atom0 and std::get<1>(bond1) == atom1)
                {
                    // we have found the matching bond!
                    bond_params_1.append(bond1);
                    found_index_0[i] = true;
                    found_index_1[j] = true;
                    found = true;
                    break;
                }
            }
        }

        if (not found)
        {
            // add a null bond with the same r0, but null k
            found_index_0[i] = true;
            bond_params_1.append(std::tuple<int, int, double, double>(atom0, atom1, std::get<2>(bond0), 0.0));
        }
    }

    for (int j = 0; j < perturbed->bond_params.count(); ++j)
    {
        if (not found_index_1[j])
        {
            // need to add a bond missing in the reference state
            const auto &bond1 = perturbed->bond_params.at(j);

            int atom0 = std::get<0>(bond1);
            int atom1 = std::get<1>(bond1);

            // add a null bond with the same r0, but null k
            bond_params.append(std::tuple<int, int, double, double>(atom0, atom1, std::get<2>(bond1), 0.0));
            bond_params_1.append(bond1);

            found_index_1[j] = true;

            unbonded_atoms.remove(atom0);
            unbonded_atoms.remove(atom1);
        }
    }

    // all of found_index_0 and found_index_1 should be true...
    if (found_index_0.indexOf(false) != -1 or found_index_1.indexOf(false) != -1)
    {
        throw SireError::program_bug(QObject::tr(
                                         "Failed to align the bonds!"),
                                     CODELOC);
    }

    perturbed->bond_params = bond_params_1;

    QVector<std::tuple<int, int, int, double, double>> ang_params_1;
    ang_params_1.reserve(ang_params.count());

    found_index_0 = QVector<bool>(ang_params.count(), false);
    found_index_1 = QVector<bool>(perturbed->ang_params.count(), false);

    for (int i = 0; i < ang_params.count(); ++i)
    {
        const auto &ang0 = ang_params.at(i);

        int atom0 = std::get<0>(ang0);
        int atom1 = std::get<1>(ang0);
        int atom2 = std::get<2>(ang0);

        bool found = false;

        for (int j = 0; j < perturbed->ang_params.count(); ++j)
        {
            if (not found_index_1[j])
            {
                const auto &ang1 = perturbed->ang_params.at(j);

                if (std::get<0>(ang1) == atom0 and std::get<1>(ang1) == atom1 and std::get<2>(ang1) == atom2)
                {
                    // we have found the matching angle!
                    ang_params_1.append(ang1);
                    found_index_0[i] = true;
                    found_index_1[j] = true;
                    found = true;
                    break;
                }
            }
        }

        if (not found)
        {
            // add a null angle with the same theta0, but null k
            found_index_0[i] = true;
            ang_params_1.append(std::tuple<int, int, int, double, double>(atom0, atom1, atom2, std::get<3>(ang0), 0.0));
        }
    }

    for (int j = 0; j < perturbed->ang_params.count(); ++j)
    {
        if (not found_index_1[j])
        {
            // need to add a bond missing in the reference state
            const auto &ang1 = perturbed->ang_params.at(j);

            int atom0 = std::get<0>(ang1);
            int atom1 = std::get<1>(ang1);
            int atom2 = std::get<2>(ang1);

            // add a null angle with the same theta0, but null k
            ang_params.append(std::tuple<int, int, int, double, double>(atom0, atom1, atom2, std::get<3>(ang1), 0.0));
            ang_params_1.append(ang1);

            found_index_1[j] = true;
        }
    }

    // all of found_index_0 and found_index_1 should be true...
    if (found_index_0.indexOf(false) != -1 or found_index_1.indexOf(false) != -1)
    {
        throw SireError::program_bug(QObject::tr(
                                         "Failed to align the angles!"),
                                     CODELOC);
    }

    perturbed->ang_params = ang_params_1;

    QVector<std::tuple<int, int, int, int, int, double, double>> dih_params_1;
    dih_params_1.reserve(dih_params.count());

    found_index_0 = QVector<bool>(dih_params.count(), false);
    found_index_1 = QVector<bool>(perturbed->dih_params.count(), false);

    for (int i = 0; i < dih_params.count(); ++i)
    {
        const auto &dih0 = dih_params.at(i);

        int atom0 = std::get<0>(dih0);
        int atom1 = std::get<1>(dih0);
        int atom2 = std::get<2>(dih0);
        int atom3 = std::get<3>(dih0);

        bool found = false;

        for (int j = 0; j < perturbed->dih_params.count(); ++j)
        {
            if (not found_index_1[j])
            {
                const auto &dih1 = perturbed->dih_params.at(j);

                // we need to match all of the atoms AND the periodicity
                if (std::get<0>(dih1) == atom0 and std::get<1>(dih1) == atom1 and
                    std::get<2>(dih1) == atom2 and std::get<3>(dih1) == atom3 and
                    std::get<4>(dih0) == std::get<4>(dih1))
                {
                    // we have found the matching torsion!
                    dih_params_1.append(dih1);
                    found_index_0[i] = true;
                    found_index_1[j] = true;
                    found = true;
                    break;
                }
            }
        }

        if (not found)
        {
            // add a null dihedral with the same periodicity and phase, but null k
            dih_params_1.append(std::tuple<int, int, int, int, int, double, double>(atom0, atom1, atom2, atom3, std::get<4>(dih0), std::get<5>(dih0), 0.0));
            found_index_0[i] = true;
        }
    }

    for (int j = 0; j < perturbed->dih_params.count(); ++j)
    {
        if (not found_index_1[j])
        {
            // need to add a dihedral missing in the reference state
            const auto &dih1 = perturbed->dih_params.at(j);

            int atom0 = std::get<0>(dih1);
            int atom1 = std::get<1>(dih1);
            int atom2 = std::get<2>(dih1);
            int atom3 = std::get<3>(dih1);

            // add a null dihedral with the same periodicity and phase, but null k
            dih_params.append(std::tuple<int, int, int, int, int, double, double>(atom0, atom1, atom2, atom3, std::get<4>(dih1), std::get<5>(dih1), 0.0));
            dih_params_1.append(dih1);
            found_index_1[j] = true;
        }
    }

    // all of found_index_0 and found_index_1 should be true...
    if (found_index_0.indexOf(false) != -1 or found_index_1.indexOf(false) != -1)
    {
        throw SireError::program_bug(QObject::tr(
                                         "Failed to align the dihedrals!"),
                                     CODELOC);
    }

    perturbed->dih_params = dih_params_1;

    // now align all of the exceptions - this should allow the bonding
    // to change during the perturbation
    QVector<std::tuple<int, int, double, double>> exception_params_1;
    exception_params_1.reserve(exception_params.count());

    found_index_0 = QVector<bool>(exception_params.count(), false);
    found_index_1 = QVector<bool>(perturbed->exception_params.count(), false);

    for (int i = 0; i < exception_params.count(); ++i)
    {
        const auto &ex0 = exception_params.at(i);

        int atom0 = std::get<0>(ex0);
        int atom1 = std::get<1>(ex0);

        bool found = false;

        for (int j = 0; j < perturbed->exception_params.count(); ++j)
        {
            const auto &ex1 = perturbed->exception_params.at(j);

            if (std::get<0>(ex1) == atom0 and std::get<1>(ex1) == atom1)
            {
                // we have found the matching exception!
                exception_params_1.append(ex1);
                found_index_0[i] = true;
                found_index_1[j] = true;
                found = true;
                break;
            }
        }

        if (not found)
        {
            // add a null exception with the same scale factors
            exception_params_1.append(std::tuple<int, int, double, double>(atom0, atom1, 1.0, 1.0));
            found_index_0[i] = true;
        }
    }

    for (int j = 0; j < perturbed->exception_params.count(); ++j)
    {
        if (not found_index_1[j])
        {
            // need to add an exception missing in the reference state
            const auto &ex1 = perturbed->exception_params.at(j);

            int atom0 = std::get<0>(ex1);
            int atom1 = std::get<1>(ex1);

            // add a null exception
            exception_params.append(std::tuple<int, int, double, double>(atom0, atom1, 1.0, 1.0));
            exception_params_1.append(ex1);
            found_index_1[j] = true;
        }
    }

    // all of found_index_0 and found_index_1 should be true...
    if (found_index_0.indexOf(false) != -1 or found_index_1.indexOf(false) != -1)
    {
        throw SireError::program_bug(QObject::tr(
                                         "Failed to align the exceptions!"),
                                     CODELOC);
    }

    perturbed->exception_params = exception_params_1;

    if (exception_params.count() != perturbed->exception_params.count())
    {
        throw SireError::program_bug(QObject::tr(
                                         "Different number of exceptions between the reference "
                                         "(%1) and perturbed (%2) states.")
                                         .arg(exception_params.count())
                                         .arg(perturbed->exception_params.count()),
                                     CODELOC);
    }
}

/** Internal function that builds all of the exceptions for all of the
    atoms in the molecule
*/
void OpenMMMolecule::buildExceptions(const Molecule &mol,
                                     QSet<qint64> &constrained_pairs,
                                     const PropertyMap &map)
{
    // we will build the complete exception list, and will not rely
    // on this list being built by an OpenMM forcefield. This is because
    // we need to allow this set to morph
    exception_params.clear();

    // save memory
    if (unbonded_atoms.isEmpty())
        unbonded_atoms = QSet<qint32>();

    const int nats = this->cljs.count();

    const auto &nbpairs = mol.property(map["intrascale"]).asA<CLJNBPairs>();

    const int ncgs = mol.nCutGroups();

    auto add_exception = [&](const AtomIdx &atom0, const AtomIdx &atom1,
                             double cscl, double ljscl)
    {
        if (cscl != 1 or ljscl != 1)
        {
            const int i = atom0.value();
            const int j = atom1.value();

            exception_params.append(std::make_tuple(i, j, cscl, ljscl));

            if (cscl == 0 and ljscl == 0)
            {
                // are any of these atoms unbonded? If so, then
                // we will need to add a constraint to hold them
                // in place
                if (unbonded_atoms.contains(i) or unbonded_atoms.contains(j))
                {
                    if (not constrained_pairs.contains(to_pair(i, j)))
                    {
                        const auto delta = coords[j] - coords[i];
                        const auto length = std::sqrt((delta[0] * delta[0]) +
                                                      (delta[1] * delta[1]) +
                                                      (delta[2] * delta[2]));
                        constraints.append(std::make_tuple(i, j, getSharedConstraintLength(length)));
                        constrained_pairs.insert(to_pair(i, j));
                    }
                }
            }
        }
    };

    // loop over all pairs of CutGroups and get the NB scale factors
    if (ncgs == 1)
    {
        if (nats == 1)
        {
            // nothing to do :-)
        }
        else if (nats == 2)
        {
            // two atoms must be excluded
            exception_params.append(std::make_tuple(0, 1, 0.0, 0.0));

            if (unbonded_atoms.contains(0) or unbonded_atoms.contains(1))
            {
                // these atoms are not bonded - we need to add a constraint
                // to keep them in place
                const auto delta = coords[1] - coords[0];
                const auto length = std::sqrt((delta[0] * delta[0]) +
                                              (delta[1] * delta[1]) +
                                              (delta[2] * delta[2]));
                constraints.append(std::make_tuple(0, 1, getSharedConstraintLength(length)));
                constrained_pairs.insert(to_pair(0, 1));
            }
        }
        else if (nats == 3)
        {
            // three atoms must be excluded
            exception_params.append(std::make_tuple(0, 1, 0.0, 0.0));
            exception_params.append(std::make_tuple(1, 2, 0.0, 0.0));
            exception_params.append(std::make_tuple(0, 2, 0.0, 0.0));

            if (unbonded_atoms.contains(0) or unbonded_atoms.contains(1) or unbonded_atoms.contains(2))
            {
                // these atoms are not bonded - we need to add a constraint
                // to keep them in place
                auto delta = coords[1] - coords[0];
                auto length = std::sqrt((delta[0] * delta[0]) +
                                        (delta[1] * delta[1]) +
                                        (delta[2] * delta[2]));
                constraints.append(std::make_tuple(0, 1, getSharedConstraintLength(length)));
                constrained_pairs.insert(to_pair(0, 1));

                delta = coords[2] - coords[0];
                length = std::sqrt((delta[0] * delta[0]) +
                                   (delta[1] * delta[1]) +
                                   (delta[2] * delta[2]));
                constraints.append(std::make_tuple(0, 2, getSharedConstraintLength(length)));
                constrained_pairs.insert(to_pair(0, 2));

                delta = coords[2] - coords[1];
                length = std::sqrt((delta[0] * delta[0]) +
                                   (delta[1] * delta[1]) +
                                   (delta[2] * delta[2]));
                constraints.append(std::make_tuple(1, 2, getSharedConstraintLength(length)));
                constrained_pairs.insert(to_pair(1, 2));
            }
        }
        else
        {
            // we only need to worry about ourselves
            const auto &cgpairs = nbpairs.get(CGIdx(0), CGIdx(0));

            if (cgpairs.isEmpty())
            {
                // all of the pairs have the same value (surprising!)
                const auto &cljscl = cgpairs.defaultValue();

                if (cljscl.coulomb() != 1 or cljscl.lj() != 1)
                {
                    for (int i = 0; i < nats - 1; ++i)
                    {
                        for (int j = i + 1; j < nats; ++j)
                        {
                            add_exception(AtomIdx(i), AtomIdx(j),
                                          cljscl.coulomb(), cljscl.lj());
                        }
                    }
                }
            }
            else
            {
                // the pairs have different values, so add these in
                for (int i = 0; i < nats - 1; ++i)
                {
                    for (int j = i + 1; j < nats; ++j)
                    {
                        const auto &scl = cgpairs.get(i, j);

                        if (scl.coulomb() != 1 or scl.lj() != 1)
                        {
                            add_exception(AtomIdx(i), AtomIdx(j),
                                          scl.coulomb(), scl.lj());
                        }
                    }
                }
            }
        }
    }
    else
    {
        const auto &molinfo = mol.data().info();

        // do all individual cutgroups
        for (int icg = 0; icg < ncgs; ++icg)
        {
            const int i_nats = molinfo.nAtoms(CGIdx(icg));

            const auto &cgpairs = nbpairs.get(CGIdx(icg), CGIdx(icg));

            if (cgpairs.isEmpty())
            {
                const auto &scl = cgpairs.defaultValue();

                if (scl.coulomb() != 1 or scl.lj() != 1)
                {
                    // all of the pairs are excluded for all atoms
                    for (int i = 0; i < i_nats - 1; ++i)
                    {
                        for (int j = i + 1; j < i_nats; ++j)
                        {
                            add_exception(molinfo.atomIdx(CGAtomIdx(CGIdx(icg), Index(i))),
                                          molinfo.atomIdx(CGAtomIdx(CGIdx(icg), Index(j))),
                                          scl.coulomb(), scl.lj());
                        }
                    }
                }
            }
            else
            {
                // all of the pairs are excluded for all atoms
                for (int i = 0; i < i_nats - 1; ++i)
                {
                    for (int j = i + 1; j < i_nats; ++j)
                    {
                        const auto &scl = cgpairs.get(i, j);

                        if (scl.coulomb() != 1 or scl.lj() != 1)
                        {
                            add_exception(molinfo.atomIdx(CGAtomIdx(CGIdx(icg), Index(i))),
                                          molinfo.atomIdx(CGAtomIdx(CGIdx(icg), Index(j))),
                                          scl.coulomb(), scl.lj());
                        }
                    }
                }
            }
        }

        // do all cutgroup pairs
        for (int icg = 0; icg < ncgs - 1; ++icg)
        {
            const int i_nats = molinfo.nAtoms(CGIdx(icg));

            for (int jcg = icg + 1; jcg < ncgs; ++jcg)
            {
                const int j_nats = molinfo.nAtoms(CGIdx(jcg));

                const auto &cgpairs = nbpairs.get(CGIdx(icg), CGIdx(jcg));

                if (cgpairs.isEmpty())
                {
                    const auto &scl = cgpairs.defaultValue();

                    if (scl.coulomb() != 1 or scl.lj() != 1)
                    {
                        // all of the pairs are excluded for all atoms
                        for (int i = 0; i < i_nats; ++i)
                        {
                            for (int j = 0; j < j_nats; ++j)
                            {
                                add_exception(molinfo.atomIdx(CGAtomIdx(CGIdx(icg), Index(i))),
                                              molinfo.atomIdx(CGAtomIdx(CGIdx(jcg), Index(j))),
                                              scl.coulomb(), scl.lj());
                            }
                        }
                    }
                }
                else
                {
                    // all of the pairs are excluded for all atoms
                    for (int i = 0; i < i_nats; ++i)
                    {
                        for (int j = 0; j < j_nats; ++j)
                        {
                            const auto &scl = cgpairs.get(i, j);

                            if (scl.coulomb() != 1 or scl.lj() != 1)
                            {
                                add_exception(molinfo.atomIdx(CGAtomIdx(CGIdx(icg), Index(i))),
                                              molinfo.atomIdx(CGAtomIdx(CGIdx(jcg), Index(j))),
                                              scl.coulomb(), scl.lj());
                            }
                        }
                    }
                }
            }
        }
    }
}

void OpenMMMolecule::copyInCoordsAndVelocities(OpenMM::Vec3 *c, OpenMM::Vec3 *v) const
{
    const auto *coords_data = coords.constData();

    if (c != 0 and coords_data != 0)
    {
        const int nats = coords.count();

        for (int i = 0; i < nats; ++i)
        {
            *c = coords_data[i];
            c += 1;
        }
    }

    const auto *vels_data = vels.constData();

    if (v != 0 and vels_data != 0)
    {
        const int nats = vels.count();

        for (int i = 0; i < nats; ++i)
        {
            *v = vels_data[i];
            v += 1;
        }
    }
}

/** Process this molecule to remove unnecessary parameters associated
 *  with field atoms. Field atoms are massless atoms that are not
 *  bonded to any other atom
 */
void OpenMMMolecule::processFieldAtoms()
{
    if (not this->hasFieldAtoms())
        return;

    if (this->isPerturbable())
        throw SireError::incompatible_error(QObject::tr(
                                                "We cannot (yet) support perturbable molecules that also "
                                                "contain field atoms. Please raise an issue if you would "
                                                "like us to write the code to do this."),
                                            CODELOC);

    if (this->isFieldMolecule())
    {
        qDebug() << "This is a field molecule";

        // very easy case - all atoms are field atoms
        field_points.reserve(field_atoms.count());

        for (int i = 0; i < this->atoms.count(); ++i)
        {
            const auto &clj = cljs.at(i);
            field_points.append(std::make_tuple(coords[i],
                                                std::get<0>(clj),
                                                std::get<1>(clj),
                                                std::get<2>(clj)));
        }

        coords.clear();
        vels.clear();
        masses.clear();
        light_atoms.clear();
        virtual_sites.clear();
        cljs.clear();
        exception_params.clear();
        bond_params.clear();
        ang_params.clear();
        dih_params.clear();
        constraints.clear();
        perturbed.reset();
        exception_idxs.clear();
        unbonded_atoms.clear();
        alphas.clear();
        to_ghost_idxs.clear();
        from_ghost_idxs.clear();

        return;
    }

    // more difficult case - a mixture of field and non-field atoms.
    // This is almost certainly a protein or similar

    // to start, the easiest is to create a bitmask to indicate
    // which atoms are field atoms - we can then go through the
    // atoms in sequence and then sort them into the appropriate
    // bins
    field_points.reserve(field_atoms.count());

    QVector<bool> is_field_atom(this->atoms.count(), false);

    for (const auto &atom : field_atoms)
    {
        is_field_atom[atom] = true;
    }
}

/** Return the alpha parameters of all atoms in atom order for
 *  this molecule
 */
QVector<double> OpenMMMolecule::getAlphas() const
{
    return this->alphas;
}

/** Return all of the atom charges in atom order for this molecule */
QVector<double> OpenMMMolecule::getCharges() const
{
    const int natoms = this->cljs.count();

    QVector<double> charges(natoms);

    auto charges_data = charges.data();
    const auto cljs_data = this->cljs.constData();

    for (int i = 0; i < natoms; ++i)
    {
        charges_data[i] = std::get<0>(cljs_data[i]);
    }

    return charges;
}

/** Return all of the LJ sigma parameters in atom order for this molecule */
QVector<double> OpenMMMolecule::getSigmas() const
{
    const int natoms = this->cljs.count();

    QVector<double> sigmas(natoms);

    auto sigmas_data = sigmas.data();
    const auto cljs_data = this->cljs.constData();

    for (int i = 0; i < natoms; ++i)
    {
        sigmas_data[i] = std::get<1>(cljs_data[i]);
    }

    return sigmas;
}

/** Return all of the LJ epsilon parameters in atom order for this molecule */
QVector<double> OpenMMMolecule::getEpsilons() const
{
    const int natoms = this->cljs.count();

    QVector<double> epsilons(natoms);

    auto epsilons_data = epsilons.data();
    const auto cljs_data = this->cljs.constData();

    for (int i = 0; i < natoms; ++i)
    {
        epsilons_data[i] = std::get<2>(cljs_data[i]);
    }

    return epsilons;
}

/** Return all of the bond k parameters in bond order for this molecule */
QVector<double> OpenMMMolecule::getBondKs() const
{
    const int nbonds = this->bond_params.count();

    QVector<double> bond_ks(nbonds);

    auto bond_ks_data = bond_ks.data();
    const auto bond_params_data = this->bond_params.constData();

    for (int i = 0; i < nbonds; ++i)
    {
        bond_ks_data[i] = std::get<3>(bond_params_data[i]);
    }

    return bond_ks;
}

/** Return all of the bond length parameters in bond order for this molecule */
QVector<double> OpenMMMolecule::getBondLengths() const
{
    const int nbonds = this->bond_params.count();

    QVector<double> bond_lengths(nbonds);

    auto bond_lengths_data = bond_lengths.data();
    const auto bond_params_data = this->bond_params.constData();

    for (int i = 0; i < nbonds; ++i)
    {
        bond_lengths_data[i] = std::get<2>(bond_params_data[i]);
    }

    return bond_lengths;
}

/** Return all of the angle K parameters in angle order for this molecule */
QVector<double> OpenMMMolecule::getAngleKs() const
{
    const int nangs = this->ang_params.count();

    QVector<double> ang_ks(nangs);

    auto ang_ks_data = ang_ks.data();
    const auto ang_params_data = this->ang_params.constData();

    for (int i = 0; i < nangs; ++i)
    {
        ang_ks_data[i] = std::get<4>(ang_params_data[i]);
    }

    return ang_ks;
}

/** Return all of the angle size parameters in angle order for this molecule */
QVector<double> OpenMMMolecule::getAngleSizes() const
{
    const int nangs = this->ang_params.count();

    QVector<double> ang_sizes(nangs);

    auto ang_sizes_data = ang_sizes.data();
    const auto ang_params_data = this->ang_params.constData();

    for (int i = 0; i < nangs; ++i)
    {
        ang_sizes_data[i] = std::get<3>(ang_params_data[i]);
    }

    return ang_sizes;
}

/** Return all of the dihedral torsion periodicities in dihedral order
 *  for this molecule
 */
QVector<int> OpenMMMolecule::getTorsionPeriodicities() const
{
    const int ndihs = this->dih_params.count();

    QVector<int> dih_periodicities(ndihs);

    auto dih_periodicities_data = dih_periodicities.data();
    const auto dih_params_data = this->dih_params.constData();

    for (int i = 0; i < ndihs; ++i)
    {
        dih_periodicities_data[i] = std::get<4>(dih_params_data[i]);
    }

    return dih_periodicities;
}

/** Return all of the torsion phase parameters for all of the dihedrals
 *  in dihedral order for this molecule
 */
QVector<double> OpenMMMolecule::getTorsionPhases() const
{
    const int ndihs = this->dih_params.count();

    QVector<double> dih_phases(ndihs);

    auto dih_phases_data = dih_phases.data();
    const auto dih_params_data = this->dih_params.constData();

    for (int i = 0; i < ndihs; ++i)
    {
        dih_phases_data[i] = std::get<5>(dih_params_data[i]);
    }

    return dih_phases;
}

/** Return all of the torsion K values for all of the dihedrals in
 *  dihedral order for this molecule
 */
QVector<double> OpenMMMolecule::getTorsionKs() const
{
    const int ndihs = this->dih_params.count();

    QVector<double> dih_ks(ndihs);

    auto dih_ks_data = dih_ks.data();
    const auto dih_params_data = this->dih_params.constData();

    for (int i = 0; i < ndihs; ++i)
    {
        dih_ks_data[i] = std::get<6>(dih_params_data[i]);
    }

    return dih_ks;
}

/** Return the atom indexes of the atoms in the exceptions, in
 *  exception order for this molecule
 */
QVector<std::pair<int, int>> OpenMMMolecule::getExceptionAtoms() const
{
    const int nexceptions = this->exception_params.count();

    QVector<std::pair<int, int>> exception_atoms(nexceptions);

    auto exception_atoms_data = exception_atoms.data();
    const auto exception_params_data = this->exception_params.constData();

    for (int i = 0; i < nexceptions; ++i)
    {
        exception_atoms_data[i] = std::make_pair(std::get<0>(exception_params_data[i]),
                                                 std::get<1>(exception_params_data[i]));
    }

    return exception_atoms;
}

/** Return all of the coulomb intramolecular scale factors (nbscl)
 *  in exception order for this molecule
 */
QVector<double> OpenMMMolecule::getChargeScales() const
{
    const int nexceptions = this->exception_params.count();

    QVector<double> charge_scales(nexceptions);

    auto charge_scales_data = charge_scales.data();
    const auto exception_params_data = this->exception_params.constData();

    for (int i = 0; i < nexceptions; ++i)
    {
        charge_scales_data[i] = std::get<2>(exception_params_data[i]);
    }

    return charge_scales;
}

/** Return all of the LJ intramolecular scale factors (nbscl)
 *  in exception order for this molecule
 */
QVector<double> OpenMMMolecule::getLJScales() const
{
    const int nexceptions = this->exception_params.count();

    QVector<double> lj_scales(nexceptions);

    auto lj_scales_data = lj_scales.data();
    const auto exception_params_data = this->exception_params.constData();

    for (int i = 0; i < nexceptions; ++i)
    {
        lj_scales_data[i] = std::get<3>(exception_params_data[i]);
    }

    return lj_scales;
}

////////
//////// Implementation of PerturbableOpenMMMolecule
////////

/** Null constructor */
PerturbableOpenMMMolecule::PerturbableOpenMMMolecule()
{
}

/** Construct from the passed OpenMMMolecule */
PerturbableOpenMMMolecule::PerturbableOpenMMMolecule(const OpenMMMolecule &mol)
{
    if (mol.perturbed.get() == 0)
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot construct a PerturbableOpenMMMolecule from an "
                                                "OpenMMMolecule that has not been perturbed!"),
                                            CODELOC);

    alpha0 = mol.getAlphas();
    alpha1 = mol.perturbed->getAlphas();

    chg0 = mol.getCharges();
    chg1 = mol.perturbed->getCharges();

    sig0 = mol.getSigmas();
    sig1 = mol.perturbed->getSigmas();

    eps0 = mol.getEpsilons();
    eps1 = mol.perturbed->getEpsilons();

    bond_k0 = mol.getBondKs();
    bond_k1 = mol.perturbed->getBondKs();

    bond_r0 = mol.getBondLengths();
    bond_r1 = mol.perturbed->getBondLengths();

    ang_k0 = mol.getAngleKs();
    ang_k1 = mol.perturbed->getAngleKs();

    ang_t0 = mol.getAngleSizes();
    ang_t1 = mol.perturbed->getAngleSizes();

    tors_k0 = mol.getTorsionKs();
    tors_k1 = mol.perturbed->getTorsionKs();

    tors_periodicity0 = mol.getTorsionPeriodicities();
    tors_periodicity1 = mol.perturbed->getTorsionPeriodicities();

    tors_phase0 = mol.getTorsionPhases();
    tors_phase1 = mol.perturbed->getTorsionPhases();

    charge_scl0 = mol.getChargeScales();
    charge_scl1 = mol.perturbed->getChargeScales();

    exception_atoms = mol.getExceptionAtoms();

    lj_scl0 = mol.getLJScales();
    lj_scl1 = mol.perturbed->getLJScales();

    to_ghost_idxs = mol.to_ghost_idxs;
    from_ghost_idxs = mol.from_ghost_idxs;
}

/** Copy constructor */
PerturbableOpenMMMolecule::PerturbableOpenMMMolecule(const PerturbableOpenMMMolecule &other)
    : alpha0(other.alpha0), alpha1(other.alpha1),
      chg0(other.chg0), chg1(other.chg1),
      sig0(other.sig0), sig1(other.sig1),
      eps0(other.eps0), eps1(other.eps1),
      bond_k0(other.bond_k0), bond_k1(other.bond_k1),
      bond_r0(other.bond_r0), bond_r1(other.bond_r1),
      ang_k0(other.ang_k0), ang_k1(other.ang_k1),
      ang_t0(other.ang_t0), ang_t1(other.ang_t1),
      tors_k0(other.tors_k0), tors_k1(other.tors_k1),
      tors_periodicity0(other.tors_periodicity0),
      tors_periodicity1(other.tors_periodicity1),
      tors_phase0(other.tors_phase0), tors_phase1(other.tors_phase1),
      charge_scl0(other.charge_scl0), charge_scl1(other.charge_scl1),
      lj_scl0(other.lj_scl0), lj_scl1(other.lj_scl1),
      to_ghost_idxs(other.to_ghost_idxs), from_ghost_idxs(other.from_ghost_idxs),
      exception_atoms(other.exception_atoms), exception_idxs(other.exception_idxs)
{
}

/** Destructor */
PerturbableOpenMMMolecule::~PerturbableOpenMMMolecule()
{
}

/** Comparison operator */
bool PerturbableOpenMMMolecule::operator==(const PerturbableOpenMMMolecule &other) const
{
    return alpha0 == alpha1 and chg0 == chg1 and sig0 == sig1 and eps0 == eps1 and
           bond_k0 == bond_k1 and bond_r0 == bond_r1 and ang_k0 == ang_k1 and
           ang_t0 == ang_t1 and tors_k0 == tors_k1 and tors_periodicity0 == tors_periodicity1 and
           tors_phase0 == tors_phase1 and charge_scl0 == charge_scl1 and lj_scl0 == lj_scl1 and
           to_ghost_idxs == other.to_ghost_idxs and from_ghost_idxs == other.from_ghost_idxs and
           exception_atoms == other.exception_atoms and exception_idxs == other.exception_idxs;
}

/** Comparison operator */
bool PerturbableOpenMMMolecule::operator!=(const PerturbableOpenMMMolecule &other) const
{
    return not this->operator==(other);
}

/** Return the alpha parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getAlphas0() const
{
    return this->alpha0;
}

/** Return the alpha parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getAlphas1() const
{
    return this->alpha1;
}

/** Return the atom charges of the reference state */
QVector<double> PerturbableOpenMMMolecule::getCharges0() const
{
    return this->chg0;
}

/** Return the atom charges of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getCharges1() const
{
    return this->chg1;
}

/** Return the LJ sigma parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getSigmas0() const
{
    return this->sig0;
}

/** Return the LJ sigma parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getSigmas1() const
{
    return this->sig1;
}

/** Return the LJ epsilon parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getEpsilons0() const
{
    return this->eps0;
}

/** Return the LJ epsilon parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getEpsilons1() const
{
    return this->eps1;
}

/** Return the bond k parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getBondKs0() const
{
    return this->bond_k0;
}

/** Return the bond k parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getBondKs1() const
{
    return this->bond_k1;
}

/** Return the bond length parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getBondLengths0() const
{
    return this->bond_r0;
}

/** Return the bond length parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getBondLengths1() const
{
    return this->bond_r1;
}

/** Return the angle k parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getAngleKs0() const
{
    return this->ang_k0;
}

/** Return the angle k parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getAngleKs1() const
{
    return this->ang_k1;
}

/** Return the angle size parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getAngleSizes0() const
{
    return this->ang_t0;
}

/** Return the angle size parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getAngleSizes1() const
{
    return this->ang_t1;
}

/** Return the torsion k parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getTorsionKs0() const
{
    return this->tors_k0;
}

/** Return the torsion k parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getTorsionKs1() const
{
    return this->tors_k1;
}

/** Return the torsion periodicity parameters of the reference state */
QVector<int> PerturbableOpenMMMolecule::getTorsionPeriodicities0() const
{
    return this->tors_periodicity0;
}

/** Return the torsion periodicity parameters of the perturbed state */
QVector<int> PerturbableOpenMMMolecule::getTorsionPeriodicities1() const
{
    return this->tors_periodicity1;
}

/** Return the torsion phase parameters of the reference state */
QVector<double> PerturbableOpenMMMolecule::getTorsionPhases0() const
{
    return this->tors_phase0;
}

/** Return the torsion phase parameters of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getTorsionPhases1() const
{
    return this->tors_phase1;
}

/** Return the coulomb intramolecular scale factors of the reference state */
QVector<double> PerturbableOpenMMMolecule::getChargeScales0() const
{
    return this->charge_scl0;
}

/** Return the coulomb intramolecular scale factors of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getChargeScales1() const
{
    return this->charge_scl1;
}

/** Return the LJ intramolecular scale factors of the reference state */
QVector<double> PerturbableOpenMMMolecule::getLJScales0() const
{
    return this->lj_scl0;
}

/** Return the LJ intramolecular scale factors of the perturbed state */
QVector<double> PerturbableOpenMMMolecule::getLJScales1() const
{
    return this->lj_scl1;
}

/** Return the indexes of the atoms that are to be ghosted in the
 *  perturbed state
 */
QSet<int> PerturbableOpenMMMolecule::getToGhostIdxs() const
{
    return this->to_ghost_idxs;
}

/** Return the indexes of the atoms that were ghosts in the
 *  reference state
 */
QSet<int> PerturbableOpenMMMolecule::getFromGhostIdxs() const
{
    return this->from_ghost_idxs;
}

/** Return true if the atom is a ghost atom in the
 *  referenece or perturbed states */
bool PerturbableOpenMMMolecule::isGhostAtom(int atom) const
{
    return from_ghost_idxs.contains(atom) or to_ghost_idxs.contains(atom);
}

/** Return the indices of the atoms in the exceptions */
QVector<std::pair<int, int>> PerturbableOpenMMMolecule::getExceptionAtoms() const
{
    return this->exception_atoms;
}

/** Return the global indexes of the exceptions in the non-bonded and
 *  ghost-14 forces
 */
QVector<std::pair<int, int>> PerturbableOpenMMMolecule::getExceptionIndicies(const QString &name) const
{
    return this->exception_idxs.value(name);
}

/** Set the global indexes of the exceptions in the non-bonded and
 *  ghost-14 forces
 */
void PerturbableOpenMMMolecule::setExceptionIndicies(const QString &name,
                                                     const QVector<std::pair<int, int>> &exception_idxs)
{
    if (exception_idxs.count() != this->exception_atoms.count())
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of exception indicies (%1) does not match the number of exceptions (%2)")
                                                .arg(exception_idxs.count())
                                                .arg(this->exception_atoms.count()),
                                            CODELOC);

    this->exception_idxs.insert(name, exception_idxs);
}
