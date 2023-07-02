
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

#include "SireMaths/vector.h"

#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "tostring.h"

#include <QSet>

#include <QDebug>

using namespace SireOpenMM;
using namespace SireMM;
using namespace SireMol;
using namespace SireBase;

OpenMMMolecule::OpenMMMolecule()
{
}

OpenMMMolecule::OpenMMMolecule(const Molecule &mol,
                               const PropertyMap &map)
{
    molinfo = mol.info();

    if (molinfo.nAtoms() == 0)
    {
        // nothing to extract
        molinfo = MoleculeInfo();
        return;
    }

    ffinfo = mol.property(map["forcefield"]).asA<MMDetail>();

    if (map.specified("constraint"))
    {
        const auto c = map["constraint"].source().toLower().simplified();

        if (c == "none")
        {
            constraint_type = CONSTRAIN_NONE;
        }
        else if (c == "h-bonds")
        {
            constraint_type = CONSTRAIN_HBONDS;
        }
        else if (c == "bonds")
        {
            constraint_type = CONSTRAIN_BONDS;
        }
        else if (c == "h-bonds-h-angles")
        {
            constraint_type = CONSTRAIN_HBONDS | CONSTRAIN_ANGLES;
        }
        else if (c == "bonds-h-angles")
        {
            constraint_type = CONSTRAIN_BONDS | CONSTRAIN_ANGLES;
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

    if (ffinfo.isAmberStyle())
    {
        this->constructFromAmber(mol, map);
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
}

OpenMMMolecule::~OpenMMMolecule()
{
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

void OpenMMMolecule::constructFromAmber(const Molecule &mol, const PropertyMap &map)
{
    const auto &moldata = mol.data();
    const int nats = molinfo.nAtoms();

    if (nats <= 0)
    {
        return;
    }

    // extract the parameters in amber format - this should work,
    // as the 'forcefield' property has promised that this is
    // an amber-style molecule
    const auto params = SireMM::AmberParams(mol, map);

    // look up the CGAtomIdx of each atom - this is because we
    // will use AtomIdx for the ordering and atom identifiers
    auto idx_to_cgatomidx = QVector<SireMol::CGAtomIdx>(nats);
    auto idx_to_cgatomidx_data = idx_to_cgatomidx.data();

    for (int i = 0; i < nats; ++i)
    {
        idx_to_cgatomidx_data[i] = molinfo.cgAtomIdx(SireMol::AtomIdx(i));
    }

    // extract the coordinates and convert to OpenMM units
    const auto coords_prop = map["coordinates"];
    coords = QVector<OpenMM::Vec3>(nats, OpenMM::Vec3(0, 0, 0));

    if (moldata.hasProperty(coords_prop))
    {
        const auto &c = moldata.property(coords_prop).asA<SireMol::AtomCoords>();
        auto coords_data = coords.data();

        for (int i = 0; i < nats; ++i)
        {
            coords_data[i] = to_vec3(c.at(idx_to_cgatomidx_data[i]));
        }
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
    masses = QVector<double>(nats, 0.0);

    auto masses_data = masses.data();

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

    // extract the charges and LJ parameters and convert to OpenMM units
    const auto params_charges = params.charges();
    const auto params_ljs = params.ljs();

    cljs = QVector<std::tuple<double, double, double>>(nats, std::make_tuple(0.0, 0.0, 0.0));
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

    bond_pairs.clear();
    bond_params.clear();
    constraints.clear();

    auto bond_param = bond_params.data();
    auto bond_pair = bond_pairs.data();

    // now the bonds
    const double bond_k_to_openmm = 2.0 * (SireUnits::kcal_per_mol / (SireUnits::angstrom * SireUnits::angstrom)).to(SireUnits::kJ_per_mol / (SireUnits::nanometer * SireUnits::nanometer));
    const double bond_r0_to_openmm = SireUnits::angstrom.to(SireUnits::nanometer);

    QSet<qint64> constrained_pairs;

    for (auto it = params.bonds().constBegin();
         it != params.bonds().constEnd();
         ++it)
    {
        const auto bondid = it.key().map(molinfo);
        const auto &bondparam = it.value().first;
        const auto includes_h = it.value().second;

        const int atom0 = bondid.get<0>().value();
        const int atom1 = bondid.get<1>().value();

        const double k = bondparam.k() * bond_k_to_openmm;
        const double r0 = bondparam.r0() * bond_r0_to_openmm;

        const bool has_light_atom = includes_h or (masses_data[atom0] < 2.5 or masses_data[atom1] < 2.5);
        const bool has_massless_atom = masses_data[atom0] < 0.5 or masses_data[atom1] < 0.5;

        bond_pairs.append(std::make_pair(atom0, atom1));

        if ((not has_massless_atom) and ((constraint_type & CONSTRAIN_BONDS) or (has_light_atom and (constraint_type & CONSTRAIN_HBONDS))))
        {
            constraints.append(std::make_tuple(atom0, atom1, r0));
            constrained_pairs.insert(to_pair(atom0, atom1));
        }
        else
        {
            bond_params.append(std::make_tuple(atom0, atom1, r0, k));
        }
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
        const auto includes_h = it.value().second;

        const int atom0 = angid.get<0>().value();
        const int atom1 = angid.get<1>().value();
        const int atom2 = angid.get<2>().value();

        const double k = angparam.k() * angle_k_to_openmm;
        const double theta0 = angparam.theta0(); // already in radians

        const bool is_h_x_h = includes_h and (masses_data[atom0] < 2.5 and masses_data[atom2] < 2.5);

        if (not constrained_pairs.contains(to_pair(atom0, atom2)))
        {
            // only include the angle X-y-Z if X-Z are not already constrained
            if ((constraint_type & CONSTRAIN_ANGLES) and is_h_x_h)
            {
                const auto delta = coords[atom2] - coords[atom0];
                const auto length = std::sqrt((delta[0] * delta[0]) +
                                              (delta[1] * delta[1]) +
                                              (delta[2] * delta[2]));
                constraints.append(std::make_tuple(atom0, atom2, length));
                constrained_pairs.insert(to_pair(atom0, atom2));
            }
            else
            {
                ang_params.append(std::make_tuple(atom0, atom1, atom2,
                                                  theta0, k));
            }
        }
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
        const auto includes_h = it.value().second;

        const int atom0 = dihid.get<0>().value();
        const int atom1 = dihid.get<1>().value();
        const int atom2 = dihid.get<2>().value();
        const int atom3 = dihid.get<3>().value();

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
        const auto includes_h = it.value().second;

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

    // now find all of the 1-4 pairs that don't have standard 1-4 values
    const double coul_14_scale = ffinfo.electrostatic14ScaleFactor();
    const double lj_14_scale = ffinfo.vdw14ScaleFactor();

    bool is_arithmetic = ffinfo.usesArithmeticCombiningRules();
    bool is_geometric = ffinfo.usesGeometricCombiningRules();

    for (auto it = params.nb14s().constBegin();
         it != params.nb14s().constEnd();
         ++it)
    {
        const double cscl = it.value().cscl();
        const double ljscl = it.value().ljscl();

        if (std::abs(cscl - coul_14_scale) > 0.0001 or
            std::abs(ljscl - lj_14_scale) > 0.0001)
        {
            const auto nbid = it.key().map(molinfo);

            const int atom0 = nbid.get<0>().value();
            const int atom1 = nbid.get<1>().value();

            const auto cgatomidx0 = idx_to_cgatomidx_data[atom0];
            const auto cgatomidx1 = idx_to_cgatomidx_data[atom1];

            const auto &clj0 = cljs[atom0];
            const auto &clj1 = cljs[atom1];

            const auto chg0 = std::get<0>(clj0);
            const auto chg1 = std::get<0>(clj1);

            const auto sig0 = std::get<1>(clj0);
            const auto sig1 = std::get<1>(clj1);

            const auto eps0 = std::get<2>(clj0);
            const auto eps1 = std::get<2>(clj1);

            double charge_pair = cscl * chg0 * chg1;
            double sig_pair = 0.5 * ljscl * (sig0 + sig1);
            double eps_pair = ljscl * std::sqrt(eps0 * eps1);

            if (is_geometric)
                sig_pair = ljscl * std::sqrt(sig0 * sig1);
            else if (not is_arithmetic)
                throw SireError::unsupported(QObject::tr(
                                                 "We only support arithmetic or geometric combining rules. "
                                                 "We cannot support the rules in %1.")
                                                 .arg(ffinfo.toString()),
                                             CODELOC);

            custom_pairs.append(std::make_tuple(
                atom0, atom1, charge_pair, sig_pair, eps_pair));
        }
    }

    // finally, add in all of the excluded atoms
    excl_pairs = ExcludedPairs(mol, map);

    if (excl_pairs.count() > 0)
    {
        for (int i = 0; i < excl_pairs.count(); ++i)
        {
            auto pair = excl_pairs[i];

            custom_pairs.append(std::make_tuple(
                std::get<0>(pair).value(), std::get<1>(pair).value(), 0.0, 0.0, 0.0));
        }

        // and finally (finally!) find any atoms that are not bonded to
        // anything else and make sure that they are constrained. These
        // atoms will be excluded atoms (by definition) so just look
        // through those
        const auto &connectivity = mol.property(map["connectivity"]).asA<Connectivity>();

        for (int i = 0; i < excl_pairs.count(); ++i)
        {
            auto pair = excl_pairs[i];

            const int atom0 = std::get<0>(pair).value();
            const int atom1 = std::get<1>(pair).value();

            if (not constrained_pairs.contains(to_pair(atom0, atom1)))
            {
                if (connectivity.nConnections(std::get<0>(pair)) == 0 or
                    connectivity.nConnections(std::get<1>(pair)) == 0)
                {
                    const auto delta = coords[atom1] - coords[atom0];
                    const auto length = std::sqrt((delta[0] * delta[0]) +
                                                  (delta[1] * delta[1]) +
                                                  (delta[2] * delta[2]));
                    constraints.append(std::make_tuple(atom0, atom1, length));
                    constrained_pairs.insert(to_pair(atom0, atom1));
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
