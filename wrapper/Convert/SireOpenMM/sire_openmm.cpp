
#include "sire_openmm.h"

#include <OpenMM.h>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

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

#include <QDebug>

using SireBase::PropertyMap;
using SireMol::Molecule;
using SireMol::SelectorMol;

namespace SireOpenMM
{
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

    /** Internal class used to hold all of the extracted information
     *  of an OpenMM Molecule
     */
    class OpenMMMolecule
    {
    public:
        OpenMMMolecule()
        {
        }

        OpenMMMolecule(const Molecule &mol, const PropertyMap &map)
        {
            const auto &moldata = mol.data();
            molinfo = mol.info();

            const int nats = molinfo.nAtoms();

            if (nats <= 0)
            {
                return;
            }

            coords = QVector<OpenMM::Vec3>(nats, OpenMM::Vec3(0, 0, 0));
            vels = QVector<OpenMM::Vec3>(nats, OpenMM::Vec3(0, 0, 0));
            masses = QVector<double>(nats, 0.0);
            charges = QVector<double>(nats, 0.0);
            sigmas = QVector<double>(nats, 0.0);
            epsilons = QVector<double>(nats, 0.0);
            indexes = QVector<int>(nats, -1);

            params = SireMM::AmberParams(mol, map);

            auto idx_to_cgatomidx = QVector<SireMol::CGAtomIdx>(nats);
            auto idx_to_cgatomidx_data = idx_to_cgatomidx.data();

            for (int i = 0; i < nats; ++i)
            {
                idx_to_cgatomidx_data[i] = molinfo.cgAtomIdx(SireMol::AtomIdx(i));
            }

            const auto coords_prop = map["coordinates"];

            if (moldata.hasProperty(coords_prop))
            {
                const auto &coords = moldata.property(coords_prop).asA<SireMol::AtomCoords>();
                auto mol_coords = this->coords.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol_coords[i] = to_vec3(coords.at(idx_to_cgatomidx_data[i]));
                }
            }

            const auto vels_prop = map["velocities"];

            if (moldata.hasProperty(vels_prop))
            {
                const auto &vels = moldata.property(vels_prop).asA<SireMol::AtomVelocities>();
                auto mol_vels = this->vels.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol_vels[i] = to_vec3(vels.at(idx_to_cgatomidx_data[i]));
                }
            }

            const auto params_masses = params.masses();
            auto mol_masses = this->masses.data();

            for (int i = 0; i < nats; ++i)
            {
                mol_masses[i] = params_masses.at(idx_to_cgatomidx_data[i]).to(SireUnits::g_per_mol);

                if (mol_masses[i] < 0.05)
                {
                    mol_masses[i] = 0.0;
                }

                if (mol_masses[i] < 0.5)
                {
                    this->virtual_sites.append(i);
                }
                else if (mol_masses[i] < 2.5)
                {
                    this->light_atoms.append(i);
                }
            }

            const auto params_charges = params.charges();
            auto mol_charges = this->charges.data();

            for (int i = 0; i < nats; ++i)
            {
                mol_charges[i] = params_charges.at(idx_to_cgatomidx_data[i]).to(SireUnits::mod_electron);
            }

            const auto params_ljs = params.ljs();
            auto mol_sigmas = this->sigmas.data();
            auto mol_epsilons = this->epsilons.data();

            for (int i = 0; i < nats; ++i)
            {
                const auto &lj = params_ljs.at(idx_to_cgatomidx_data[i]);

                mol_sigmas[i] = lj.sigma().to(SireUnits::nanometer);
                mol_epsilons[i] = lj.epsilon().to(SireUnits::kJ_per_mol);
            }
        }

        ~OpenMMMolecule()
        {
        }

        SireMol::MoleculeInfo molinfo;
        SireMM::AmberParams params;

        QVector<OpenMM::Vec3> coords;
        QVector<OpenMM::Vec3> vels;
        QVector<double> masses;
        QVector<double> charges;
        QVector<double> sigmas;
        QVector<double> epsilons;
        QVector<int> indexes;
        QList<int> light_atoms;
        QList<int> virtual_sites;
    };

    bool use_parallel(int n, const SireBase::PropertyMap &map)
    {
        if (n <= 16)
            return false;

        if (map["parallel"].hasValue())
        {
            return map["parallel"].value().asABoolean();
        }

        return true;
    }

    SelectorMol openmm_to_sire(const OpenMM::System &mols,
                               const PropertyMap &map)
    {
        throw SireError::incomplete_code(QObject::tr(
                                             "Still need to write openmm_to_sire"),
                                         CODELOC);

        return SelectorMol();
    }

    void sire_to_openmm(OpenMM::System &system,
                        const SelectorMol &mols,
                        const PropertyMap &map)
    {
        // we can assume that an empty system has been passed to us

        // we will get parameters from the map (e.g. cutoffs etc)

        // extract the data from all of the molecules
        const int nmols = mols.count();
        QVector<OpenMMMolecule> openmm_mols(nmols);
        auto openmm_mols_data = openmm_mols.data();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, mols.count()), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    openmm_mols_data[i] = OpenMMMolecule(mols[i], map);
                } });
        }
        else
        {
            for (int i = 0; i < nmols; ++i)
            {
                openmm_mols_data[i] = OpenMMMolecule(mols[i], map);
            }
        }

        // create all the forcefields - ownership is taken by 'system'
        OpenMM::NonbondedForce *cljff = new OpenMM::NonbondedForce();
        system.addForce(cljff);

        OpenMM::HarmonicBondForce *bondff = new OpenMM::HarmonicBondForce();
        system.addForce(bondff);

        OpenMM::HarmonicAngleForce *angff = new OpenMM::HarmonicAngleForce();
        system.addForce(angff);

        OpenMM::PeriodicTorsionForce *dihff = new OpenMM::PeriodicTorsionForce();
        system.addForce(dihff);

        // set the default 1-4 scaling factors. Eventually we may
        // want to loop over the molecules and find the most popular
        // scaling factors. For now, we are assuming everything is amber

        double coul_14_scale = 1.0 / 1.2; // the default amber value
        double lj_14_scale = 1.0 / 1.2;   // the default amber value

        // now place this data into the OpenMM system - this looks like
        // it has to be done in serial?
        int openmm_index = 0;
        std::vector<std::pair<int, int>> bond_pairs;
        QVector<std::tuple<int, int, double, double, double>> custom_pairs;

        // assume arithmetic combining rules...
        auto combining_rules = SireMM::LJParameter::ARITHMETIC;

        // (will deal with restraints, light atoms and virtual sites later)
        for (int i = 0; i < nmols; ++i)
        {
            auto &mol = openmm_mols_data[i];

            auto mol_indexes = mol.indexes.data();
            auto mol_masses = mol.masses.constData();
            auto mol_charges = mol.charges.constData();
            auto mol_sigmas = mol.sigmas.constData();
            auto mol_epsilons = mol.epsilons.constData();

            // first the atom parameters
            for (int j = 0; j < mol.molinfo.nAtoms(); ++j)
            {
                mol_indexes[j] = openmm_index;
                openmm_index += 1;
                system.addParticle(mol_masses[j]);
                cljff->addParticle(mol_charges[j], mol_sigmas[j], mol_epsilons[j]);
            }

            const auto &params = mol.params;

            // now the bonds
            const double bond_k_to_openmm = (SireUnits::kcal_per_mol / (SireUnits::angstrom * SireUnits::angstrom)).to(SireUnits::kJ_per_mol / (SireUnits::nanometer * SireUnits::nanometer));
            const double bond_r0_to_openmm = 2.0 * (SireUnits::angstrom).to(SireUnits::nanometer);

            for (auto it = params.bonds().constBegin();
                 it != params.bonds().constEnd();
                 ++it)
            {
                const auto bondid = it.key().map(mol.molinfo);
                const auto &bondparam = it.value().first;
                const auto includes_h = it.value().second;

                const int atom0 = bondid.get<0>().value();
                const int atom1 = bondid.get<1>().value();

                const int atom0_index = mol_indexes[atom0];
                const int atom1_index = mol_indexes[atom1];

                const double k = bondparam.k() * bond_k_to_openmm;
                const double r0 = bondparam.r0() * bond_r0_to_openmm;

                // add this as a standard bond - will deal with this
                // being a bond with hydrogen (or a light atom) later
                bondff->addBond(atom0_index, atom1_index, r0, k);

                bond_pairs.push_back(std::make_pair(atom0_index, atom1_index));
            }

            // now the angles
            const double angle_k_to_openmm = 2.0 * (SireUnits::kcal_per_mol).to(SireUnits::kJ_per_mol);

            for (auto it = params.angles().constBegin();
                 it != params.angles().constEnd();
                 ++it)
            {
                const auto angid = it.key().map(mol.molinfo);
                const auto &angparam = it.value().first;
                const auto includes_h = it.value().second;

                const int atom0 = angid.get<0>().value();
                const int atom1 = angid.get<1>().value();
                const int atom2 = angid.get<2>().value();

                const int atom0_index = mol_indexes[atom0];
                const int atom1_index = mol_indexes[atom1];
                const int atom2_index = mol_indexes[atom2];

                const double k = angparam.k() * angle_k_to_openmm;
                const double theta0 = angparam.theta0(); // already in radians

                // add this as a standard angle - will deal with H angle
                // constraints later
                angff->addAngle(atom0_index, atom1_index, atom2_index,
                                theta0, k);
            }

            // now the dihedrals
            const double dihedral_k_to_openmm = (SireUnits::kcal_per_mol).to(SireUnits::kJ_per_mol);

            for (auto it = params.dihedrals().constBegin();
                 it != params.dihedrals().constEnd();
                 ++it)
            {
                const auto dihid = it.key().map(mol.molinfo);
                const auto &dihparam = it.value().first;
                const auto includes_h = it.value().second;

                const int atom0 = dihid.get<0>().value();
                const int atom1 = dihid.get<1>().value();
                const int atom2 = dihid.get<2>().value();
                const int atom3 = dihid.get<3>().value();

                const int atom0_index = mol_indexes[atom0];
                const int atom1_index = mol_indexes[atom1];
                const int atom2_index = mol_indexes[atom2];
                const int atom3_index = mol_indexes[atom3];

                for (const auto &term : dihparam.terms())
                {
                    const double v = term.k() * dihedral_k_to_openmm;
                    const double phase = term.phase();             // already in radians
                    const double periodicity = term.periodicity(); // already in correct units

                    dihff->addTorsion(atom0_index, atom1_index,
                                      atom2_index, atom3_index,
                                      periodicity, phase, v);
                }
            }

            // now the impropers
            for (auto it = params.impropers().constBegin();
                 it != params.impropers().constEnd();
                 ++it)
            {
                const auto impid = it.key().map(mol.molinfo);
                const auto &impparam = it.value().first;
                const auto includes_h = it.value().second;

                const int atom0 = impid.get<0>().value();
                const int atom1 = impid.get<0>().value();
                const int atom2 = impid.get<0>().value();
                const int atom3 = impid.get<0>().value();

                const int atom0_index = mol_indexes[atom0];
                const int atom1_index = mol_indexes[atom1];
                const int atom2_index = mol_indexes[atom2];
                const int atom3_index = mol_indexes[atom3];

                for (const auto &term : impparam.terms())
                {
                    const double v = term.k() * dihedral_k_to_openmm;
                    const double phase = term.phase();             // already in radians
                    const double periodicity = term.periodicity(); // already in correct units

                    dihff->addTorsion(atom0_index, atom1_index,
                                      atom2_index, atom3_index,
                                      periodicity, phase, v);
                }
            }

            // now find all of the 1-4 pairs that don't have these values
            for (auto it = params.nb14s().constBegin();
                 it != params.nb14s().constEnd();
                 ++it)
            {
                const double cscl = it.value().cscl();
                const double ljscl = it.value().ljscl();

                if (std::abs(cscl - coul_14_scale) > 0.0001 or
                    std::abs(ljscl - lj_14_scale) > 0.0001)
                {
                    const auto nbid = it.key().map(mol.molinfo);

                    const int atom0 = nbid.get<0>().value();
                    const int atom1 = nbid.get<1>().value();

                    const int atom0_index = mol_indexes[atom0];
                    const int atom1_index = mol_indexes[atom1];

                    const auto cgatomidx0 = mol.molinfo.cgAtomIdx(nbid.get<0>());
                    const auto cgatomidx1 = mol.molinfo.cgAtomIdx(nbid.get<1>());

                    const auto chg0 = params.charges().at(cgatomidx0).to(SireUnits::mod_electron);
                    const auto chg1 = params.charges().at(cgatomidx1).to(SireUnits::mod_electron);

                    const auto lj0 = params.ljs().at(cgatomidx0);
                    const auto lj1 = params.ljs().at(cgatomidx1);

                    double charge_pair = cscl * chg0 * chg1;

                    // Amber, so assume arithmetic combining rules - this should really
                    // be got from the molecule / forcefield or something...
                    auto lj_pair = lj0.combine(lj1, combining_rules);

                    double sig_pair = lj_pair.sigma().to(SireUnits::nanometer);
                    double eps_pair = lj_pair.epsilon().to(SireUnits::kJ_per_mol);

                    custom_pairs.append(std::make_tuple(
                        atom0_index, atom1_index, charge_pair, sig_pair, eps_pair));
                }
            }
        }

        // add exclusions based on the bonding of this molecule
        cljff->createExceptionsFromBonds(bond_pairs, coul_14_scale, lj_14_scale);

        for (const auto &p : custom_pairs)
        {
            cljff->addException(std::get<0>(p), std::get<1>(p),
                                std::get<2>(p), std::get<3>(p),
                                std::get<4>(p), true);
        }
    }

} // end of namespace SireOpenMM
