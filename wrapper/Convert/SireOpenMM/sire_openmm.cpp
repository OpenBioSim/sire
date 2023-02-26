
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
    /** Internal class used to hold all of the extracted information
     *  of an OpenMM Molecule
     */
    class OpenMMMolecule
    {
    public:
        OpenMMMolecule() : nats(0)
        {
        }

        OpenMMMolecule(int num_atoms)
        {
            if (num_atoms <= 0)
            {
                nats = 0;
                return;
            }

            nats = num_atoms;

            coords = QVector<OpenMM::Vec3>(nats);
            vels = QVector<OpenMM::Vec3>(nats);
            masses = QVector<double>(nats);
            charges = QVector<double>(nats);
            sigmas = QVector<double>(nats);
            epsilons = QVector<double>(nats);
            indexes = QVector<int>(nats);
        }

        ~OpenMMMolecule()
        {
        }

        int nats;
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

    void sire_to_openmm(OpenMM::System &system,
                        const SelectorMol &mols,
                        const PropertyMap &map)
    {
        // we can assume that an empty system has been passed to us

        // we will get parameters from the map (e.g. cutoffs etc)

        // get all the names of the atom properties
        const auto chg_prop = map["charge"];
        const auto lj_prop = map["LJ"];
        const auto elem_prop = map["element"];
        const auto mass_prop = map["mass"];
        const auto coords_prop = map["coordinates"];
        const auto vels_prop = map["velocities"];
        const auto conn_prop = map["connectivity"];

        // extract all of the properties into an intermediate format
        auto to_openmm_mol = [&](const Molecule &molecule)
        {
            const auto &moldata = molecule.data();
            const auto &molinfo = moldata.info();

            const int nats = molinfo.nAtoms();

            QVector<SireMol::CGAtomIdx> idx_to_cgatomidx(nats);
            auto idx_to_cgatomidx_data = idx_to_cgatomidx.data();

            for (int i = 0; i < nats; ++i)
            {
                idx_to_cgatomidx_data[i] = molinfo.cgAtomIdx(SireMol::AtomIdx(i));
            }

            OpenMMMolecule mol(nats);

            try
            {
                const auto &coords = moldata.property(coords_prop).asA<SireMol::AtomCoords>();
                auto mol_coords = mol.coords.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol_coords[i] = to_vec3(coords.at(idx_to_cgatomidx_data[i]));
                }
            }
            catch (...)
            {
                auto mol_coords = mol.coords.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol_coords[i] = OpenMM::Vec3(0, 0, 0);
                }
            }

            try
            {
                const auto &vels = moldata.property(vels_prop).asA<SireMol::AtomVelocities>();
                auto mol_vels = mol.vels.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol_vels[i] = to_vec3(vels.at(idx_to_cgatomidx_data[i]));
                }
            }
            catch (...)
            {
                auto mol_vels = mol.vels.data();
                for (int i = 0; i < nats; ++i)
                {
                    mol_vels[i] = OpenMM::Vec3(0, 0, 0);
                }
            }

            auto mol_masses = mol.masses.data();

            try
            {
                const auto &masses = moldata.property(mass_prop).asA<SireMol::AtomMasses>();

                for (int i = 0; i < nats; ++i)
                {
                    mol_masses[i] = masses.at(idx_to_cgatomidx_data[i]).to(SireUnits::g_per_mol);
                }
            }
            catch (...)
            {
                try
                {
                    const auto &elems = moldata.property(elem_prop).asA<SireMol::AtomElements>();

                    for (int i = 0; i < nats; ++i)
                    {
                        mol_masses[i] = elems.at(idx_to_cgatomidx_data[i]).mass().to(SireUnits::g_per_mol);
                    }
                }
                catch (...)
                {
                    for (int i = 0; i < nats; ++i)
                        mol_masses[i] = 0.0;
                }
            }

            // check for light atoms or virtual sites
            for (int i = 0; i < nats; ++i)
            {
                if (mol_masses[i] < 0.05)
                {
                    mol_masses[i] = 0.0;
                }

                if (mol_masses[i] < 0.5)
                {
                    mol.virtual_sites.append(i);
                }
                else if (mol_masses[i] < 2.5)
                {
                    mol.light_atoms.append(i);
                }
            }

            try
            {
                const auto &chgs = moldata.property(chg_prop).asA<SireMol::AtomCharges>();
                auto mol_charges = mol.charges.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol_charges[i] = chgs.at(idx_to_cgatomidx_data[i]).to(SireUnits::mod_electron);
                }
            }
            catch (...)
            {
                auto mol_charges = mol.charges.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol.charges[i] = 0.0;
                }
            }

            try
            {
                const auto &ljs = moldata.property(lj_prop).asA<SireMM::AtomLJs>();
                auto mol_sigmas = mol.sigmas.data();
                auto mol_epsilons = mol.epsilons.data();

                for (int i = 0; i < nats; ++i)
                {
                    const auto &lj = ljs.at(idx_to_cgatomidx_data[i]);

                    mol_sigmas[i] = lj.sigma().to(SireUnits::nanometer);
                    mol_epsilons[i] = lj.epsilon().to(SireUnits::kJ_per_mol);
                }
            }
            catch (...)
            {
                auto mol_sigmas = mol.sigmas.data();
                auto mol_epsilons = mol.epsilons.data();

                for (int i = 0; i < nats; ++i)
                {
                    mol_sigmas[i] = 0;
                    mol_epsilons[i] = 0;
                }
            }

            return mol;
        };

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
                    openmm_mols_data[i] = to_openmm_mol(mols[i]);
                } });
        }
        else
        {
            for (int i = 0; i < nmols; ++i)
            {
                openmm_mols_data[i] = to_openmm_mol(mols[i]);
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

        // now place this data into the OpenMM system - this looks like
        // it has to be done in serial?
        int openmm_index = 0;

        for (int i = 0; i < nmols; ++i)
        {
            auto &mol = openmm_mols_data[i];

            auto mol_indexes = mol.indexes.data();
            auto mol_masses = mol.masses.constData();
            auto mol_charges = mol.charges.constData();
            auto mol_sigmas = mol.sigmas.constData();
            auto mol_epsilons = mol.epsilons.constData();

            for (int j = 0; j < mol.nats; ++j)
            {
                mol_indexes[j] = openmm_index;
                openmm_index += 1;
                system.addParticle(mol_masses[j]);
                cljff->addParticle(mol_charges[j], mol_sigmas[j], mol_epsilons[j]);
            }
        }
    }

} // end of namespace SireOpenMM
