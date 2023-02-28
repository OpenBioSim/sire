
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

#include "openmmmolecule.h"

#include <QDebug>

using SireBase::PropertyMap;
using SireMol::Molecule;
using SireMol::SelectorMol;

namespace SireOpenMM
{
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

    CoordsAndVelocities::CoordsAndVelocities()
    {
    }

    CoordsAndVelocities::CoordsAndVelocities(std::shared_ptr<std::vector<OpenMM::Vec3>> c,
                                             std::shared_ptr<std::vector<OpenMM::Vec3>> v)
        : coords(c), vels(v)
    {
    }

    CoordsAndVelocities::~CoordsAndVelocities()
    {
    }

    bool CoordsAndVelocities::hasCoordinates() const
    {
        return coords.get() != 0;
    }

    bool CoordsAndVelocities::hasVelocities() const
    {
        return vels.get() != 0;
    }

    const std::vector<OpenMM::Vec3> &CoordsAndVelocities::coordinates() const
    {
        if (coords.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no coordinates available!"),
                                                CODELOC);

        return *coords;
    }

    const std::vector<OpenMM::Vec3> &CoordsAndVelocities::velocities() const
    {
        if (vels.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no velocities available!"),
                                                CODELOC);

        return *vels;
    }

    SelectorMol openmm_system_to_sire(const OpenMM::System &mols,
                                      const PropertyMap &map)
    {
        throw SireError::incomplete_code(QObject::tr(
                                             "Still need to write openmm_to_sire"),
                                         CODELOC);

        return SelectorMol();
    }

    CoordsAndVelocities sire_to_openmm_system(OpenMM::System &system,
                                              const SelectorMol &mols,
                                              const PropertyMap &map)
    {
        // we can assume that an empty system has been passed to us

        // we will get parameters from the map (e.g. cutoffs etc)

        // extract the data from all of the molecules
        const int nmols = mols.count();

        if (nmols == 0)
        {
            // nothing to do
            return CoordsAndVelocities();
        }

        // Extract all of the data needed by OpenMM from the Sire
        // molecules into some temporary OpenMMMolecule objects
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

        // create all the OpenMM forcefields - ownership is taken by 'system'
        OpenMM::NonbondedForce *cljff = new OpenMM::NonbondedForce();
        // system.addForce(cljff);

        OpenMM::HarmonicBondForce *bondff = new OpenMM::HarmonicBondForce();
        // system.addForce(bondff);

        OpenMM::HarmonicAngleForce *angff = new OpenMM::HarmonicAngleForce();
        // system.addForce(angff);

        OpenMM::PeriodicTorsionForce *dihff = new OpenMM::PeriodicTorsionForce();
        system.addForce(dihff);

        // Now copy data from the temporary OpenMMMolecule objects
        // into these forcefields
        // (will deal with restraints, light atoms and virtual sites later)

        int start_index = 0;

        std::vector<std::pair<int, int>> bond_pairs;
        std::vector<std::tuple<int, int, double, double, double>> custom_pairs;

        // get the 1-4 scaling factors from the first molecule
        const double coul_14_scl = openmm_mols_data[0].ffinfo.electrostatic14ScaleFactor();
        const double lj_14_scl = openmm_mols_data[0].ffinfo.vdw14ScaleFactor();

        QVector<int> start_indexes(nmols);

        for (int i = 0; i < nmols; ++i)
        {
            start_indexes[i] = start_index;
            const auto &mol = openmm_mols_data[i];

            if (std::abs(mol.ffinfo.electrostatic14ScaleFactor() - coul_14_scl) > 0.001 or
                std::abs(mol.ffinfo.vdw14ScaleFactor() - lj_14_scl) > 0.001)
            {
                throw SireError::incompatible_error(QObject::tr(
                                                        "We cannot create the OpenMM system because the forcefields "
                                                        "of %1 and %2 (%3 and %4) are not compatible. They have "
                                                        "different 1-4 non-bonded scale factors.")
                                                        .arg(mols[0].toString())
                                                        .arg(mols[i].toString())
                                                        .arg(openmm_mols_data[0].ffinfo.toString())
                                                        .arg(mol.ffinfo.toString()),
                                                    CODELOC);
            }

            // first the atom parameters
            auto masses_data = mol.masses.constData();
            auto cljs_data = mol.cljs.constData();

            for (int j = 0; j < mol.molinfo.nAtoms(); ++j)
            {
                system.addParticle(masses_data[j]);
                const auto &clj = cljs_data[j];
                cljff->addParticle(std::get<0>(clj), std::get<1>(clj), std::get<2>(clj));
            }

            // now the connectivity
            for (const auto &bond : mol.bond_pairs)
            {
                bond_pairs.push_back(std::make_pair(std::get<0>(bond) + start_index,
                                                    std::get<1>(bond) + start_index));
            }

            // now any custom pairs
            for (const auto &pair : mol.custom_pairs)
            {
                custom_pairs.push_back(std::make_tuple(std::get<0>(pair) + start_index,
                                                       std::get<1>(pair) + start_index,
                                                       std::get<2>(pair),
                                                       std::get<3>(pair),
                                                       std::get<4>(pair)));
            }

            // now bond parameters
            for (const auto &bond : mol.bond_params)
            {
                bondff->addBond(std::get<0>(bond) + start_index,
                                std::get<1>(bond) + start_index,
                                std::get<2>(bond), std::get<3>(bond));
            }

            // now the angles
            for (const auto &ang : mol.ang_params)
            {
                angff->addAngle(std::get<0>(ang) + start_index,
                                std::get<1>(ang) + start_index,
                                std::get<2>(ang) + start_index,
                                std::get<3>(ang), std::get<4>(ang));
            }

            // now the dihedrals and impropers
            for (const auto &dih : mol.dih_params)
            {
                dihff->addTorsion(std::get<0>(dih) + start_index,
                                  std::get<1>(dih) + start_index,
                                  std::get<2>(dih) + start_index,
                                  std::get<3>(dih) + start_index,
                                  std::get<4>(dih), std::get<5>(dih), std::get<6>(dih));
            }

            start_index += mol.masses.count();
        }

        const int natoms = start_index;

        // add exclusions based on the bonding of this molecule
        cljff->createExceptionsFromBonds(bond_pairs, coul_14_scl, lj_14_scl);

        for (const auto &p : custom_pairs)
        {
            cljff->addException(std::get<0>(p), std::get<1>(p),
                                std::get<2>(p), std::get<3>(p),
                                std::get<4>(p), false);
        }

        // now get the coordinates and velocities
        std::shared_ptr<std::vector<OpenMM::Vec3>> coords, vels;

        coords.reset(new std::vector<OpenMM::Vec3>(natoms));
        vels.reset(new std::vector<OpenMM::Vec3>(natoms));

        auto coords_data = coords->data();
        auto vels_data = vels->data();

        const int *start_indexes_data = start_indexes.constData();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const int start_index = start_indexes_data[i];
                    const auto &mol = openmm_mols_data[i];
                    mol.copyInCoordsAndVelocities(coords_data + start_index,
                                                vels_data + start_index);
                } });
        }
        else
        {
            for (int i = 0; i < nmols; ++i)
            {
                const int start_index = start_indexes_data[i];
                const auto &mol = openmm_mols_data[i];
                mol.copyInCoordsAndVelocities(coords_data + start_index,
                                              vels_data + start_index);
            }
        }

        return CoordsAndVelocities(coords, vels);
    }

    void set_openmm_coordinates_and_velocities(OpenMM::Context &context,
                                               const CoordsAndVelocities &coords_and_velocities)
    {
        if (coords_and_velocities.hasCoordinates())
        {
            context.setPositions(coords_and_velocities.coordinates());
        }

        if (coords_and_velocities.hasVelocities())
        {
            context.setVelocities(coords_and_velocities.velocities());
        }
    }

} // end of namespace SireOpenMM
