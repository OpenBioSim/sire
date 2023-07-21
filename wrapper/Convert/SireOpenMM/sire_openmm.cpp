
#include "sire_openmm.h"

#include <OpenMM.h>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireSystem/forcefieldinfo.h"

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

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

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
using SireSystem::ForceFieldInfo;

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
                                             std::shared_ptr<std::vector<OpenMM::Vec3>> v,
                                             std::shared_ptr<std::vector<OpenMM::Vec3>> b)
        : coords(c), vels(v), boxvecs(b)
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

    bool CoordsAndVelocities::hasBoxVectors() const
    {
        return boxvecs.get() != 0;
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

    const std::vector<OpenMM::Vec3> &CoordsAndVelocities::boxVectors() const
    {
        if (boxvecs.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no box vectors available!"),
                                                CODELOC);

        return *boxvecs;
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
        ForceFieldInfo ffinfo(mols, map);

        // extract the data from all of the molecules
        const int nmols = mols.count();

        if (nmols == 0)
        {
            // nothing to do
            return CoordsAndVelocities();
        }

        // whether or not to ignore perturbations
        bool ignore_perturbations = false;
        bool any_perturbable = false;

        if (map.specified("ignore_perturbations"))
        {
            ignore_perturbations = map["ignore_perturbations"].value().asABoolean();
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

        // check to see if there are any perturbable molecules
        if (not ignore_perturbations)
        {
            for (int i = 0; i < nmols; ++i)
            {
                if (openmm_mols_data[i].isPerturbable())
                {
                    any_perturbable = true;
                    break;
                }
            }
        }

        // create all the OpenMM forcefields - ownership is taken by 'system'
        // (memory leak until we give this to 'system')

        // CLJ energy between all non-perturbable atoms
        OpenMM::NonbondedForce *cljff = new OpenMM::NonbondedForce();

        bool use_dispersion_correction = false;

        if (map.specified("use_dispersion_correction"))
        {
            use_dispersion_correction = map["use_dispersion_correction"].value().asABoolean();
        }

        cljff->setUseDispersionCorrection(use_dispersion_correction);

        // create the periodic box vectors
        std::shared_ptr<std::vector<OpenMM::Vec3>> boxvecs;

        if (ffinfo.space().isPeriodic())
        {
            boxvecs.reset(new std::vector<OpenMM::Vec3>(3));
            auto boxvecs_data = boxvecs->data();

            if (ffinfo.space().isA<SireVol::PeriodicBox>())
            {
                const auto &space = ffinfo.space().asA<SireVol::PeriodicBox>();

                const double x = space.dimensions()[0] * OpenMM::NmPerAngstrom;
                const double y = space.dimensions()[1] * OpenMM::NmPerAngstrom;
                const double z = space.dimensions()[2] * OpenMM::NmPerAngstrom;

                boxvecs_data[0] = OpenMM::Vec3(x, 0, 0);
                boxvecs_data[1] = OpenMM::Vec3(0, y, 0);
                boxvecs_data[2] = OpenMM::Vec3(0, 0, z);
            }
            else if (ffinfo.space().isA<SireVol::TriclinicBox>())
            {
                const auto &space = ffinfo.space().asA<SireVol::TriclinicBox>();

                // Get the three triclinic box vectors.
                const auto v0 = space.vector0();
                const auto v1 = space.vector1();
                const auto v2 = space.vector2();

                // Get cell matrix components in nm.
                const double xx = v0.x() * OpenMM::NmPerAngstrom;
                const double xy = v0.y() * OpenMM::NmPerAngstrom;
                const double xz = v0.z() * OpenMM::NmPerAngstrom;
                const double yx = v1.x() * OpenMM::NmPerAngstrom;
                const double yy = v1.y() * OpenMM::NmPerAngstrom;
                const double yz = v1.z() * OpenMM::NmPerAngstrom;
                const double zx = v2.x() * OpenMM::NmPerAngstrom;
                const double zy = v2.y() * OpenMM::NmPerAngstrom;
                const double zz = v2.z() * OpenMM::NmPerAngstrom;

                boxvecs_data[0] = OpenMM::Vec3(xx, xy, xz);
                boxvecs_data[1] = OpenMM::Vec3(yx, yy, yz);
                boxvecs_data[2] = OpenMM::Vec3(zx, zy, zz);
            }

            system.setDefaultPeriodicBoxVectors(boxvecs_data[0],
                                                boxvecs_data[1],
                                                boxvecs_data[2]);
        }

        if (ffinfo.hasCutoff())
        {
            const auto typ = ffinfo.cutoffType();

            // need to set cutoff for perturbable forcefields

            if (typ == "PME" or typ == "EWALD")
            {
                if (not ffinfo.space().isPeriodic())
                {
                    throw SireError::incompatible_error(QObject::tr(
                                                            "You cannot use Ewald or PME with the non-periodic space %1.")
                                                            .arg(ffinfo.space().toString()),
                                                        CODELOC);
                }

                auto nbmethod = OpenMM::NonbondedForce::PME;

                if (typ != "PME")
                    nbmethod = OpenMM::NonbondedForce::Ewald;

                cljff->setNonbondedMethod(nbmethod);

                double tolerance = ffinfo.getParameter("tolerance").value();

                if (tolerance <= 0)
                    tolerance = 0.001;

                cljff->setEwaldErrorTolerance(tolerance);
            }
            else if (typ == "REACTION_FIELD")
            {
                auto nbmethod = OpenMM::NonbondedForce::CutoffPeriodic;

                if (not ffinfo.space().isPeriodic())
                {
                    nbmethod = OpenMM::NonbondedForce::CutoffNonPeriodic;
                }

                auto dielectric = ffinfo.getParameter("dielectric").value();

                if (dielectric <= 0)
                    dielectric = 78.3;

                cljff->setNonbondedMethod(nbmethod);
                cljff->setReactionFieldDielectric(dielectric);
            }
            else if (typ == "CUTOFF")
            {
                // use reaction field for non-periodic spaces, and PME for periodic
                if (ffinfo.space().isPeriodic())
                {
                    const auto nbmethod = OpenMM::NonbondedForce::PME;
                    cljff->setNonbondedMethod(nbmethod);

                    double tolerance = ffinfo.getParameter("tolerance").value();

                    if (tolerance <= 0)
                        tolerance = 0.001;

                    cljff->setEwaldErrorTolerance(tolerance);
                }
                else
                {
                    const auto nbmethod = OpenMM::NonbondedForce::CutoffNonPeriodic;
                    cljff->setNonbondedMethod(nbmethod);

                    double dielectric = ffinfo.getParameter("dielectric").value();

                    if (dielectric <= 0)
                        dielectric = 78.3;

                    cljff->setReactionFieldDielectric(dielectric);
                }
            }

            const auto cutoff = ffinfo.cutoff().to(SireUnits::nanometers);
            cljff->setCutoffDistance(cutoff);
        }

        system.addForce(cljff);

        OpenMM::HarmonicBondForce *bondff = new OpenMM::HarmonicBondForce();
        system.addForce(bondff);

        OpenMM::HarmonicAngleForce *angff = new OpenMM::HarmonicAngleForce();
        system.addForce(angff);

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

            // is this a perturbable molecule (and we haven't disabled perturbations)?
            // (not used now, but we will need to use this when we add softening)
            // bool is_perturbable_mol = any_perturbable and mol.isPerturbable();

            // first the atom parameters
            auto masses_data = mol.masses.constData();
            auto cljs_data = mol.cljs.constData();

            for (int j = 0; j < mol.molinfo.nAtoms(); ++j)
            {
                system.addParticle(masses_data[j]);
                const int atom_index = start_index + j;

                const auto &clj = cljs_data[j];

                cljff->addParticle(std::get<0>(clj), std::get<1>(clj), std::get<2>(clj));
            }

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

            // now constraints
            for (const auto &constraint : mol.constraints)
            {
                system.addConstraint(std::get<0>(constraint) + start_index,
                                     std::get<1>(constraint) + start_index,
                                     std::get<2>(constraint));
            }

            start_index += mol.masses.count();
        }

        const int natoms = start_index;

        // add exclusions based on the bonding of the molecules
        cljff->createExceptionsFromBonds(bond_pairs, coul_14_scl, lj_14_scl);

        for (const auto &p : custom_pairs)
        {
            cljff->addException(std::get<0>(p), std::get<1>(p),
                                std::get<2>(p), std::get<3>(p),
                                std::get<4>(p), false);
        }

        // will have to add exceptions for perturbable forces

        // see if we want to remove COM motion
        const auto com_remove_prop = map["com_reset_frequency"];

        if (com_remove_prop.hasValue())
        {
            const int freq = com_remove_prop.value().asAnInteger();

            if (freq > 0)
            {
                OpenMM::CMMotionRemover *com_remover = new OpenMM::CMMotionRemover(freq);
                system.addForce(com_remover);
            }
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

        return CoordsAndVelocities(coords, vels, boxvecs);
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

        if (coords_and_velocities.hasBoxVectors())
        {
            const auto boxvecs = coords_and_velocities.boxVectors();

            context.setPeriodicBoxVectors(boxvecs[0], boxvecs[1], boxvecs[2]);
        }
    }

    inline void _populate_coords(QVector<SireMaths::Vector> &coords,
                                 const OpenMM::Vec3 *omm_coords,
                                 int natoms)
    {
        coords.resize(natoms);
        auto coords_data = coords.data();

        const double nm_to_internal = (1 * SireUnits::nanometer).to(SireUnits::angstrom);

        for (int i = 0; i < natoms; ++i)
        {
            const auto &omm = omm_coords[i];
            coords_data[i] = SireMaths::Vector(omm[0] * nm_to_internal,
                                               omm[1] * nm_to_internal,
                                               omm[2] * nm_to_internal);
        }
    }

    inline void _populate_vels(QVector<SireMol::Velocity3D> &vels,
                               const OpenMM::Vec3 *omm_vels,
                               int natoms)
    {
        vels.resize(natoms);
        auto vels_data = vels.data();

        for (int i = 0; i < natoms; ++i)
        {
            const auto &omm = omm_vels[i];
            vels_data[i] = SireMol::Velocity3D(omm[0] * SireUnits::nanometers_per_ps,
                                               omm[1] * SireUnits::nanometers_per_ps,
                                               omm[2] * SireUnits::nanometers_per_ps);
        }
    }

    SireVol::SpacePtr extract_space(const OpenMM::State &state)
    {
        OpenMM::Vec3 a, b, c;

        try
        {
            state.getPeriodicBoxVectors(a, b, c);
        }
        catch (...)
        {
            return SireVol::SpacePtr(SireVol::Cartesian());
        }

        const double nm_to_internal = (1 * SireUnits::nanometer).to(SireUnits::angstrom);

        SireMaths::Vector x(a[0] * nm_to_internal,
                            a[1] * nm_to_internal,
                            a[2] * nm_to_internal);

        SireMaths::Vector y(b[0] * nm_to_internal,
                            b[1] * nm_to_internal,
                            b[2] * nm_to_internal);

        SireMaths::Vector z(c[0] * nm_to_internal,
                            c[1] * nm_to_internal,
                            c[2] * nm_to_internal);

        SireVol::TriclinicBox triclinic;

        try
        {
            triclinic = SireVol::TriclinicBox(x, y, z);
        }
        catch (...)
        {
            // this is not a valid space - could be an infinite space
            return SireVol::SpacePtr(SireVol::Cartesian());
        }

        if (triclinic.alpha() == 90 and triclinic.beta() == 90 and triclinic.gamma() == 90)
        {
            // this is a PeriodicBox?
            SireVol::PeriodicBox pbox(x.max(y).max(z) - x.min(y).min(z));

            if (std::abs(pbox.volume().value() - triclinic.volume().value()) < 0.001)
            {
                // yes - periodic box
                return SireVol::SpacePtr(pbox);
            }
        }

        return SireVol::SpacePtr(triclinic);
    }

    SelectorMol extract_coordinates_and_velocities(const OpenMM::State &state,
                                                   const SelectorMol &mols,
                                                   const PropertyMap &map)
    {
        const auto positions = state.getPositions();
        const auto velocities = state.getVelocities();

        const int natoms = positions.size();
        const auto positions_data = positions.data();
        const auto velocities_data = velocities.data();

        if (mols.nAtoms() != natoms)
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "Different number of atoms from OpenMM and sire. "
                                                    "%1 versus %2. Cannot extract the coordinates.")
                                                    .arg(natoms)
                                                    .arg(mols.nAtoms()),
                                                CODELOC);
        }

        const auto coords_prop = map["coordinates"];
        const auto vels_prop = map["velocity"];

        const int nmols = mols.count();

        QVector<Molecule> ret(nmols);
        auto ret_data = ret.data();

        QVector<int> offsets(nmols);

        int offset = 0;

        for (int i = 0; i < nmols; ++i)
        {
            offsets[i] = offset;
            offset += mols[i].nAtoms();
        }

        const auto offsets_data = offsets.constData();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                QVector<SireMaths::Vector> converted_coords;
                QVector<SireMol::Velocity3D> converted_vels;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    auto mol = mols[i].edit();
                    const int mol_natoms = mol.nAtoms();

                    _populate_coords(converted_coords, positions_data+offsets_data[i], mol_natoms);
                    _populate_vels(converted_vels, velocities_data+offsets_data[i], mol_natoms);

                    SireMol::AtomCoords mol_coords;
                    SireMol::AtomVelocities mol_vels;

                    if (mol.data().hasProperty(coords_prop))
                    {
                        mol_coords = mol.data().property(coords_prop).asA<SireMol::AtomCoords>();
                    }
                    else
                    {
                        mol_coords = SireMol::AtomCoords(mol.data().info());
                    }

                    if (mol.data().hasProperty(vels_prop))
                    {
                        mol_vels = mol.data().property(vels_prop).asA<SireMol::AtomVelocities>();
                    }
                    else
                    {
                        mol_vels = SireMol::AtomVelocities(mol.data().info());
                    }

                    mol_coords.copyFrom(converted_coords);
                    mol_vels.copyFrom(converted_vels);

                    mol.setProperty(coords_prop.source(), mol_coords);
                    mol.setProperty(vels_prop.source(), mol_vels);

                    ret_data[i] = mol.commit();
                } });
        }
        else
        {
            QVector<SireMaths::Vector> converted_coords;
            QVector<SireMol::Velocity3D> converted_vels;

            for (int i = 0; i < nmols; ++i)
            {
                auto mol = mols[i].edit();
                const int mol_natoms = mol.nAtoms();

                _populate_coords(converted_coords, positions_data + offsets_data[i], mol_natoms);
                _populate_vels(converted_vels, velocities_data + offsets_data[i], mol_natoms);

                SireMol::AtomCoords mol_coords;
                SireMol::AtomVelocities mol_vels;

                if (mol.data().hasProperty(coords_prop))
                {
                    mol_coords = mol.data().property(coords_prop).asA<SireMol::AtomCoords>();
                }
                else
                {
                    mol_coords = SireMol::AtomCoords(mol.data().info());
                }

                if (mol.data().hasProperty(vels_prop))
                {
                    mol_vels = mol.data().property(vels_prop).asA<SireMol::AtomVelocities>();
                }
                else
                {
                    mol_vels = SireMol::AtomVelocities(mol.data().info());
                }

                mol_coords.copyFrom(converted_coords);
                mol_vels.copyFrom(converted_vels);

                mol.setProperty(coords_prop.source(), mol_coords);
                mol.setProperty(vels_prop.source(), mol_vels);

                ret_data[i] = mol.commit();
            }
        }

        return SelectorMol(ret);
    }

    SelectorMol extract_coordinates(const OpenMM::State &state,
                                    const SelectorMol &mols,
                                    const PropertyMap &map)
    {
        const auto positions = state.getPositions();

        const int natoms = positions.size();
        const auto positions_data = positions.data();

        if (mols.nAtoms() != natoms)
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "Different number of atoms from OpenMM and sire. "
                                                    "%1 versus %2. Cannot extract the coordinates.")
                                                    .arg(natoms)
                                                    .arg(mols.nAtoms()),
                                                CODELOC);
        }

        const auto coords_prop = map["coordinates"];

        const int nmols = mols.count();

        QVector<Molecule> ret(nmols);
        auto ret_data = ret.data();

        QVector<int> offsets(nmols);

        int offset = 0;

        for (int i = 0; i < nmols; ++i)
        {
            offsets[i] = offset;
            offset += mols[i].nAtoms();
        }

        const auto offsets_data = offsets.constData();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                QVector<SireMaths::Vector> converted_coords;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    auto mol = mols[i].edit();
                    const int mol_natoms = mol.nAtoms();

                    _populate_coords(converted_coords, positions_data+offsets_data[i], mol_natoms);

                    SireMol::AtomCoords mol_coords;

                    if (mol.data().hasProperty(coords_prop))
                    {
                        mol_coords = mol.data().property(coords_prop).asA<SireMol::AtomCoords>();
                    }
                    else
                    {
                        mol_coords = SireMol::AtomCoords(mol.data().info());
                    }

                    mol_coords.copyFrom(converted_coords);

                    mol.setProperty(coords_prop.source(), mol_coords);

                    ret_data[i] = mol.commit();
                } });
        }
        else
        {
            QVector<SireMaths::Vector> converted_coords;

            for (int i = 0; i < nmols; ++i)
            {
                auto mol = mols[i].edit();
                const int mol_natoms = mol.nAtoms();

                _populate_coords(converted_coords, positions_data + offsets_data[i], mol_natoms);

                SireMol::AtomCoords mol_coords;

                if (mol.data().hasProperty(coords_prop))
                {
                    mol_coords = mol.data().property(coords_prop).asA<SireMol::AtomCoords>();
                }
                else
                {
                    mol_coords = SireMol::AtomCoords(mol.data().info());
                }

                mol_coords.copyFrom(converted_coords);

                mol.setProperty(coords_prop.source(), mol_coords);

                ret_data[i] = mol.commit();
            }
        }

        return SelectorMol(ret);
    }

} // end of namespace SireOpenMM
