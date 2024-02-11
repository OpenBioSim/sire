
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

#include "SireCAS/lambdaschedule.h"

#include "SireMaths/vector.h"

#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "tostring.h"

#include "openmmmolecule.h"

#include <QDebug>

using SireBase::PropertyMap;
using SireCAS::LambdaSchedule;
using SireMol::Molecule;
using SireMol::MolNum;
using SireMol::SelectorMol;
using SireSystem::ForceFieldInfo;

namespace SireOpenMM
{
    ////
    //// Implementation of OpenMMMetaData
    ////

    OpenMMMetaData::OpenMMMetaData()
    {
    }

    OpenMMMetaData::OpenMMMetaData(const SireMol::SelectorM<SireMol::Atom> &i,
                                   std::shared_ptr<std::vector<OpenMM::Vec3>> c,
                                   std::shared_ptr<std::vector<OpenMM::Vec3>> v,
                                   std::shared_ptr<std::vector<OpenMM::Vec3>> b,
                                   const LambdaLever &l)
        : atom_index(i), coords(c), vels(v), boxvecs(b), lambda_lever(l)
    {
    }

    OpenMMMetaData::~OpenMMMetaData()
    {
    }

    const char *OpenMMMetaData::typeName()
    {
        return "SireOpenMM::OpenMMMetaData";
    }

    const char *OpenMMMetaData::what() const
    {
        return OpenMMMetaData::typeName();
    }

    QString OpenMMMetaData::toString() const
    {
        return QString("OpenMMMetaData()");
    }

    SireMol::SelectorM<SireMol::Atom> OpenMMMetaData::index() const
    {
        return atom_index;
    }

    LambdaLever OpenMMMetaData::lambdaLever() const
    {
        return lambda_lever;
    }

    bool OpenMMMetaData::hasCoordinates() const
    {
        return coords.get() != 0;
    }

    bool OpenMMMetaData::hasVelocities() const
    {
        return vels.get() != 0;
    }

    bool OpenMMMetaData::hasBoxVectors() const
    {
        return boxvecs.get() != 0;
    }

    const std::vector<OpenMM::Vec3> &OpenMMMetaData::coordinates() const
    {
        if (coords.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no coordinates available!"),
                                                CODELOC);

        return *coords;
    }

    const std::vector<OpenMM::Vec3> &OpenMMMetaData::velocities() const
    {
        if (vels.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no velocities available!"),
                                                CODELOC);

        return *vels;
    }

    const std::vector<OpenMM::Vec3> &OpenMMMetaData::boxVectors() const
    {
        if (boxvecs.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no box vectors available!"),
                                                CODELOC);

        return *boxvecs;
    }

    ////
    //// Implementation of standalone functions
    ////

    SelectorMol openmm_system_to_sire(const OpenMM::System &mols,
                                      const PropertyMap &map)
    {
        throw SireError::incomplete_code(QObject::tr(
                                             "Still need to write openmm_to_sire"),
                                         CODELOC);

        return SelectorMol();
    }

    void set_openmm_coordinates_and_velocities(OpenMM::Context &context,
                                               const OpenMMMetaData &metadata)
    {
        if (metadata.hasCoordinates())
        {
            context.setPositions(metadata.coordinates());
        }

        if (metadata.hasVelocities())
        {
            context.setVelocities(metadata.velocities());
        }

        if (metadata.hasBoxVectors())
        {
            const auto boxvecs = metadata.boxVectors();

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

    void set_context_platform_property(OpenMM::Context &context,
                                       const QString &key,
                                       const QString &value)
    {
        OpenMM::Platform &platform = context.getPlatform();

        platform.setPropertyValue(context,
                                  key.toStdString(),
                                  value.toStdString());

        QString new_value = QString::fromStdString(platform.getPropertyValue(context, key.toStdString()));

        if (new_value != value)
            throw SireError::incompatible_error(QObject::tr(
                                                    "Unable to change the value of property %1 to `%2` in the "
                                                    "platform %3. The property value is still '%4'.")
                                                    .arg(key)
                                                    .arg(value)
                                                    .arg(QString::fromStdString(platform.getName()))
                                                    .arg(new_value),
                                                CODELOC);
    }

    SelectorMol extract_coordinates(const OpenMM::State &state,
                                    const SireMol::SelectorMol &mols,
                                    const QHash<SireMol::MolNum, SireBase::PropertyMap> &perturbable_maps,
                                    const SireBase::PropertyMap &map)
    {
        const auto positions = state.getPositions();

        const int natoms = positions.size();
        const auto positions_data = positions.data();

        if (mols.nAtoms() > natoms)
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

        if (SireBase::should_run_in_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                QVector<SireMaths::Vector> converted_coords;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    auto mol = mols[i].edit();
                    const int mol_natoms = mol.nAtoms();

                    _populate_coords(converted_coords, positions_data+offsets_data[i], mol_natoms);

                    auto my_coords_prop = coords_prop.source();

                    if (perturbable_maps.contains(mol.number()))
                    {
                        my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                    }

                    if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                        converted_coords, false))
                    {
                        SireMol::AtomCoords c(mol.data().info());
                        c.copyFrom(converted_coords);
                        mol.setProperty(my_coords_prop, c);
                    }

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

                auto my_coords_prop = coords_prop.source();

                if (perturbable_maps.contains(mol.number()))
                {
                    my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                }

                if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                    converted_coords, false))
                {
                    SireMol::AtomCoords c(mol.data().info());
                    c.copyFrom(converted_coords);
                    mol.setProperty(my_coords_prop, c);
                }

                ret_data[i] = mol.commit();
            }
        }

        return SelectorMol(ret);
    }

    SelectorMol extract_coordinates_and_velocities(const OpenMM::State &state,
                                                   const SireMol::SelectorMol &mols,
                                                   const QHash<SireMol::MolNum, SireBase::PropertyMap> &perturbable_maps,
                                                   const SireBase::PropertyMap &map)
    {
        const auto positions = state.getPositions();
        const auto velocities = state.getVelocities();

        const int natoms = positions.size();
        const auto positions_data = positions.data();
        const auto velocities_data = velocities.data();

        if (mols.nAtoms() > natoms)
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

        if (SireBase::should_run_in_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                QVector<SireMaths::Vector> converted_coords;
                QVector<SireMol::Velocity3D> converted_vels;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    auto mol = mols[i].edit();
                    const int mol_natoms = mol.nAtoms();

                    auto my_coords_prop = coords_prop.source();
                    auto my_vels_prop = vels_prop.source();

                    if (perturbable_maps.contains(mol.number()))
                    {
                        my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                        my_vels_prop = perturbable_maps[mol.number()]["velocity"].source();
                    }

                    _populate_coords(converted_coords, positions_data+offsets_data[i], mol_natoms);
                    _populate_vels(converted_vels, velocities_data+offsets_data[i], mol_natoms);

                    if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                        converted_coords, false))
                    {
                        SireMol::AtomCoords c(mol.data().info());
                        c.copyFrom(converted_coords);
                        mol.setProperty(my_coords_prop, c);
                    }

                    if (not mol.updatePropertyFrom<SireMol::AtomVelocities>(my_vels_prop,
                                                                            converted_vels, false))
                    {
                        SireMol::AtomVelocities v(mol.data().info());
                        v.copyFrom(converted_vels);
                        mol.setProperty(my_vels_prop, v);
                    }

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

                auto my_coords_prop = coords_prop.source();
                auto my_vels_prop = vels_prop.source();

                if (perturbable_maps.contains(mol.number()))
                {
                    my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                    my_vels_prop = perturbable_maps[mol.number()]["velocity"].source();
                }

                _populate_coords(converted_coords, positions_data + offsets_data[i], mol_natoms);
                _populate_vels(converted_vels, velocities_data + offsets_data[i], mol_natoms);

                if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                    converted_coords, false))
                {
                    SireMol::AtomCoords c(mol.data().info());
                    c.copyFrom(converted_coords);
                    mol.setProperty(my_coords_prop, c);
                }

                if (not mol.updatePropertyFrom<SireMol::AtomVelocities>(my_vels_prop,
                                                                        converted_vels, false))
                {
                    SireMol::AtomVelocities v(mol.data().info());
                    v.copyFrom(converted_vels);
                    mol.setProperty(my_vels_prop, v);
                }

                ret_data[i] = mol.commit();
            }
        }

        return SelectorMol(ret);
    }

    SireUnits::Dimension::MolarEnergy get_potential_energy(OpenMM::Context &context)
    {
        return context.getState(OpenMM::State::Energy).getPotentialEnergy() * SireUnits::kJ_per_mol;
    }

} // end of namespace SireOpenMM
