#ifndef SIRE_OPENMM_H
#define SIRE_OPENMM_H

#include <OpenMM.h>

#include "SireMol/selectormol.h"
#include "SireVol/space.h"

#include "SireMol/core.h"
#include "SireMol/selectorm.hpp"
#include "SireMol/atom.h"
#include "SireMol/molnum.h"

#include "SireBase/propertymap.h"

#include "openmmmolecule.h"
#include "lambdalever.h"

namespace SireOpenMM
{
    /** This is a read-only container for the extra information
     *  we need to create an OpenMM system and make it editable
     *
     *  This is used for a number of things:
     *
     *  1. For an array of
     *     OpenMM::Vec3 objects (e.g. coordinates or velocities).
     *     This is used internally to return coordinate and velocity
     *     data up to Python, so that it can be passed back down
     *     again to populate the openmm.Context
     *
     *  2. Index - this is a SelectorM<Atom> in the order that the
     *     atoms appear in the system. This lets us easily locate
     *     atoms by searching
     *
     *  3. LambdaLever - the simplest LambdaLever which can be
     *     used to update the parameters of the system containing
     *     perturbable atoms between the two end states
     */
    class OpenMMMetaData
    {
    public:
        OpenMMMetaData();
        OpenMMMetaData(const SireMol::SelectorM<SireMol::Atom> &atoms,
                       std::shared_ptr<std::vector<OpenMM::Vec3>> coords,
                       std::shared_ptr<std::vector<OpenMM::Vec3>> vels,
                       std::shared_ptr<std::vector<OpenMM::Vec3>> boxvecs,
                       const LambdaLever &lever);
        ~OpenMMMetaData();

        static const char *typeName();
        const char *what() const;

        QString toString() const;

        bool hasCoordinates() const;
        bool hasVelocities() const;
        bool hasBoxVectors() const;

        const std::vector<OpenMM::Vec3> &coordinates() const;
        const std::vector<OpenMM::Vec3> &velocities() const;
        const std::vector<OpenMM::Vec3> &boxVectors() const;

        SireMol::SelectorM<SireMol::Atom> index() const;

        LambdaLever lambdaLever() const;

    private:
        std::shared_ptr<std::vector<OpenMM::Vec3>> coords;
        std::shared_ptr<std::vector<OpenMM::Vec3>> vels;
        std::shared_ptr<std::vector<OpenMM::Vec3>> boxvecs;

        SireMol::SelectorM<SireMol::Atom> atom_index;
        LambdaLever lambda_lever;
    };

    SireMol::SelectorMol openmm_system_to_sire(const OpenMM::System &system,
                                               const SireBase::PropertyMap &map);

    OpenMMMetaData sire_to_openmm_system(OpenMM::System &system,
                                         const SireMol::SelectorMol &mols,
                                         const SireBase::PropertyMap &map);

    void set_openmm_coordinates_and_velocities(OpenMM::Context &context,
                                               const OpenMMMetaData &coords_and_velocities);

    SireUnits::Dimension::MolarEnergy get_potential_energy(OpenMM::Context &context);

    SireMol::SelectorMol extract_coordinates(const OpenMM::State &state,
                                             const SireMol::SelectorMol &mols,
                                             const QHash<SireMol::MolNum, SireBase::PropertyMap> &perturbable_maps,
                                             const SireBase::PropertyMap &map);

    SireMol::SelectorMol extract_coordinates_and_velocities(const OpenMM::State &state,
                                                            const SireMol::SelectorMol &mols,
                                                            const QHash<SireMol::MolNum, SireBase::PropertyMap> &perturbable_maps,
                                                            const SireBase::PropertyMap &map);

    SireVol::SpacePtr extract_space(const OpenMM::State &state);

    void set_context_platform_property(OpenMM::Context &context,
                                       const QString &key,
                                       const QString &value);
}

SIRE_EXPOSE_CLASS(SireOpenMM::OpenMMMetaData)

SIRE_EXPOSE_FUNCTION(SireOpenMM::openmm_system_to_sire)
SIRE_EXPOSE_FUNCTION(SireOpenMM::sire_to_openmm_system)
SIRE_EXPOSE_FUNCTION(SireOpenMM::set_openmm_coordinates_and_velocities)
SIRE_EXPOSE_FUNCTION(SireOpenMM::get_potential_energy)
SIRE_EXPOSE_FUNCTION(SireOpenMM::extract_coordinates)
SIRE_EXPOSE_FUNCTION(SireOpenMM::extract_coordinates_and_velocities)
SIRE_EXPOSE_FUNCTION(SireOpenMM::extract_space)
SIRE_EXPOSE_FUNCTION(SireOpenMM::set_context_platform_property)

#endif
