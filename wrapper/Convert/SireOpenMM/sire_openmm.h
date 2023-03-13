#ifndef SIRE_OPENMM_H
#define SIRE_OPENMM_H

#include <OpenMM.h>

#include "SireMol/selectormol.h"

#include "SireBase/propertymap.h"

namespace SireOpenMM
{
    /** This is a read-only container for an array of
     *  OpenMM::Vec3 objects (e.g. coordinates or velocities).
     *
     *  This is used internally to return coordinate and velocity
     *  data up to Python, so that it can be passed back down
     *  again to populate the openmm.Context
     */
    class CoordsAndVelocities
    {
    public:
        CoordsAndVelocities();
        CoordsAndVelocities(std::shared_ptr<std::vector<OpenMM::Vec3>> coords,
                            std::shared_ptr<std::vector<OpenMM::Vec3>> vels,
                            std::shared_ptr<std::vector<OpenMM::Vec3>> boxvecs);
        ~CoordsAndVelocities();

        bool hasCoordinates() const;
        bool hasVelocities() const;
        bool hasBoxVectors() const;

        const std::vector<OpenMM::Vec3> &coordinates() const;
        const std::vector<OpenMM::Vec3> &velocities() const;
        const std::vector<OpenMM::Vec3> &boxVectors() const;

    private:
        std::shared_ptr<std::vector<OpenMM::Vec3>> coords;
        std::shared_ptr<std::vector<OpenMM::Vec3>> vels;
        std::shared_ptr<std::vector<OpenMM::Vec3>> boxvecs;
    };

    SireMol::SelectorMol openmm_system_to_sire(const OpenMM::System &mols,
                                               const SireBase::PropertyMap &map);

    CoordsAndVelocities sire_to_openmm_system(OpenMM::System &system,
                                              const SireMol::SelectorMol &mols,
                                              const SireBase::PropertyMap &map);

    void set_openmm_coordinates_and_velocities(OpenMM::Context &context,
                                               const CoordsAndVelocities &coords_and_velocities);

    SireUnits::Dimension::MolarEnergy get_potential_energy(OpenMM::Context &context);

    SireMol::SelectorMol extract_coordinates(const OpenMM::State &state,
                                             const SireMol::SelectorMol &mols,
                                             const SireBase::PropertyMap &map);

    SireMol::SelectorMol extract_coordinates_and_velocities(const OpenMM::State &state,
                                                            const SireMol::SelectorMol &mols,
                                                            const SireBase::PropertyMap &map);
}

#endif
