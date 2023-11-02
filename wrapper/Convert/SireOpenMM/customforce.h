#ifndef SIRE_OPENMM_CUSTOMFORCE_H
#define SIRE_OPENMM_CUSTOMFORCE_H

#include "openmm.h"
#include "openmm/Force.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    class GridForce : public OpenMM::Force
    {
    public:
        GridForce();
        ~GridForce();

        int addFixedAtom(const OpenMM::Vec3 &position,
                         float charge, float sigma, float epsilon);

        void getFixedAtom(int idx, OpenMM::Vec3 &position,
                          float &charge, float &sigma, float &epsilon) const;

        void updateFixedAtom(int idx, const OpenMM::Vec3 &position,
                             float charge, float sigma, float epsilon);

    protected:
        OpenMM::ForceImpl *createImpl() const;

    private:
        /** All of the data for the fixed atoms */
        std::vector<std::tuple<OpenMM::Vec3, float, float, float>> fixed_atoms;

        /** The coulomb potential grid */
        std::vector<double> coulomb_grid;

        /** The center of the grid */
        OpenMM::Vec3 grid_center;

        /** The half-extents of this grid (plus or minus this to
            get the grid boundaries) */
        OpenMM::Vec3 grid_half_extents;

        /** The coulomb cutoff (applies from the center of the grid) */
        double coulomb_cutoff;

        /** The grid spacing */
        double grid_spacing;

        /** The number of grid points along x, y and z */
        unsigned int dimx, dimy, dimz;
    };
}

SIRE_END_HEADER

#endif
