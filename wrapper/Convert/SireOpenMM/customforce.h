#ifndef SIRE_OPENMM_CUSTOMFORCE_H
#define SIRE_OPENMM_CUSTOMFORCE_H

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

    protected:
        OpenMM::ForceImpl *createImpl() const;
    };
}

SIRE_END_HEADER

#endif
