#ifndef SIREOPENMM_OPENMMMINIMISE_H
#define SIREOPENMM_OPENMMMINIMISE_H

#include <OpenMM.h>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    /** This is a minimiser heavily inspired by the
     *  LocalEnergyMinimizer included in OpenMM. This is re-written
     *  for sire to;
     *
     *  1. Better integrate minimisation into the sire progress
     *     monitoring / interupting framework.
     *  2. Avoid errors caused by OpenMM switching from the desired
     *     context to the CPU context, thus triggering spurious exceptions
     *     related to exclusions / exceptions not matching
     */
    void minimise_openmm_context(OpenMM::Context &context,
                                 double tolerance = 10.0,
                                 int max_iterations = -1,
                                 int max_restarts = 10,
                                 int max_ratchets = 20,
                                 int ratchet_frequency = 500,
                                 double starting_k = 100.0,
                                 double ratchet_scale = 2.0);

}

SIRE_EXPOSE_FUNCTION(SireOpenMM::minimise_openmm_context)

SIRE_END_HEADER

#endif
