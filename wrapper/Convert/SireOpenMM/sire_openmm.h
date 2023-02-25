#ifndef SIRE_OPENMM_H
#define SIRE_OPENMM_H

#include <OpenMM.h>

#include "SireMol/selectormol.h"

#include "SireBase/propertymap.h"

namespace SireOpenMM
{

    SireMol::SelectorMol openmm_to_sire(const OpenMM::System &mols,
                                        const SireBase::PropertyMap &map);

    OpenMM::System sire_to_openmm(const SireMol::SelectorMol &mols,
                                  const SireBase::PropertyMap &map);
}

#endif
