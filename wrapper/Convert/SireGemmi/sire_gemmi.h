#ifndef SIRE_GEMMI_H
#define SIRE_GEMMI_H

#include "gemmi/mmcif.hpp"

#include "SireSystem/system.h"

#include "SireBase/propertymap.h"

#include <memory>

namespace SireGemmi
{

    SireSystem::System gemmi_to_sire(const gemmi::Structure &structure,
                                     const SireBase::PropertyMap &map);

    gemmi::Structure sire_to_gemmi(const SireSystem::System &system,
                                   const SireBase::PropertyMap &map);

    void register_pdbx_loader();
}

#endif
