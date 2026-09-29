#ifndef SIRE_GEMMI_H
#define SIRE_GEMMI_H

namespace gemmi
{
    struct Structure;
}

#include "SireSystem/system.h"

#include "SireBase/propertymap.h"

#include <memory>
#include <string>

namespace SireGemmi
{

    SireSystem::System gemmi_to_sire(const gemmi::Structure &structure,
                                     const SireBase::PropertyMap &map);

    gemmi::Structure sire_to_gemmi(const SireSystem::System &system,
                                   const SireBase::PropertyMap &map);

    // Structures cross to Python as gemmi's pickle state (zpp serialized)
    SireSystem::System gemmi_state_to_sire(const std::string &state,
                                           const SireBase::PropertyMap &map);

    std::string sire_to_gemmi_state(const SireSystem::System &system,
                                    const SireBase::PropertyMap &map);

    void register_pdbx_loader();
}

#endif
