
#include "sire_openmm.h"

#include <OpenMM.h>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

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

#include "SireMM/selectorbond.h"

#include "SireMaths/vector.h"

#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"

#include "SireUnits/units.h"

#include "tostring.h"

#include <QDebug>

using SireBase::PropertyMap;
using SireMol::Molecule;
using SireMol::SelectorMol;

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

    SelectorMol openmm_to_sire(const OpenMM::System &mols,
                               const PropertyMap &map)
    {
        return SelectorMol();
    }

    void sire_to_openmm(OpenMM::System &system,
                        const SelectorMol &mols,
                        const PropertyMap &map)
    {
        system.addParticle(3.0);
    }

} // end of namespace SireOpenMM
