/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
  *
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published
  *  by the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  *
  *  For full details of the license please see the LICENSE file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the website
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIRESYSTEM_MERGE_H
#define SIRESYSTEM_MERGE_H

#include "SireMol/atommapping.h"
#include "SireMol/molecule.h"

#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
    SIRESYSTEM_EXPORT SireMol::Molecule merge(const SireMol::AtomMapping &mols,
                                              const QStringList &properties = QStringList(),
                                              bool as_new_molecule = true,
                                              bool allow_ring_breaking = false,
                                              bool allow_ring_size_change = false,
                                              bool force = false,
                                              const SireBase::PropertyMap &map = SireBase::PropertyMap());
}

SIRE_EXPOSE_FUNCTION(SireSystem::merge);

SIRE_END_HEADER;

#endif
