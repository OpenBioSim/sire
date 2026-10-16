/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2023  Christopher Woods
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

#ifndef SIREMOL_ISWATER_H
#define SIREMOL_ISWATER_H

#include "SireMol/core.h"
#include "SireMol/select.h"
#include "SireMol/selectormol.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    SIREMOL_EXPORT
    bool is_water(const MoleculeView &molecule,
                  const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SIREMOL_EXPORT
    QVector<bool> is_water(const SelectorMol &molecules,
                           const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SIREMOL_EXPORT
    QVector<bool> is_water(const SelectResult &molecules,
                           const SireBase::PropertyMap &map = SireBase::PropertyMap());
}

SIRE_EXPOSE_FUNCTION(SireMol::is_water)

SIRE_END_HEADER

#endif
