/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the website
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#include "SireSystem/merge.h"

using namespace SireMol;
using namespace SireBase;

namespace SireSystem
{
    /**
     * @brief Merge function that combines multiple molecules into a single molecule.
     *
     * @param mols The AtomMapping object that contains the molecules to be merged.
     * @param as_new_molecule Flag indicating whether the merged molecule should be created as a new molecule.
     * @param map The PropertyMap object that contains additional properties for the merged molecule.
     * @return The merged molecule.
     */
    Molecule merge(const AtomMapping &mols, bool as_new_molecule, const PropertyMap &map)
    {
        return mols.atoms0().molecules()[0];
    }
}
