/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREIO_GETCOORDSARRAY_H
#define SIREIO_GETCOORDSARRAY_H

#include "SireBase/propertymap.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/moleculeview.h"
#include "SireSystem/system.h"
#include "SireUnits/dimensions.h"

#include <QVector>

SIRE_BEGIN_HEADER

namespace SireIO
{
    SIREIO_EXPORT QVector<float> getCoordsArray(const SireMol::MoleculeView &mol,
                                                const SireUnits::Dimension::Length &to_unit,
                                                const SireBase::PropertyMap &map);

    SIREIO_EXPORT QVector<float> getCoordsArray(const SireMol::MoleculeGroup &mols,
                                                const SireUnits::Dimension::Length &to_unit,
                                                const SireBase::PropertyMap &map);

    SIREIO_EXPORT QVector<float> getCoordsArray(const SireSystem::System &system,
                                                const SireUnits::Dimension::Length &to_unit,
                                                const SireBase::PropertyMap &map);
} // namespace SireIO

SIRE_EXPOSE_FUNCTION(SireIO::getCoordsArray)

SIRE_END_HEADER

#endif
