/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2023  Christopher Woods
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
  *  at http://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREMOL_GETRMSD_H
#define SIREMOL_GETRMSD_H

#include "core.h"
#include "trajectoryaligner.h"
#include "selectorm.hpp"
#include "atom.h"

#include "SireMaths/vector.h"
#include "SireBase/propertymap.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    SIREMOL_EXPORT QVector<SireUnits::Dimension::Length>
    get_rmsd(const SelectorM<Atom> &atoms,
             const QVector<SireMaths::Vector> &coords,
             const TrajectoryAligner &aligner,
             const QList<qint64> &frames,
             const PropertyMap &map);
}

SIRE_EXPOSE_FUNCTION(SireMol::get_rmsd)

SIRE_END_HEADER

#endif
