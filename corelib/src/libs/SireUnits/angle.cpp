/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "angle.h"
#include "units.h"

using namespace SireUnits;
using namespace SireUnits::Dimension;

Angle SireUnits::acute(Angle angle)
{
    const Angle three_sixty = 360 * degrees;
    const Angle one_eighty = 180 * degrees;

    // unwind the angle...
    while (angle < 0)
    {
        angle += three_sixty;
    }

    while (angle >= three_sixty)
    {
        angle -= three_sixty;
    }

    if (angle > one_eighty)
        return three_sixty - angle;
    else
        return angle;
}
