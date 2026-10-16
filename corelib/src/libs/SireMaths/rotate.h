/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMATHS_ROTATE_H
#define SIREMATHS_ROTATE_H

#include "matrix.h"
#include "vector.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

    /** Function that rotates the input vector 'input' by the rotation
        matrix 'rotmat' about the point 'point'. The output is returned. */
    SIRE_ALWAYS_INLINE Vector rotate(const Vector &input, const Matrix &rotmat, const Vector &point)
    {
        return point + rotmat * (input - point);
    }

} // namespace SireMaths

SIRE_EXPOSE_FUNCTION(SireMaths::rotate)

SIRE_END_HEADER

#endif
