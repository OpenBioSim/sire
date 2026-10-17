/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "SireMM/errors.h"

using namespace SireMM;

const char *missing_bond::typeName()
{
    return QMetaType::typeName(qMetaTypeId<missing_bond>());
}

const char *missing_angle::typeName()
{
    return QMetaType::typeName(qMetaTypeId<missing_angle>());
}

const char *missing_dihedral::typeName()
{
    return QMetaType::typeName(qMetaTypeId<missing_dihedral>());
}

static const RegisterMetaType<missing_bond> r_mbond;
static const RegisterMetaType<missing_angle> r_mang;
static const RegisterMetaType<missing_dihedral> r_mdih;
