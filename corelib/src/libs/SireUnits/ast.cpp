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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-builtins"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "SireUnits/ast.h"

#include <QDebug>

namespace SireUnits
{
    namespace AST
    {
        QString Unit::toString() const
        {
            return this->toUnit().toString();
        }

        GeneralUnit Unit::toUnit() const
        {
            return unit * prefix;
        }

        QString Expression::toString() const
        {
            return this->toUnit().toString();
        }

        GeneralUnit Expression::toUnit() const
        {
            return unit.toUnit();
        }

        QString Node::toString() const
        {
            return this->toUnit().toString();
        }

        GeneralUnit Node::toUnit() const
        {
            return values.toUnit();
        }

        QString FullUnit::toString() const
        {
            return this->toUnit().toString();
        }

        GeneralUnit FullUnit::toUnit() const
        {
            GeneralUnit ret = unit.toUnit();

            if (power != 1)
            {
                int pow = std::abs(power);

                if (pow == 0)
                {
                    return GeneralUnit(1);
                }

                GeneralUnit multiple = ret;

                for (int i = 1; i < pow; ++i)
                {
                    multiple *= ret;
                }

                if (power < 0)
                {
                    return multiple.invert();
                }
                else
                {
                    return multiple;
                }
            }
            else
            {
                return ret;
            }
        }

    } // namespace AST
} // namespace SireUnits

#ifdef __clang__
#pragma clang diagnostic pop
#endif
