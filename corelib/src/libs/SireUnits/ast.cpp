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

#include "SireUnits/ast.h"

namespace SireUnits
{
    namespace AST
    {
        QString Unit::toString() const
        {
            return unit.toString();
        }

        GeneralUnit Unit::toUnit() const
        {
            return unit;
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
            GeneralUnit ret = scale * unit.toUnit();

            if (power != 1.0)
            {
                int pow = std::abs(int(power));

                if (pow == 0)
                {
                    return GeneralUnit(1);
                }

                GeneralUnit multiple = ret;

                for (int i = 0; i < pow; ++i)
                {
                    multiple *= ret;
                }

                if (power < 0)
                {
                    return 1.0 / multiple;
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
