/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#ifndef SIREMOL_BEADID_H
#define SIREMOL_BEADID_H

#include "SireID/id.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

    /** This is the base class of all identifiers that are used
        to identify a Bead

        @author Christopher Woods
    */
    class SIREMOL_EXPORT BeadID : public SireID::ID
    {

    public:
        BeadID();
        BeadID(const BeadID &other);

        virtual ~BeadID();

        static const char *typeName()
        {
            return "SireMol::BeadID";
        }

        virtual BeadID *clone() const = 0;
    };

} // namespace SireMol

SIRE_EXPOSE_CLASS(SireMol::BeadID)

SIRE_END_HEADER

#endif
