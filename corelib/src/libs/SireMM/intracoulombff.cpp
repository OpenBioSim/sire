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

#include "intracoulombff.h"

#include "SireMol/partialmolecule.h"

#include "SireMol/mover.hpp"

using namespace SireMM;
using namespace SireFF;

namespace SireMM
{
    template class CoulombPotentialInterface<IntraCoulombPotential>;
}

namespace SireFF
{
    template class Intra2BFF<CoulombPotentialInterface<IntraCoulombPotential>>;
    template class Intra2B3DFF<CoulombPotentialInterface<IntraCoulombPotential>>;

    template class Intra2B2GFF<CoulombPotentialInterface<IntraCoulombPotential>>;
    template class Intra2B2G3DFF<CoulombPotentialInterface<IntraCoulombPotential>>;
} // namespace SireFF

static const RegisterMetaType<IntraCoulombFF> r_intracoulff;
static const RegisterMetaType<IntraGroupCoulombFF> r_intragroupcoulff;
