/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "intrasoftcljff.h"

#include "SireMol/partialmolecule.h"

#include "SireMol/mover.hpp"

using namespace SireMM;
using namespace SireFF;

namespace SireMM
{
    template class SoftCLJPotentialInterface<IntraSoftCLJPotential>;
}

namespace SireFF
{
    template class Intra2BFF<SoftCLJPotentialInterface<IntraSoftCLJPotential>>;
    template class Intra2B3DFF<SoftCLJPotentialInterface<IntraSoftCLJPotential>>;

    template class Intra2B2GFF<SoftCLJPotentialInterface<IntraSoftCLJPotential>>;
    template class Intra2B2G3DFF<SoftCLJPotentialInterface<IntraSoftCLJPotential>>;
} // namespace SireFF

static const RegisterMetaType<IntraSoftCLJFF> r_intrasoftcljff;
static const RegisterMetaType<IntraGroupSoftCLJFF> r_intragroupsoftcljff;
