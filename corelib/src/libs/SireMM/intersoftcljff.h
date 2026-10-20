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

#ifndef SIREMM_INTERSOFTCLJFF_H
#define SIREMM_INTERSOFTCLJFF_H

#include "softcljpotential.h"

#include "SireFF/inter2b2g3dff.hpp"
#include "SireFF/inter2b3dff.hpp"

SIRE_BEGIN_HEADER

namespace SireMM
{

    using SireFF::Inter2B2G3DFF;
    using SireFF::Inter2B2GFF;
    using SireFF::Inter2B3DFF;
    using SireFF::Inter2BFF;

    typedef Inter2BFF<SoftCLJPotentialInterface<InterSoftCLJPotential>> InterSoftCLJFFBase;
    typedef Inter2B3DFF<SoftCLJPotentialInterface<InterSoftCLJPotential>> InterSoftCLJFF;

    typedef Inter2B2GFF<SoftCLJPotentialInterface<InterSoftCLJPotential>> InterGroupSoftCLJFFBase;

    typedef Inter2B2G3DFF<SoftCLJPotentialInterface<InterSoftCLJPotential>> InterGroupSoftCLJFF;

} // namespace SireMM

Q_DECLARE_METATYPE(SireMM::InterSoftCLJFFBase);
Q_DECLARE_METATYPE(SireMM::InterSoftCLJFF);

Q_DECLARE_METATYPE(SireMM::InterGroupSoftCLJFFBase);
Q_DECLARE_METATYPE(SireMM::InterGroupSoftCLJFF);

SIRE_EXPOSE_ALIAS(SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential>,
                  SireMM::CLJPotentialInterface_InterSoftCLJPotential_)

SIRE_EXPOSE_ALIAS(SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>,
                  SireMM::SoftCLJPotentialInterface_InterSoftCLJPotential_)

SIRE_EXPOSE_ALIAS(SireFF::Inter2BFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>,
                  SireMM::InterSoftCLJFFBase)

SIRE_EXPOSE_ALIAS(SireFF::Inter2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>,
                  SireMM::InterSoftCLJFF)

SIRE_EXPOSE_ALIAS(SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>,
                  SireMM::InterGroupSoftCLJFFBase)

SIRE_EXPOSE_ALIAS(SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>,
                  SireMM::InterGroupSoftCLJFF)

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>;

template class SireFF::Inter2BFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>;

template class SireFF::Inter2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>;

template class SireFF::Inter2B2GFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>;

template class SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>>;
#endif

SIRE_END_HEADER

#endif
