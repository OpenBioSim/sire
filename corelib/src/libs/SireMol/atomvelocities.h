/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREMOL_ATOMVELOCITIES_H
#define SIREMOL_ATOMVELOCITIES_H

#include "atomproperty.hpp"

#include "SireMaths/vector3d.hpp"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

    using SireUnits::Dimension::Velocity;

    typedef SireMaths::Vector3D<Velocity> Velocity3D;

    typedef AtomProperty<Velocity3D> AtomVelocities;

} // namespace SireMol

SIRE_EXPOSE_ALIAS((SireMaths::Vector3D<SireUnits::Dimension::PhysUnit<0, 1, -1, 0, 0, 0, 0>>), SireMol::Velocity3D)

Q_DECLARE_METATYPE(SireMol::AtomVelocities);
Q_DECLARE_METATYPE(SireMol::Velocity3D);

SIRE_EXPOSE_ATOM_PROPERTY(SireMaths::Vector3D<SireUnits::Dimension::Velocity>, SireMol::AtomVelocities)

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMaths::Vector3D<SireUnits::Dimension::Velocity>;
template class SireMol::AtomProperty<SireMol::Velocity3D>;
#endif

SIRE_END_HEADER

#endif
