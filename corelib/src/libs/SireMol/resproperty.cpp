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

#include "resproperty.hpp"

using namespace SireMol;

namespace SireMol
{
    template class ResProperty<QString>;
    template class ResProperty<qint64>;
    template class ResProperty<double>;
    template class ResProperty<QVariant>;
    template class ResProperty<SireBase::PropertyPtr>;
} // namespace SireMol

static const RegisterMetaType<ResStringProperty> r_resstring;
static const RegisterMetaType<ResIntProperty> r_resint;
static const RegisterMetaType<ResFloatProperty> r_resfloat;
static const RegisterMetaType<ResVariantProperty> r_resvariant;
static const RegisterMetaType<ResPropertyProperty> r_respropprop;
