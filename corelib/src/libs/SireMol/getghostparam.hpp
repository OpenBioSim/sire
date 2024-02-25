/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREMOL_GETGHOSTPARAM_HPP
#define SIREMOL_GETGHOSTPARAM_HPP

#include "SireBase/booleanproperty.h"
#include "SireBase/stringproperty.h"

#include "SireMol/element.h"

#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    template <class T>
    inline T getGhostParam(const QString &ghost)
    {
        return T();
    }

    template <>
    inline QString getGhostParam(const QString &ghost)
    {
        return ghost;
    }

    template <>
    inline qint64 getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return 0;
        else
            return ghost.toLongLong();
    }

    template <>
    inline double getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return 0.0;
        else
            return ghost.toDouble();
    }

    template <>
    inline bool getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return false;
        else
            return SireBase::BooleanProperty(ghost).value();
    }

    template <>
    inline SireUnits::Dimension::GeneralUnit getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireUnits::Dimension::GeneralUnit(0);
        else
            return SireUnits::Dimension::GeneralUnit(ghost);
    }

    template <>
    inline QVariant getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return QVariant();
        else
            return QVariant(ghost);
    }

    template <>
    inline SireUnits::Dimension::MolarMass getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireUnits::Dimension::MolarMass(0.0);
        else
            return SireUnits::Dimension::MolarMass(ghost);
    }

    template <>
    inline SireUnits::Dimension::MolarEnergy getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireUnits::Dimension::MolarEnergy(0.0);
        else
            return SireUnits::Dimension::MolarEnergy(ghost);
    }

    template <>
    inline SireUnits::Dimension::Charge getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireUnits::Dimension::Charge(0.0);
        else
            return SireUnits::Dimension::Charge(ghost);
    }

    template <>
    inline SireMol::Element getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireMol::Element(0);
        else
            return SireMol::Element(ghost);
    }

    template <>
    inline SireUnits::Dimension::Volume getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireUnits::Dimension::Volume(0.0);
        else
            return SireUnits::Dimension::Volume(ghost);
    }

    template <>
    SireUnits::Dimension::Length getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireUnits::Dimension::Length(0.0);
        else
            return SireUnits::Dimension::Length(ghost);
    }

    template <>
    SireBase::PropertyPtr getGhostParam(const QString &ghost)
    {
        if (ghost.isEmpty())
            return SireBase::PropertyPtr();
        else
            return SireBase::PropertyPtr(SireBase::StringProperty(ghost));
    }

}

SIRE_END_HEADER

#endif
