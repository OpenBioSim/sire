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

#ifndef SIREMM_RESTRAINTS_H
#define SIREMM_RESTRAINTS_H

#include "SireBase/property.h"

namespace SireMM
{
    class Restraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &ds, const SireMM::Restraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &ds, SireMM::Restraints &);

namespace SireMM
{
    /** This is the base class of all of the `Restraints` collections,
     *  e.g. PositionalRestraints, BoreschRestraints etc.
     *
     *  These are not related to the legacy SireMM::Restraint classes,
     *  which calculate restraints directly. Instead, these classes
     *  provide information about the restraints that can be picked
     *  up and used by different engines (e.g. OpenMM)
     */
    class SIREMM_EXPORT Restraints : public SireBase::Property
    {
        friend QDataStream & ::operator<<(QDataStream &, const Restraints &);
        friend QDataStream & ::operator>>(QDataStream &, Restraints &);

    public:
        Restraints();
        Restraints(const QString &name);

        Restraints(const Restraints &other);

        virtual ~Restraints();

        virtual Restraints *clone() const = 0;

        static const char *typeName()
        {
            return "SireMM::Restraints";
        }

        QString name() const;

        void setName(const QString &name);

        bool isRestraints() const;

    protected:
        Restraints &operator=(const Restraints &other);

        bool operator==(const Restraints &other) const;
        bool operator!=(const Restraints &other) const;

    private:
        /** The name of this group of restraints */
        QString nme;
    };
}

SIRE_EXPOSE_CLASS(SireMM::Restraints)

#endif
