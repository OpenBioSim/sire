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

#ifndef SIRESYSTEM_FORCEFIELDINFO_H
#define SIRESYSTEM_FORCEFIELDINFO_H

#include "SireBase/property.h"
#include "SireBase/properties.h"
#include "SireBase/propertymap.h"

#include "SireVol/space.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
    class ForceFieldInfo;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireSystem::ForceFieldInfo &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireSystem::ForceFieldInfo &);

namespace SireMol
{
    class MoleculeView;
    class SelectorMol;
}

namespace SireSystem
{
    class System;

    /** This class can generate and hold all of the system-level information
     *  that is needed to specify a forcefield, e.g. space, cutoff type,
     *  cutoff lengths etc.
     */
    class SIRESYSTEM_EXPORT ForceFieldInfo
        : public SireBase::ConcreteProperty<ForceFieldInfo, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const ForceFieldInfo &);
        friend QDataStream & ::operator>>(QDataStream &, ForceFieldInfo &);

    public:
        ForceFieldInfo();
        ForceFieldInfo(const System &system,
                       const SireBase::PropertyMap &map = SireBase::PropertyMap());
        ForceFieldInfo(const SireMol::MoleculeView &mol,
                       const SireBase::PropertyMap &map = SireBase::PropertyMap());
        ForceFieldInfo(const SireMol::SelectorMol &mols,
                       const SireBase::PropertyMap &map = SireBase::PropertyMap());

        ForceFieldInfo(const ForceFieldInfo &other);

        ~ForceFieldInfo();

        ForceFieldInfo &operator=(const ForceFieldInfo &other);

        bool operator==(const ForceFieldInfo &other) const;
        bool operator!=(const ForceFieldInfo &other) const;

        QString toString() const;

        const char *what() const;

        static const char *typeName();

        ForceFieldInfo *clone() const;

        const SireVol::Space &space() const;
        void setSpace(const SireVol::Space &space);

        SireUnits::Dimension::Length cutoff() const;
        void setCutoff(SireUnits::Dimension::Length length);
        void setLargestCutoff();

        bool hasCutoff() const;

        bool hasNoCutoff() const;
        void setNoCutoff();

        QString cutoffType() const;
        void setCutoffType(QString cutoff_type);

        void setCutoffType(QString cutoff_type,
                           const SireBase::PropertyMap &parameters);

        static QStringList cutoffTypes();

        SireUnits::Dimension::GeneralUnit getParameter(const QString &parameter) const;
        void setParameter(const QString &parameter,
                          const SireUnits::Dimension::GeneralUnit &value);

        SireBase::Properties parameters() const;

    private:
        enum CUTOFF_TYPE
        {
            NO_CUTOFF = 0,
            CUTOFF = 1,
            EWALD = 2,
            PME = 3,
            REACTION_FIELD = 4,
            SHIFT_ELECTROSTATICS = 5
        };

        static CUTOFF_TYPE string_to_cutoff_type(QString s);
        static QString cutoff_type_to_string(CUTOFF_TYPE type);

        /** The system space */
        SireVol::SpacePtr spc;

        /** Any extra parameters */
        SireBase::Properties params;

        /** The cutoff - this is zero if there is no cutoff */
        SireUnits::Dimension::Length ctff;

        /** The type of cutoff, e.g. PME etc. */
        CUTOFF_TYPE ctff_typ;
    };

}

Q_DECLARE_METATYPE(SireSystem::ForceFieldInfo)

SIRE_EXPOSE_CLASS(SireSystem::ForceFieldInfo)

SIRE_END_HEADER

#endif
