/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMOL_MGIDENTIFIER_H
#define SIREMOL_MGIDENTIFIER_H

#include "mgid.h"

#include <memory>

namespace SireMol
{
    class MGIdentifier;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::MGIdentifier &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::MGIdentifier &);

namespace SireMol
{

    /** This is a generic holder for any MGID class!

        @author Christopher Woods
    */
    class SIREMOL_EXPORT MGIdentifier : public MGID
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const MGIdentifier &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, MGIdentifier &);

    public:
        MGIdentifier();
        MGIdentifier(const MGID &atomid);
        MGIdentifier(const MGIdentifier &other);

        ~MGIdentifier();

        static const char *typeName();

        const char *what() const
        {
            return MGIdentifier::typeName();
        }

        MGIdentifier *clone() const;

        bool isNull() const;

        uint hash() const;

        QString toString() const;

        const MGID &base() const;

        MGIdentifier &operator=(const MGIdentifier &other);
        MGIdentifier &operator=(const MGID &other);

        bool operator==(const SireID::ID &other) const;

        bool operator!=(const SireID::ID &other) const
        {
            return not this->operator==(other);
        }

        bool operator==(const MGIdentifier &other) const;
        bool operator!=(const MGIdentifier &other) const;

        bool operator==(const MGID &other) const;
        bool operator!=(const MGID &other) const;

        QList<MGNum> map(const MolGroupsBase &molgroups) const;

    private:
        /** Pointer to the MGID */
        std::shared_ptr<MGID> d;
    };

    SIRE_ALWAYS_INLINE uint qHash(const MGIdentifier &mgid)
    {
        return mgid.hash();
    }

} // namespace SireMol

#include "mgnum.h"

Q_DECLARE_METATYPE(SireID::IDAndSet<SireMol::MGID>)
Q_DECLARE_METATYPE(SireID::IDOrSet<SireMol::MGID>)
Q_DECLARE_METATYPE(SireID::Specify<SireMol::MGID>)

Q_DECLARE_METATYPE(SireMol::MGIdentifier);

#endif
