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

#ifndef SIREMOL_CGIDX_H
#define SIREMOL_CGIDX_H

#include "SireID/index.h"

#include "cgid.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class CGIdx;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::CGIdx &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::CGIdx &);

namespace SireMol
{

    class CGAtomIdx;

    /** This is an ID object that is used to index CutGroups

        @author Christopher Woods
    */
    class SIREMOL_EXPORT CGIdx : public SireID::Index_T_<CGIdx>, public CGID
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const CGIdx &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, CGIdx &);

    public:
        CGIdx();

        explicit CGIdx(qint32 idx);

        CGIdx(const CGIdx &other);

        ~CGIdx();

        static const char *typeName();

        const char *what() const
        {
            return SireID::Index_T_<CGIdx>::what();
        }

        CGIdx *clone() const;

        static CGIdx null();

        using CGID::operator+;

        CGAtomIdx operator+(const SireID::Index &other) const;

        bool isNull() const;

        uint hash() const;

        QString toString() const;

        CGIdx &operator=(const CGIdx &other);

        bool operator==(const SireID::ID &other) const;

        bool operator!=(const SireID::ID &other) const
        {
            return not this->operator==(other);
        }

        using SireID::Index_T_<CGIdx>::operator=;

        using SireID::Index_T_<CGIdx>::operator==;
        using SireID::Index_T_<CGIdx>::operator!=;

        using SireID::Index_T_<CGIdx>::operator+=;
        using SireID::Index_T_<CGIdx>::operator++;
        using SireID::Index_T_<CGIdx>::operator-=;
        using SireID::Index_T_<CGIdx>::operator--;

        using SireID::Index_T_<CGIdx>::map;

        QList<CGIdx> map(const MolInfo &molinfo) const;
    };

    SIREMOL_EXPORT CGAtomIdx operator+(const SireID::Index &index, const CGIdx &cgidx);

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::CGIdx);

SIRE_EXPOSE_CLASS(SireMol::CGIdx)

SIRE_END_HEADER

#endif
