/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2025  Christopher Woods
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

#ifndef SIREMOL_CMAPID_H
#define SIREMOL_CMAPID_H

#include "atomidentifier.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
    class CMAPID;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::CMAPID &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::CMAPID &);

namespace SireMol
{

    class MoleculeData;
    class MoleculeInfoData;
    class AtomIdx;

    using boost::tuple;

    /** This class provides a generic ID for a cmap between
        five atoms

        @author Christopher Woods
    */
    class SIREMOL_EXPORT CMAPID : public SireID::ID
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const CMAPID &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, CMAPID &);

    public:
        CMAPID();
        CMAPID(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2,
               const AtomID &atom3, const AtomID &atom4);

        CMAPID(const CMAPID &other);

        ~CMAPID();

        static const char *typeName();

        const char *what() const
        {
            return CMAPID::typeName();
        }

        CMAPID *clone() const;

        uint hash() const;

        QString toString() const;

        bool isNull() const;

        CMAPID &operator=(const CMAPID &other);

        bool operator==(const SireID::ID &other) const;

        bool operator!=(const SireID::ID &other) const
        {
            return not this->operator==(other);
        }

        bool operator==(const CMAPID &other) const;
        bool operator!=(const CMAPID &other) const;

        CMAPID mirror() const;

        const AtomID &operator[](int i) const;

        tuple<AtomIdx, AtomIdx, AtomIdx, AtomIdx, AtomIdx> map(const MoleculeInfoData &molinfo) const;

        tuple<AtomIdx, AtomIdx, AtomIdx, AtomIdx, AtomIdx> map(const MoleculeInfoData &mol0info,
                                                               const MoleculeInfoData &mol1info,
                                                               const MoleculeInfoData &mol2info,
                                                               const MoleculeInfoData &mol3info,
                                                               const MoleculeInfoData &mol4info) const;

        const AtomID &atom0() const;
        const AtomID &atom1() const;
        const AtomID &atom2() const;
        const AtomID &atom3() const;
        const AtomID &atom4() const;

    private:
        /** The identifiers of the five atoms */
        AtomIdentifier atm0, atm1, atm2, atm3, atm4;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

    SIRE_ALWAYS_INLINE uint qHash(const CMAPID &CMAPID)
    {
        return CMAPID.hash();
    }

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::CMAPID);

SIRE_EXPOSE_CLASS(SireMol::CMAPID)

SIRE_END_HEADER

#endif
