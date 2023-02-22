/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMOL_BONDORDER_H
#define SIREMOL_BONDORDER_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class BondOrder;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::BondOrder &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::BondOrder &);

namespace SireMol
{

    /** This class represents a bond type (e.g. single, double etc.)

        @author Christopher Woods
    */
    class SIREMOL_EXPORT BondOrder : public SireBase::ConcreteProperty<BondOrder, SireBase::Property>
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const BondOrder &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, BondOrder &);

    public:
        BondOrder();

        BondOrder(const QString &s);

        static BondOrder fromSDF(int value);
        static BondOrder fromRDKit(const QString &value);

        BondOrder(const BondOrder &other);

        ~BondOrder();

        BondOrder &operator=(const BondOrder &other);

        bool operator==(const BondOrder &other) const;
        bool operator!=(const BondOrder &other) const;

        static const char *typeName();

        static BondOrder singleBond();
        static BondOrder doubleBond();
        static BondOrder tripleBond();
        static BondOrder aromaticBond();
        static BondOrder undefinedBond();

        QString toString() const;

        int value() const;

        double valueAsDouble() const;

        int toSDF() const;
        QString toRDKit() const;

        bool isDefined() const;
        bool isSingle() const;
        bool isDouble() const;
        bool isTriple() const;
        bool isAromatic() const;

    private:
        /** The bond type. We use an integer in SDF format,
         *  plus some RDKit additions...
         */
        qint32 bond_type;
    };

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::BondOrder)

SIRE_EXPOSE_CLASS(SireMol::BondOrder)

SIRE_END_HEADER

#endif
