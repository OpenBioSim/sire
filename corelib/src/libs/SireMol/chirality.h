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

#ifndef SIREMOL_CHIRALITY_H
#define SIREMOL_CHIRALITY_H

#include "SireBase/property.h"

#include "SireMol/atomproperty.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class Chirality;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::Chirality &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::Chirality &);

namespace SireMol
{

    /** This class represents an atom's chirality

        @author Christopher Woods
    */
    class SIREMOL_EXPORT Chirality : public SireBase::ConcreteProperty<Chirality, SireBase::Property>
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const Chirality &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, Chirality &);

    public:
        Chirality();

        Chirality(const Chirality &other);

        ~Chirality();

        static Chirality fromSDF(int val);
        static Chirality fromRDKit(const QString &value);

        Chirality &operator=(const Chirality &other);

        bool operator==(const Chirality &other) const;
        bool operator!=(const Chirality &other) const;

        static const char *typeName();

        static Chirality clockwise();
        static Chirality counterClockwise();
        static Chirality other();
        static Chirality undefined();

        QString toString() const;

        int toSDF() const;
        QString toRDKit() const;

        bool isClockwise() const;
        bool isCounterClockwise() const;
        bool isOther() const;
        bool isUndefined() const;

    private:
        /** The chiral type */
        qint32 chiral_type;
    };

    typedef AtomProperty<Chirality> AtomChiralities;

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::Chirality)
Q_DECLARE_METATYPE(SireMol::AtomChiralities);

SIRE_EXPOSE_CLASS(SireMol::Chirality)

SIRE_EXPOSE_ATOM_PROPERTY(SireMol::Chirality, SireMol::AtomChiralities)

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::AtomProperty<SireMol::Chirality>;
#endif

SIRE_END_HEADER

#endif
