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

#ifndef SIREMOL_HYBRIDIZATION_H
#define SIREMOL_HYBRIDIZATION_H

#include "SireBase/property.h"

#include "SireMol/atomproperty.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class Hybridization;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::Hybridization &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::Hybridization &);

namespace SireMol
{

    /** This class represents an atom's hybridization (e.g. SP2, SP3 etc.)

        @author Christopher Woods
    */
    class SIREMOL_EXPORT Hybridization : public SireBase::ConcreteProperty<Hybridization, SireBase::Property>
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const Hybridization &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, Hybridization &);

    public:
        Hybridization();

        Hybridization(const Hybridization &other);

        ~Hybridization();

        static Hybridization fromSDF(int val);
        static Hybridization fromRDKit(const QString &value);

        Hybridization &operator=(const Hybridization &other);

        bool operator==(const Hybridization &other) const;
        bool operator!=(const Hybridization &other) const;

        static const char *typeName();

        static Hybridization s();
        static Hybridization sp();
        static Hybridization sp2();
        static Hybridization sp3();
        static Hybridization other();
        static Hybridization unspecified();
        static Hybridization unknown();

        QString toString() const;

        int toSDF() const;
        QString toRDKit() const;

        bool is_s() const;
        bool is_sp() const;
        bool is_sp2() const;
        bool is_sp3() const;
        bool isOther() const;
        bool isUnspecified() const;
        bool isUnknown() const;

    private:
        /** The hybridization type */
        qint32 hybrid_type;
    };

    typedef AtomProperty<Hybridization> AtomHybridizations;

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::Hybridization)
Q_DECLARE_METATYPE(SireMol::AtomHybridizations);

SIRE_EXPOSE_CLASS(SireMol::Hybridization)

SIRE_EXPOSE_ATOM_PROPERTY(SireMol::Hybridization, SireMol::AtomHybridizations)

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireMol::AtomProperty<SireMol::Hybridization>;
#endif

SIRE_END_HEADER

#endif
