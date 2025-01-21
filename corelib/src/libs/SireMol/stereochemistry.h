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

#ifndef SIREMOL_STEREOCHEMISTRY_H
#define SIREMOL_STEREOCHEMISTRY_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class Stereochemistry;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::Stereochemistry &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::Stereochemistry &);

namespace SireMol
{

    /** This class represents a bond's stereochemistry

        @author Christopher Woods
    */
    class SIREMOL_EXPORT Stereochemistry : public SireBase::ConcreteProperty<Stereochemistry, SireBase::Property>
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const Stereochemistry &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, Stereochemistry &);

    public:
        Stereochemistry();

        Stereochemistry(const QString &s);

        Stereochemistry(const Stereochemistry &other);

        ~Stereochemistry();

        static Stereochemistry fromSDF(int val);
        static Stereochemistry fromRDKit(const QString &value);

        Stereochemistry &operator=(const Stereochemistry &other);

        bool operator==(const Stereochemistry &other) const;
        bool operator!=(const Stereochemistry &other) const;

        static const char *typeName();

        static Stereochemistry up();
        static Stereochemistry cisOrTrans();
        static Stereochemistry down();
        static Stereochemistry notStereo();
        static Stereochemistry undefined();

        QString toString() const;

        int value() const;
        int toSDF() const;
        QString toRDKit() const;

        bool isDefined() const;
        bool isUp() const;
        bool isCisOrTrans() const;
        bool isDown() const;
        bool isNotStereo() const;

    private:
        /** The stereo type. We use an integer in SDF format */
        qint32 stereo_type;
    };

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::Stereochemistry)

SIRE_EXPOSE_CLASS(SireMol::Stereochemistry)

SIRE_END_HEADER

#endif
