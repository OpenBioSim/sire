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
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREMM_LJ1264PARAMETER_H
#define SIREMM_LJ1264PARAMETER_H

#include "ljparameter.h"
#include "ljpair.h"

#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class LJ1264Parameter;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::LJ1264Parameter &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::LJ1264Parameter &);

namespace SireMM
{
    /** This class holds a LJ 12-6-4 parameter. These are used
     *  between some atom-pairs to define a 12-6-4 potential
     *  with functional form
     *
     *  E_LJ = A / r^12 - B / r^6 - C / r^4
     */
    class SIREMM_EXPORT LJ1264Parameter
    {
        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const LJ1264Parameter &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, LJ1264Parameter &);

    public:
        LJ1264Parameter();
        LJ1264Parameter(double a, double b, double c);

        LJ1264Parameter(const SireUnits::Dimension::GeneralUnit &a,
                        const SireUnits::Dimension::GeneralUnit &b,
                        const SireUnits::Dimension::GeneralUnit &c);

        LJ1264Parameter(const LJParameter &ljparam);
        LJ1264Parameter(const LJPair &ljpair);

        LJ1264Parameter(const LJ1264Parameter &other);

        ~LJ1264Parameter();

        LJ1264Parameter &operator=(const LJ1264Parameter &other);

        bool operator==(const LJ1264Parameter &other) const;
        bool operator!=(const LJ1264Parameter &other) const;

        static const char *typeName();

        const char *what() const;

        LJ1264Parameter *clone() const;

        bool isDummy() const;
        bool zeroLJ() const;

        double A() const;
        double B() const;
        double C() const;

        static SireUnits::Dimension::GeneralUnit AUnit();
        static SireUnits::Dimension::GeneralUnit BUnit();
        static SireUnits::Dimension::GeneralUnit CUnit();

        QString toString() const;

        static LJ1264Parameter dummy();

        bool hasC() const;
        bool isLJParameter() const;

        void assertIsLJParameter() const;

        LJParameter toLJParameter() const;
        LJPair toLJPair() const;

    private:
        double a, b, c;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

    /** Hash a LJ parameter */
    SIRE_ALWAYS_INLINE uint qHash(const LJ1264Parameter &ljparam)
    {
        return uint(1000000.0 * ljparam.A() + 10000.0 * ljparam.B() + 10.0 * ljparam.C());
    }

    /** Return whether or not two LJ1264Parameters are equal */
    SIRE_ALWAYS_INLINE bool LJ1264Parameter::operator==(const LJ1264Parameter &other) const
    {
        return a == other.a and b == other.b and c == other.c;
    }

    /** Return whether or not two LJ1264Parameters are different */
    SIRE_ALWAYS_INLINE bool LJ1264Parameter::operator!=(const LJ1264Parameter &other) const
    {
        return not operator==(other);
    }

    /** Return whether or not this is a dummy LJ parameter */
    SIRE_ALWAYS_INLINE bool LJ1264Parameter::isDummy() const
    {
        return a == 0.0 and b == 0.0 and c == 0.0;
    }

    /** Return whether or not this is a dummy LJ parameter */
    SIRE_ALWAYS_INLINE bool LJ1264Parameter::zeroLJ() const
    {
        return isDummy();
    }

    /** Return the A parameter */
    SIRE_ALWAYS_INLINE double LJ1264Parameter::A() const
    {
        return a;
    }

    /** Return the B parameter */
    SIRE_ALWAYS_INLINE double LJ1264Parameter::B() const
    {
        return b;
    }

    /** Return the C parameter */
    SIRE_ALWAYS_INLINE double LJ1264Parameter::C() const
    {
        return c;
    }

#endif

} // namespace SireMM

Q_DECLARE_TYPEINFO(SireMM::LJ1264Parameter, Q_MOVABLE_TYPE);
Q_DECLARE_METATYPE(SireMM::LJ1264Parameter);

SIRE_EXPOSE_CLASS(SireMM::LJ1264Parameter)

SIRE_END_HEADER

#endif
