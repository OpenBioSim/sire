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

#ifndef SIREUNITS_TEMPERATURE
#define SIREUNITS_TEMPERATURE

#include "dimensions.h"

SIRE_BEGIN_HEADER

namespace SireUnits
{

    class Celsius;
    class Fahrenheit;

    namespace Dimension
    {

        class GeneralUnit;

// skip this completely when parsing with gccxml as it is broken!
#ifdef SKIP_BROKEN_GCCXML_PARTS

        class Temperature
        {
        public:
            Temperature();
            Temperature(double);
            ~Temperature();

            operator double() const;
        };

#endif // end of 'ifdef SKIP_BROKEN_GCCXML_PARTS'

        class SIREUNITS_EXPORT TempBase
        {
            friend class SireUnits::Celsius;
            friend class SireUnits::Fahrenheit;

        public:
            TempBase(double value = 0);
            TempBase(const TempBase &other);
            TempBase(const Temperature &temp);

            virtual ~TempBase();

            TempBase &operator=(const TempBase &other);
            TempBase &operator=(const Temperature &temp);

            bool operator==(const TempBase &other) const;
            bool operator!=(const TempBase &other) const;
            bool operator==(const Temperature &temp) const;
            bool operator!=(const Temperature &temp) const;

            double value() const;

            QString toString() const;

            operator Temperature() const;

            operator double() const;

            double in(const TempBase &other) const;
            double in(const Temperature &temp) const;

            double to(const TempBase &other) const;
            double to(const GeneralUnit &other) const;
            double to(const QString &other) const;

            virtual double convertToInternal(double value) const = 0;
            virtual double convertFromInternal(double value) const = 0;

            double convertFromInternal() const;

        protected:
            virtual QString unitString() const;

            /** This holds the temperature in internal units (K) */
            double val;
        };

    } // end of namespace Dimension

    class SIREUNITS_EXPORT Celsius : public Dimension::TempBase
    {

    public:
        Celsius();
        explicit Celsius(double value);
        Celsius(const Dimension::Temperature &temp);
        Celsius(const Dimension::TempBase &other);
        Celsius(const Celsius &other);

        ~Celsius();

        double convertToInternal(double value) const;

        double convertFromInternal(double value) const;
        double convertFromInternal() const;

        Celsius &operator=(const Celsius &other);
        Celsius &operator=(const Dimension::Temperature &temp);

        Celsius operator-() const;

        Celsius operator+(const Celsius &other) const;
        Celsius operator-(const Celsius &other) const;

        Celsius &operator+=(const Celsius &other);
        Celsius &operator-=(const Celsius &other);

        Celsius operator+(const Dimension::Temperature &other) const;
        Celsius operator-(const Dimension::Temperature &other) const;

        Celsius &operator+=(const Dimension::Temperature &other);
        Celsius &operator-=(const Dimension::Temperature &other);

        Celsius operator*(double value) const;
        Celsius operator/(double value) const;
        Celsius operator*(int value) const;
        Celsius operator/(int value) const;

        Dimension::GeneralUnit operator+(const Dimension::GeneralUnit &other) const;
        Dimension::GeneralUnit operator-(const Dimension::GeneralUnit &other) const;
        Dimension::GeneralUnit operator*(const Dimension::GeneralUnit &other) const;
        Dimension::GeneralUnit operator/(const Dimension::GeneralUnit &other) const;

    protected:
        QString unitString() const;
    };

#ifndef SKIP_BROKEN_GCCXML_PARTS
    SIREUNITS_EXPORT Celsius operator*(double value, const Celsius &temp);
    SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> operator/(double value, const Celsius &temp);
    SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> operator/(int value, const Celsius &temp);
    SIREUNITS_EXPORT Celsius operator*(int value, const Celsius &temp);
#endif

    class SIREUNITS_EXPORT Fahrenheit : public Dimension::TempBase
    {

    public:
        Fahrenheit();
        explicit Fahrenheit(double value);
        Fahrenheit(const Dimension::Temperature &temp);
        Fahrenheit(const Dimension::TempBase &other);
        Fahrenheit(const Fahrenheit &other);

        ~Fahrenheit();

        double convertToInternal(double value) const;

        double convertFromInternal(double value) const;
        double convertFromInternal() const;

        Fahrenheit &operator=(const Fahrenheit &other);
        Fahrenheit &operator=(const Dimension::Temperature &temp);

        Fahrenheit operator-() const;

        Fahrenheit operator+(const Fahrenheit &other) const;
        Fahrenheit operator-(const Fahrenheit &other) const;

        Fahrenheit &operator+=(const Fahrenheit &other);
        Fahrenheit &operator-=(const Fahrenheit &other);

        Fahrenheit operator+(const Dimension::Temperature &other) const;
        Fahrenheit operator-(const Dimension::Temperature &other) const;

        Fahrenheit &operator+=(const Dimension::Temperature &other);
        Fahrenheit &operator-=(const Dimension::Temperature &other);

        Fahrenheit operator*(double value) const;
        Fahrenheit operator/(double value) const;

        Fahrenheit operator*(int value) const;
        Fahrenheit operator/(int value) const;

        Dimension::GeneralUnit operator+(const Dimension::GeneralUnit &other) const;
        Dimension::GeneralUnit operator-(const Dimension::GeneralUnit &other) const;
        Dimension::GeneralUnit operator*(const Dimension::GeneralUnit &other) const;
        Dimension::GeneralUnit operator/(const Dimension::GeneralUnit &other) const;

    protected:
        QString unitString() const;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTIONS
    SIREUNITS_EXPORT Fahrenheit operator*(double value, const Fahrenheit &temp);
    SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> operator/(double value, const Fahrenheit &temp);
    SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> operator/(int value, const Fahrenheit &temp);
    SIREUNITS_EXPORT Fahrenheit operator*(int value, const Fahrenheit &temp);
#endif // SIRE_SKIP_INLINE_FUNCTIONS

    const Celsius celsius(1);
    const Fahrenheit fahrenheit(1);

} // namespace SireUnits

SIRE_EXPOSE_CLASS(SireUnits::Dimension::TempBase)
SIRE_EXPOSE_CLASS(SireUnits::Celsius)
SIRE_EXPOSE_CLASS(SireUnits::Fahrenheit)

SIRE_END_HEADER

#endif
