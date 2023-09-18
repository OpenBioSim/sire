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

#include "temperature.h"
#include "generalunit.h"

using namespace SireUnits;
using namespace SireUnits::Dimension;

////////
//////// Implementation of TempBase
////////

TempBase::TempBase(double value) : val(value)
{
}

TempBase::TempBase(const TempBase &other) : val(other.val)
{
}

TempBase::TempBase(const Temperature &temp) : val(temp)
{
}

TempBase::~TempBase()
{
}

TempBase &TempBase::operator=(const TempBase &other)
{
    val = other.val;
    return *this;
}

TempBase &TempBase::operator=(const Temperature &temp)
{
    val = double(temp);
    return *this;
}

bool TempBase::operator==(const TempBase &other) const
{
    return val == other.val;
}

bool TempBase::operator!=(const TempBase &other) const
{
    return val != other.val;
}

bool TempBase::operator==(const Temperature &temp) const
{
    return val == double(temp);
}

bool TempBase::operator!=(const Temperature &temp) const
{
    return val != double(temp);
}

double TempBase::value() const
{
    return val;
}

QString TempBase::toString() const
{
    return QString("%1°%2").arg(this->convertFromInternal()).arg(this->unitString());
}

/** Convert this into a temperature object */
TempBase::operator Temperature() const
{
    return Temperature(val);
}

TempBase::operator double() const
{
    return val;
}

double TempBase::in(const TempBase &other) const
{
    return other.convertFromInternal(val) / other.convertFromInternal();
}

double TempBase::in(const Temperature &temp) const
{
    return val * temp;
}

double TempBase::to(const TempBase &other) const
{
    return this->in(other);
}

double TempBase::to(const GeneralUnit &other) const
{
    return GeneralUnit(*this).to(other);
}

double TempBase::to(const QString &other) const
{
    return GeneralUnit(*this).to(other);
}

double TempBase::convertFromInternal() const
{
    return this->convertFromInternal(val);
}

QString TempBase::unitString() const
{
    return "K";
}

//////
////// Implementation of Celsius
//////

/** Construct a Unit from a TempBase */
Unit::Unit(const TempBase &temperature) : sclfac(temperature)
{
}

Celsius::Celsius() : Dimension::TempBase(1)
{
}

Celsius::Celsius(double value) : Dimension::TempBase()
{
    val = convertToInternal(value);
}

Celsius::Celsius(const Dimension::Temperature &temp) : Dimension::TempBase(temp)
{
}

Celsius::Celsius(const Dimension::TempBase &other) : Dimension::TempBase(other)
{
}

Celsius::Celsius(const Celsius &other) : Dimension::TempBase(other)
{
}

Celsius::~Celsius()
{
}

double Celsius::convertToInternal(double value) const
{
    return value + 273.15;
}

double Celsius::convertFromInternal(double value) const
{
    return value - 273.15;
}

double Celsius::convertFromInternal() const
{
    return Dimension::TempBase::convertFromInternal();
}

Celsius &Celsius::operator=(const Celsius &other)
{
    Dimension::TempBase::operator=(other);
    return *this;
}

Celsius &Celsius::operator=(const Dimension::Temperature &temp)
{
    Dimension::TempBase::operator=(temp);
    return *this;
}

Celsius Celsius::operator-() const
{
    return Celsius(-convertFromInternal());
}

Celsius Celsius::operator+(const Celsius &other) const
{
    return Celsius(convertFromInternal() + other.convertFromInternal());
}

Celsius Celsius::operator-(const Celsius &other) const
{
    return Celsius(convertFromInternal() - other.convertFromInternal());
}

Celsius &Celsius::operator+=(const Celsius &other)
{
    convertToInternal(convertFromInternal() + other.convertFromInternal());
    return *this;
}

Celsius &Celsius::operator-=(const Celsius &other)
{
    convertToInternal(convertFromInternal() - other.convertFromInternal());
    return *this;
}

Celsius Celsius::operator+(const Dimension::Temperature &other) const
{
    return *this + Celsius(other);
}

Celsius Celsius::operator-(const Dimension::Temperature &other) const
{
    return *this - Celsius(other);
}

Celsius &Celsius::operator+=(const Dimension::Temperature &other)
{
    return this->operator+=(Celsius(other));
}

Celsius &Celsius::operator-=(const Dimension::Temperature &other)
{
    return this->operator-=(Celsius(other));
}

Celsius Celsius::operator*(double value) const
{
    return Celsius(value * convertFromInternal());
}

Celsius Celsius::operator/(double value) const
{
    return Celsius(value / convertFromInternal());
}

Celsius Celsius::operator*(int value) const
{
    return Celsius(value * convertFromInternal());
}

Celsius Celsius::operator/(int value) const
{
    return Celsius(value / convertFromInternal());
}

QString Celsius::unitString() const
{
    return "C";
}

SIREUNITS_EXPORT Celsius SireUnits::operator*(double value, const Celsius &temp)
{
    return temp * value;
}

SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> SireUnits::operator/(double value, const Celsius &temp)
{
    return Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0>(value / temp.convertFromInternal());
}

SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> SireUnits::operator/(int value, const Celsius &temp)
{
    return Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0>(value / temp.convertFromInternal());
}

SIREUNITS_EXPORT Celsius SireUnits::operator*(int value, const Celsius &temp)
{
    return temp * value;
}

///////
/////// Implementation of Fahrenheit
///////

Fahrenheit::Fahrenheit() : Dimension::TempBase(1)
{
}

Fahrenheit::Fahrenheit(double value) : Dimension::TempBase()
{
    val = convertToInternal(value);
}

Fahrenheit::Fahrenheit(const Dimension::Temperature &temp) : Dimension::TempBase(temp)
{
}

Fahrenheit::Fahrenheit(const Dimension::TempBase &other) : Dimension::TempBase(other)
{
}

Fahrenheit::Fahrenheit(const Fahrenheit &other) : Dimension::TempBase(other)
{
}

Fahrenheit::~Fahrenheit()
{
}

double Fahrenheit::convertToInternal(double value) const
{
    return (value + 459.67) / 1.8;
}

double Fahrenheit::convertFromInternal(double value) const
{
    return (value * 1.8) - 459.67;
}

double Fahrenheit::convertFromInternal() const
{
    return Dimension::TempBase::convertFromInternal();
}

Fahrenheit &Fahrenheit::operator=(const Fahrenheit &other)
{
    Dimension::TempBase::operator=(other);
    return *this;
}

Fahrenheit &Fahrenheit::operator=(const Dimension::Temperature &temp)
{
    Dimension::TempBase::operator=(temp);
    return *this;
}

Fahrenheit Fahrenheit::operator-() const
{
    return Fahrenheit(-convertFromInternal());
}

Fahrenheit Fahrenheit::operator+(const Fahrenheit &other) const
{
    return Fahrenheit(convertFromInternal() + other.convertFromInternal());
}

Fahrenheit Fahrenheit::operator-(const Fahrenheit &other) const
{
    return Fahrenheit(convertFromInternal() - other.convertFromInternal());
}

Fahrenheit &Fahrenheit::operator+=(const Fahrenheit &other)
{
    convertToInternal(convertFromInternal() + other.convertFromInternal());
    return *this;
}

Fahrenheit &Fahrenheit::operator-=(const Fahrenheit &other)
{
    convertToInternal(convertFromInternal() - other.convertFromInternal());
    return *this;
}

Fahrenheit Fahrenheit::operator+(const Dimension::Temperature &other) const
{
    return *this + Fahrenheit(other);
}

Fahrenheit Fahrenheit::operator-(const Dimension::Temperature &other) const
{
    return *this - Fahrenheit(other);
}

Fahrenheit &Fahrenheit::operator+=(const Dimension::Temperature &other)
{
    return this->operator+=(Fahrenheit(other));
}

Fahrenheit &Fahrenheit::operator-=(const Dimension::Temperature &other)
{
    return this->operator-=(Fahrenheit(other));
}

Fahrenheit Fahrenheit::operator*(double value) const
{
    return Fahrenheit(value * convertFromInternal());
}

Fahrenheit Fahrenheit::operator/(double value) const
{
    return Fahrenheit(value / convertFromInternal());
}

Fahrenheit Fahrenheit::operator*(int value) const
{
    return Fahrenheit(value * convertFromInternal());
}

Fahrenheit Fahrenheit::operator/(int value) const
{
    return Fahrenheit(value / convertFromInternal());
}

QString Fahrenheit::unitString() const
{
    return "F";
}

SIREUNITS_EXPORT Fahrenheit SireUnits::operator*(double value, const Fahrenheit &temp)
{
    return temp * value;
}

SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> SireUnits::operator/(double value, const Fahrenheit &temp)
{
    return Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0>(value / temp.convertFromInternal());
}

SIREUNITS_EXPORT Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0> SireUnits::operator/(int value, const Fahrenheit &temp)
{
    return Dimension::PhysUnit<0, 0, 0, 0, -1, 0, 0>(value / temp.convertFromInternal());
}

SIREUNITS_EXPORT Fahrenheit SireUnits::operator*(int value, const Fahrenheit &temp)
{
    return temp * value;
}
