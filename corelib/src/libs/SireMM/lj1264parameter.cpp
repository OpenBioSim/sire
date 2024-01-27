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

#include "lj1264parameter.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;

static const RegisterMetaType<LJ1264Parameter> r_ljparam(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const LJ1264Parameter &ljparam)
{
    writeHeader(ds, r_ljparam, 1) << ljparam.a << ljparam.b << ljparam.c;

    return ds;
}


/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, LJ1264Parameter &ljparam)
{
    VersionID v = readHeader(ds, r_ljparam);

    if (v == 1)
    {
        ds >> ljparam.a >> ljparam.b >> ljparam.c;
    }
    else
        throw version_error(v, "1", r_ljparam, CODELOC);

    return ds;
}

/** Construct a dummy LJ parameter */
LJ1264Parameter::LJ1264Parameter() : a(0), b(0), c(0)
{
}

/** Construct from an LJPair */
LJ1264Parameter::LJ1264Parameter(const LJPair &ljpair)
    : a(ljpair.A()), b(ljpair.B()), c(0)
{
}

/** Construct from a LJParameter */
LJ1264Parameter::LJ1264Parameter(const LJParameter &ljparam)
    : a(ljparam.A()), b(ljparam.B()), c(0)
{
}

/** Construct from specified parameters - note that you must have
 *  converted these into the correct internal units (kcal mol-1, Angstrom)
 */
LJ1264Parameter::LJ1264Parameter(double _a, double _b, double _c)
    : a(_a), b(_b), c(_c)
{
}

/** Return the correct physical units for the A parameter */
GeneralUnit LJ1264Parameter::AUnit()
{
    return GeneralUnit(kcal_per_mol * angstrom3 * angstrom3 * angstrom3 * angstrom3);
}

/** Return the correct physical units for the B parameter */
GeneralUnit LJ1264Parameter::BUnit()
{
    return GeneralUnit(kcal_per_mol * angstrom3 * angstrom3);
}

/** Return the correct physical units for the C parameter */
GeneralUnit LJ1264Parameter::CUnit()
{
    return GeneralUnit(kcal_per_mol * angstrom3 * angstrom);
}

/** Construct from the passed (dimensioned) values */
LJ1264Parameter::LJ1264Parameter(const GeneralUnit &_a,
                                 const GeneralUnit &_b,
                                 const GeneralUnit &_c)
    : a(_a.value()), b(_b.value()), c(_c.value())
{
    if (not(_a.hasSameUnits(AUnit()) and _b.hasSameUnits(BUnit()) and _c.hasSameUnits(CUnit())))
        throw SireError::incompatible_error(QObject::tr(
                                                "LJ1264Parameter::LJ1264Parameter: incompatible units for LJ parameters: %1, %2, %3."
                                                "Units should be %4, %5, %6.")
                                                .arg(_a.toString())
                                                .arg(_b.toString())
                                                .arg(_c.toString())
                                                .arg(AUnit().unitString())
                                                .arg(BUnit().unitString())
                                                .arg(CUnit().unitString()));
}

/** Copy constructor */
LJ1264Parameter::LJ1264Parameter(const LJ1264Parameter &other)
    : a(other.a), b(other.b), c(other.c)
{
}

/** Destructor */
LJ1264Parameter::~LJ1264Parameter()
{
}

/** Copy assignment operator */
LJ1264Parameter &LJ1264Parameter::operator=(const LJ1264Parameter &other)
{
    if (this != &other)
    {
        a = other.a;
        b = other.b;
        c = other.c;
    }

    return *this;
}

/** Return a dummy LJ1264Parameter */
LJ1264Parameter LJ1264Parameter::dummy()
{
    return LJ1264Parameter();
}

/** Return a string representation of the LJ parameter */
QString LJ1264Parameter::toString() const
{
    if (c == 0)
        return this->toLJParameter().toString();
    else
        return QString("LJ( A = %1, B = %2, C = %3 )")
            .arg((this->A() * AUnit()).toString())
            .arg((this->B() * BUnit()).toString())
            .arg((this->C() * CUnit()).toString());
}

const char *LJ1264Parameter::typeName()
{
    return QMetaType::typeName(qMetaTypeId<LJ1264Parameter>());
}

const char *LJ1264Parameter::what() const
{
    return LJ1264Parameter::typeName();
}

/** Return a clone of this LJ1264Parameter */
LJ1264Parameter *LJ1264Parameter::clone() const
{
    return new LJ1264Parameter(*this);
}

/** Return whether or not this has a C parameter. If so,
 *  then it is a 12-6-4 parameter. If not, then it is a standard
 *  12-6 parameter.
 */
bool LJ1264Parameter::hasC() const
{
    return c != 0;
}

/** Return whether or not this is a standard 12-6 LJ parameter */
bool LJ1264Parameter::isLJParameter() const
{
    return c == 0;
}

/** Assert that this is a standard 12-6 LJ parameter (C is 0) */
void LJ1264Parameter::assertIsLJParameter() const
{
    if (not isLJParameter())
        throw SireError::incompatible_error(QObject::tr(
            "LJ1264Parameter::assertIsLJParameter: this is not a standard 12-6 LJ parameter."));
}

/** Return this converted to a standard 12-6 LJ parameter */
LJParameter LJ1264Parameter::toLJParameter() const
{
    this->assertIsLJParameter();

    if (a == 0 or b == 0)
        return LJParameter::dummy();
    else
        return LJParameter::fromAAndB(a, b);
}

/** Return this converted to a standard 12-6 LJ parameter */
LJPair LJ1264Parameter::toLJPair() const
{
    this->assertIsLJParameter();

    if (a == 0 or b == 0)
        return LJPair::dummy();
    else
        return LJPair::fromAAndB(a, b);
}
