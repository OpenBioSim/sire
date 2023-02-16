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

#include "bondorder.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<BondOrder> r_bondtype;

QDataStream &operator<<(QDataStream &ds, const BondOrder &b)
{
    writeHeader(ds, r_bondtype, 1);

    ds << b.bond_type;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BondOrder &b)
{
    VersionID v = readHeader(ds, r_bondtype);

    if (v == 1)
    {
        ds >> b.bond_type;
    }
    else
        throw version_error(v, "1", r_bondtype, CODELOC);

    return ds;
}

/** Constructor (default is an undefined bond) */
BondOrder::BondOrder() : ConcreteProperty<BondOrder, Property>(), bond_type(0)
{
}

/** Construct from the passed string */
BondOrder::BondOrder(const QString &str) : ConcreteProperty<BondOrder, Property>()
{
    auto s = str.trimmed().toLower();

    if (s == "single")
        this->bond_type = 1;
    else if (s == "double")
        this->bond_type = 2;
    else if (s == "triple")
        this->bond_type = 3;
    else if (s == "aromatic")
        this->bond_type = 4;
    else if (s == "undefined")
        this->bond_type = 0;
    else
        throw SireError::invalid_arg(QObject::tr("Cannot interpret bond type '%1'. Should be one of "
                                                 "'single', 'double', 'triple', 'aromatic' or 'undefined'.")
                                         .arg(str),
                                     CODELOC);
}

/** Construct from the the passed SDF number */
BondOrder BondOrder::fromSDF(int value)
{
    BondOrder ret;

    if (value < 0 or value > 4)
        throw SireError::invalid_arg(QObject::tr("Invalid bond type '%1'. Should be an integer between "
                                                 "0 and 4.")
                                         .arg(value),
                                     CODELOC);

    ret.bond_type = value;
    return ret;
}

/** Construct from a string representation of the RDKit type */
BondOrder BondOrder::fromRDKit(const QString &value)
{
    BondOrder ret;
    return ret;
}

/** Copy constructor */
BondOrder::BondOrder(const BondOrder &other) : ConcreteProperty<BondOrder, Property>(other), bond_type(other.bond_type)
{
}

/** Destructor */
BondOrder::~BondOrder()
{
}

/** Copy assignment operator */
BondOrder &BondOrder::operator=(const BondOrder &other)
{
    bond_type = other.bond_type;
    return *this;
}

/** Comparison operator */
bool BondOrder::operator==(const BondOrder &other) const
{
    return bond_type == other.bond_type;
}

/** Comparison operator */
bool BondOrder::operator!=(const BondOrder &other) const
{
    return not this->operator==(other);
}

const char *BondOrder::typeName()
{
    return QMetaType::typeName(qMetaTypeId<BondOrder>());
}

QString BondOrder::toString() const
{
    switch (this->bond_type)
    {
    case 0:
        return "undefined";
    case 1:
        return "single";
    case 2:
        return "double";
    case 3:
        return "triple";
    case 4:
        return "aromatic";
    default:
        throw SireError::program_bug(QObject::tr("Should not get here: %1").arg(this->bond_type), CODELOC);
    }
}

/** Return the bond type (uses SDF values, e.g. 0 is undefined,
    1 is single, 2 is double, 3 is triple and 4 is aromatic)
*/
int BondOrder::value() const
{
    return this->bond_type;
}

/** Return the bond order as a double precision number. This matches
 *  the value that would be returned by RDKit
 */
double BondOrder::valueAsDouble() const
{
    return 0.0;
}

/** Return the SDF-format value for this bond */
int BondOrder::toSDF() const
{
    return this->bond_type;
}

/** Return as a string representation of a RDKit bond order */
QString BondOrder::toRDKit() const
{
    return QString();
}

/** Return a single bond */
BondOrder BondOrder::singleBond()
{
    return BondOrder::fromSDF(1);
}

/** Return a double bond */
BondOrder BondOrder::doubleBond()
{
    return BondOrder::fromSDF(2);
}

/** Return a triple bond */
BondOrder BondOrder::tripleBond()
{
    return BondOrder::fromSDF(3);
}

/** Return an aromatic bond */
BondOrder BondOrder::aromaticBond()
{
    return BondOrder::fromSDF(4);
}

/** Return an undefined bond */
BondOrder BondOrder::undefinedBond()
{
    return BondOrder::fromSDF(0);
}

/** Return whether or not the bond type is defined */
bool BondOrder::isDefined() const
{
    return this->bond_type != 0;
}

/** Return whether or not this is a single bond */
bool BondOrder::isSingle() const
{
    return this->bond_type == 1;
}

/** Return whether or not this is a double bond */
bool BondOrder::isDouble() const
{
    return this->bond_type == 2;
}

/** Return whether or not this is a triple bond */
bool BondOrder::isTriple() const
{
    return this->bond_type == 3;
}

/** Return whether or not this is an aromatic bond */
bool BondOrder::isAromatic() const
{
    return this->bond_type == 4;
}
