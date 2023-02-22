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
        this->operator=(BondOrder::fromRDKit(str.toUpper()));
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

    const static std::map<QString, qint32> types = {
        {"UNSPECIFIED", 0},
        {"SINGLE", 1},
        {"DOUBLE", 2},
        {"TRIPLE", 3},
        {"QUADRUPLE", 104},
        {"QUINTUPLE", 105},
        {"HEXTUPLE", 106},
        {"ONEANDAHALF", 107},
        {"TWOANDAHALF", 108},
        {"THREEANDAHALF", 109},
        {"FOURANDAHALF", 110},
        {"FIVEANDAHALF", 111},
        {"AROMATIC", 4},
        {"IONIC", 112},
        {"HYDROGEN", 113},
        {"THREECENTER", 114},
        {"DATIVEONE", 115},
        {"DATIVE", 116},
        {"DATIVEL", 117},
        {"DATIVER", 118},
        {"OTHER", 119},
        {"ZERO", 100}};

    auto it = types.find(value);

    if (it == types.end())
        ret.bond_type = 0;
    else
        ret.bond_type = it->second;

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
    case 104:
        return "quadruple";
    case 105:
        return "quintuple";
    case 106:
        return "hextuple";
    case 107:
        return "oneandahalf";
    case 108:
        return "twoandahalf";
    case 109:
        return "threeandahalf";
    case 110:
        return "fourandahalf";
    case 111:
        return "fiveandahalf";
    case 112:
        return "ionic";
    case 113:
        return "hydrogen";
    case 114:
        return "threecenter";
    case 115:
        return "dativeone";
    case 116:
        return "dative";
    case 117:
        return "dativel";
    case 118:
        return "dativer";
    case 119:
        return "other";
    case 100:
        return "zero";
    default:
        throw "other";
    }
}

/** Return the bond type (0 to number of bonds)
 */
int BondOrder::value() const
{
    switch (this->bond_type)
    {
    case 0:
        return 0;
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
        return 3;
    case 4:
        return 1;
    case 104:
        return 4;
    case 105:
        return 5;
    case 106:
        return 6;
    case 107:
        return 1;
    case 108:
        return 2;
    case 109:
        return 3;
    case 110:
        return 4;
    case 111:
        return 5;
    case 112:
        return 0;
    case 113:
        return 0;
    case 114:
        return 0;
    case 115:
        return 1;
    case 116:
        return 1;
    case 117:
        return 1;
    case 118:
        return 1;
    case 119:
        return 1;
    case 100:
        return 0;
    default:
        return 1;
    }
}

/** Return the bond order as a double precision number. This matches
 *  the value that would be returned by RDKit
 */
double BondOrder::valueAsDouble() const
{
    switch (this->bond_type)
    {
    case 0:
        return 0.0;
    case 1:
        return 1.0;
    case 2:
        return 2.0;
    case 3:
        return 3.0;
    case 4:
        return 1.5;
    case 104:
        return 4.0;
    case 105:
        return 5.0;
    case 106:
        return 6.0;
    case 107:
        return 1.5;
    case 108:
        return 2.5;
    case 109:
        return 3.5;
    case 110:
        return 4.5;
    case 111:
        return 5.5;
    case 112:
        return 0.5;
    case 113:
        return 0.5;
    case 114:
        return 0.6;
    case 115:
        return 1.0;
    case 116:
        return 1.0;
    case 117:
        return 1.0;
    case 118:
        return 1.0;
    case 119:
        return 1.0;
    case 100:
        return 0.0;
    default:
        return 1.0;
    }
}

/** Return the SDF-format value for this bond */
int BondOrder::toSDF() const
{
    switch (this->bond_type)
    {
    case 0:
        return 0;
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
        return 3;
    case 4:
        return 4;
    case 104:
        return 3;
    case 105:
        return 3;
    case 106:
        return 3;
    case 107:
        return 4;
    case 108:
        return 4;
    case 109:
        return 4;
    case 110:
        return 4;
    case 111:
        return 4;
    case 112:
        return 0;
    case 113:
        return 0;
    case 114:
        return 0;
    case 115:
        return 1;
    case 116:
        return 1;
    case 117:
        return 1;
    case 118:
        return 1;
    case 119:
        return 0;
    case 100:
        return 0;
    default:
        return 1;
    }
}

/** Return as a string representation of a RDKit bond order */
QString BondOrder::toRDKit() const
{
    switch (this->bond_type)
    {
    case 0:
        return "UNSPECIFIED";
    case 1:
        return "SINGLE";
    case 2:
        return "DOUBLE";
    case 3:
        return "TRIPLE";
    case 4:
        return "AROMATIC";
    case 104:
        return "QUADRUPLE";
    case 105:
        return "QUINTUPLE";
    case 106:
        return "HEXTUPLE";
    case 107:
        return "ONEANDAHALF";
    case 108:
        return "TWOANDAHALF";
    case 109:
        return "THREEANDAHALF";
    case 110:
        return "FOURANDAHALF";
    case 111:
        return "FIVEANDAHALF";
    case 112:
        return "IONIC";
    case 113:
        return "HYDROGEN";
    case 114:
        return "THREECENTER";
    case 115:
        return "DATIVEONE";
    case 116:
        return "DATIVE";
    case 117:
        return "DATIVEL";
    case 118:
        return "DATIVER";
    case 119:
        return "OTHER";
    case 100:
        return "ZERO";
    default:
        return "OTHER";
    }
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
