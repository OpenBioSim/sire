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

#include "chirality.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Chirality> r_chiral;
static const RegisterMetaType<AtomChiralities> r_atomchiralities;

QDataStream &operator<<(QDataStream &ds, const Chirality &c)
{
    writeHeader(ds, r_chiral, 1);

    ds << c.chiral_type;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Chirality &c)
{
    VersionID v = readHeader(ds, r_chiral);

    if (v == 1)
    {
        ds >> c.chiral_type;
    }
    else
        throw version_error(v, "1", r_chiral, CODELOC);

    return ds;
}

/** Constructor (default is an undefined chirality) */
Chirality::Chirality() : ConcreteProperty<Chirality, Property>(), chiral_type(0)
{
}

/** Construct from the the passed SDF number */
Chirality Chirality::fromSDF(int value)
{
    Chirality ret;

    return ret;
}

/** Construct from a string representation of a RDKit chirality */
Chirality Chirality::fromRDKit(const QString &value)
{
    Chirality ret;

    auto v = value.toUpper();

    if (v == "CHI_UNSPECIFIED")
    {
        ret.chiral_type = 0;
    }
    else if (v == "CHI_TETRAHEDRAL_CW")
    {
        ret.chiral_type = 1;
    }
    else if (v == "CHI_TETRAHEDRAL_CCW")
    {
        ret.chiral_type = 2;
    }
    else if (v == "CHI_OTHER")
    {
        ret.chiral_type = 3;
    }
    else if (v == "CHI_UNKNOWN")
    {
        ret.chiral_type = 4;
    }
    else if (v == "CHI_TETRAHEDRAL")
    {
        ret.chiral_type = 5;
    }
    else if (v == "CHI_ALLENE")
    {
        ret.chiral_type = 6;
    }
    else if (v == "CHI_SQUAREPLANAR")
    {
        ret.chiral_type = 7;
    }
    else if (v == "CHI_TRIGONALBIPYRAMIDAL")
    {
        ret.chiral_type = 8;
    }
    else if (v == "CHI_OCTAHEDRAL")
    {
        ret.chiral_type = 9;
    }
    else
    {
        ret.chiral_type = 3;
    }

    return ret;
}

/** Copy constructor */
Chirality::Chirality(const Chirality &other)
    : ConcreteProperty<Chirality, Property>(other), chiral_type(other.chiral_type)
{
}

/** Destructor */
Chirality::~Chirality()
{
}

/** Copy assignment operator */
Chirality &Chirality::operator=(const Chirality &other)
{
    chiral_type = other.chiral_type;
    return *this;
}

/** Comparison operator */
bool Chirality::operator==(const Chirality &other) const
{
    return chiral_type == other.chiral_type;
}

/** Comparison operator */
bool Chirality::operator!=(const Chirality &other) const
{
    return not this->operator==(other);
}

const char *Chirality::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Chirality>());
}

QString Chirality::toString() const
{
    switch (this->chiral_type)
    {
    case 0:
        return "unspecified";
    case 1:
        return "clockwise";
    case 2:
        return "counter_clockwise";
    case 3:
        return "other";
    case 4:
        return "unknown";
    case 5:
        return "tetrahedral";
    case 6:
        return "allene";
    case 7:
        return "square_planar";
    case 8:
        return "triganol_bipyramidal";
    case 9:
        return "octahedral";
    default:
        return "other";
    }
}

/** Return the SDF-format value for this bond. This returns
    '0' if the stereoscopy is undefined
 */
int Chirality::toSDF() const
{
    return 0;
}

/** Return a string representation of the RDKit stereo value */
QString Chirality::toRDKit() const
{
    switch (this->chiral_type)
    {
    case 0:
        return "CHI_UNSPECIFIED";
    case 1:
        return "CHI_TETRAHEDRAL_CW";
    case 2:
        return "CHI_TETRAHEDRAL_CCW";
    case 3:
        return "CHI_OTHER";
    case 4:
        return "CHI_UNKNOWN";
    case 5:
        return "CHI_TETRAHEDRAL";
    case 6:
        return "CHI_ALLENE";
    case 7:
        return "CHI_SQUAREPLANAR";
    case 8:
        return "CHI_TRIGONALBIPYRAMIDAL";
    case 9:
        return "CHI_OCTAHEDRAL";
    default:
        return "CHI_OTHER";
    }
}

/** Return a clockwise Chirality */
Chirality Chirality::clockwise()
{
    return Chirality::fromRDKit("CHI_TETRAHEDRAL_CW");
}

/** Return a counter-clockwise Chirality */
Chirality Chirality::counterClockwise()
{
    return Chirality::fromRDKit("CHI_TETRAHEDRAL_CCW");
}

/** Return an "other" Chirality */
Chirality Chirality::other()
{
    return Chirality::fromRDKit("CHI_OTHER");
}

/** Return an undefined Chirality */
Chirality Chirality::undefined()
{
    return Chirality::fromRDKit("CHI_UNSPECIFIED");
}

/** Return whether or not the stereoscopy is undefined */
bool Chirality::isUndefined() const
{
    return this->chiral_type == 0;
}

/** Return whether or not this is a clockwise chirality */
bool Chirality::isClockwise() const
{
    return this->chiral_type == 1;
}

/** Return whether or not this is a counter-clockwise chirality */
bool Chirality::isCounterClockwise() const
{
    return this->chiral_type == 2;
}

/** Return whether or not this is an "other" chirality */
bool Chirality::isOther() const
{
    return this->chiral_type == 3;
}
