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

#include "stereochemistry.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Stereochemistry> r_stereo;

QDataStream &operator<<(QDataStream &ds, const Stereochemistry &s)
{
    writeHeader(ds, r_stereo, 1);

    ds << s.stereo_type;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Stereochemistry &s)
{
    VersionID v = readHeader(ds, r_stereo);

    if (v == 1)
    {
        ds >> s.stereo_type;
    }
    else
        throw version_error(v, "1", r_stereo, CODELOC);

    return ds;
}

/** Constructor (default is an undefined stereoscopy) */
Stereochemistry::Stereochemistry() : ConcreteProperty<Stereochemistry, Property>(), stereo_type(-1)
{
}

/** Construct from the passed string */
Stereochemistry::Stereochemistry(const QString &str) : ConcreteProperty<Stereochemistry, Property>()
{
    auto s = str.trimmed().toLower();

    if (s == "up")
        this->stereo_type = 1;
    else if (s == "cis or trans")
        this->stereo_type = 3;
    else if (s == "down")
        this->stereo_type = 6;
    else if (s == "not stereo")
        this->stereo_type = 0;
    else if (s == "undefined")
        this->stereo_type = -1;
    else
        this->operator=(Stereochemistry::fromRDKit(str));
}

/** Construct from the the passed SDF number */
Stereochemistry Stereochemistry::fromSDF(int value)
{
    Stereochemistry ret;

    if (value == 0 or value == 1 or value == -1 or value == 3 or value == 6)
    {
        ret.stereo_type = value;
    }
    else
    {
        throw SireError::invalid_arg(QObject::tr("Invalid stereo type '%1'. Should be an integer in "
                                                 "[-1, 0, 1, 3, 6]")
                                         .arg(value),
                                     CODELOC);
    }

    return ret;
}

/** Construct from a string representation of a RDKit stereochemistry */
Stereochemistry Stereochemistry::fromRDKit(const QString &value)
{
    Stereochemistry ret;

    auto v = value.toLower();

    if (v == "stereonone")
    {
        // not stereo
        ret.stereo_type = 0;
    }
    else if (v == "stereoany")
    {
        // undefined
        ret.stereo_type = -1;
    }
    else if (v == "stereoz")
    {
        // RDKit
        ret.stereo_type = 1000;
    }
    else if (v == "stereoe")
    {
        // RDKit
        ret.stereo_type = 1001;
    }
    else if (v == "stereocis")
    {
        // RDKit
        ret.stereo_type = 1002;
    }
    else if (v == "stereotrans")
    {
        // RDKit
        ret.stereo_type = 1003;
    }
    else
    {
        // undefined
        ret.stereo_type = -1;
    }

    return ret;
}

/** Copy constructor */
Stereochemistry::Stereochemistry(const Stereochemistry &other)
    : ConcreteProperty<Stereochemistry, Property>(other), stereo_type(other.stereo_type)
{
}

/** Destructor */
Stereochemistry::~Stereochemistry()
{
}

/** Copy assignment operator */
Stereochemistry &Stereochemistry::operator=(const Stereochemistry &other)
{
    stereo_type = other.stereo_type;
    return *this;
}

/** Comparison operator */
bool Stereochemistry::operator==(const Stereochemistry &other) const
{
    return stereo_type == other.stereo_type;
}

/** Comparison operator */
bool Stereochemistry::operator!=(const Stereochemistry &other) const
{
    return not this->operator==(other);
}

const char *Stereochemistry::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Stereochemistry>());
}

QString Stereochemistry::toString() const
{
    switch (this->stereo_type)
    {
    case 0:
        return "not stereo";
    case 1:
        return "up";
    case 3:
        return "cis or trans";
    case 6:
        return "down";
    case -1:
        return "undefined";
    default:
        return this->toRDKit();
    }
}

/** Return the stereo type (uses SDF values, e.g. 0 is not stereo,
    1 is up, 6 is down. We have added -1 to mean undefined)
*/
int Stereochemistry::value() const
{
    return this->stereo_type;
}

/** Return the SDF-format value for this bond. This returns
    '0' if the stereoscopy is undefined
 */
int Stereochemistry::toSDF() const
{
    switch (this->stereo_type)
    {
    case 1:
        return 1;
    case 3:
        return 3;
    case 6:
        return 6;
    default:
        // either undefined or an unsupported RDKit stereo value
        return 0;
    }
}

/** Return a string representation of the RDKit stereo value */
QString Stereochemistry::toRDKit() const
{
    switch (this->stereo_type)
    {
    case -1:
        return "STEREOANY";
    case 0:
        return "STEREONONE";
    case 1000:
        return "STEREOZ";
    case 1001:
        return "STEREOE";
    case 1002:
        return "STEREOCIS";
    case 1003:
        return "STEREOTRANS";
    default:
        return "STEREOANY";
    }
}

/** Return an up Stereochemistry */
Stereochemistry Stereochemistry::up()
{
    return Stereochemistry::fromSDF(1);
}

/** Return a cis or trans Stereochemistry */
Stereochemistry Stereochemistry::cisOrTrans()
{
    return Stereochemistry::fromSDF(3);
}

/** Return a down Stereochemistry */
Stereochemistry Stereochemistry::down()
{
    return Stereochemistry::fromSDF(6);
}

/** Return a "not stereo" Stereochemistry */
Stereochemistry Stereochemistry::notStereo()
{
    return Stereochemistry::fromSDF(0);
}

/** Return an undefined Stereochemistry */
Stereochemistry Stereochemistry::undefined()
{
    return Stereochemistry::fromSDF(-1);
}

/** Return whether or not the stereoscopy is defined */
bool Stereochemistry::isDefined() const
{
    return this->stereo_type != -1;
}

/** Return whether or not this is an up bond */
bool Stereochemistry::isUp() const
{
    return this->stereo_type == 1;
}

/** Return whether or not this is a cis or trans bond */
bool Stereochemistry::isCisOrTrans() const
{
    return this->stereo_type == 3;
}

/** Return whether or not this is a down bond */
bool Stereochemistry::isDown() const
{
    return this->stereo_type == 6;
}

/** Return whether or not this is a "not stereo" bond */
bool Stereochemistry::isNotStereo() const
{
    return this->stereo_type == 0;
}
