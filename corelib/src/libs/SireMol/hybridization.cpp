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

#include "hybridization.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Hybridization> r_hybrid;
static const RegisterMetaType<AtomHybridizations> r_atomhybridizations;

QDataStream &operator<<(QDataStream &ds, const Hybridization &h)
{
    writeHeader(ds, r_hybrid, 1);

    ds << h.hybrid_type;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Hybridization &h)
{
    VersionID v = readHeader(ds, r_hybrid);

    if (v == 1)
    {
        ds >> h.hybrid_type;
    }
    else
        throw version_error(v, "1", r_hybrid, CODELOC);

    return ds;
}

/** Constructor (default is an undefined chirality) */
Hybridization::Hybridization() : ConcreteProperty<Hybridization, Property>(), hybrid_type(0)
{
}

/** Construct from the the passed SDF number */
Hybridization Hybridization::fromSDF(int value)
{
    Hybridization ret;

    return ret;
}

/** Construct from a string representation of a RDKit hybridization */
Hybridization Hybridization::fromRDKit(const QString &value)
{
    Hybridization ret;

    auto v = value.toUpper();

    if (v == "UNSPECIFIED")
    {
        ret.hybrid_type = 0;
    }
    else if (v == "S")
    {
        ret.hybrid_type = 1;
    }
    else if (v == "SP")
    {
        ret.hybrid_type = 2;
    }
    else if (v == "SP2")
    {
        ret.hybrid_type = 3;
    }
    else if (v == "SP3")
    {
        ret.hybrid_type = 4;
    }
    else if (v == "SP2D")
    {
        ret.hybrid_type = 5;
    }
    else if (v == "SP3D")
    {
        ret.hybrid_type = 6;
    }
    else if (v == "SP3D2")
    {
        ret.hybrid_type = 7;
    }
    else if (v == "OTHER")
    {
        ret.hybrid_type = 8;
    }
    else if (v == "UNKNOWN")
    {
        ret.hybrid_type = -1;
    }

    return ret;
}

/** Copy constructor */
Hybridization::Hybridization(const Hybridization &other)
    : ConcreteProperty<Hybridization, Property>(other), hybrid_type(other.hybrid_type)
{
}

/** Destructor */
Hybridization::~Hybridization()
{
}

/** Copy assignment operator */
Hybridization &Hybridization::operator=(const Hybridization &other)
{
    hybrid_type = other.hybrid_type;
    return *this;
}

/** Comparison operator */
bool Hybridization::operator==(const Hybridization &other) const
{
    return hybrid_type == other.hybrid_type;
}

/** Comparison operator */
bool Hybridization::operator!=(const Hybridization &other) const
{
    return not this->operator==(other);
}

const char *Hybridization::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Hybridization>());
}

QString Hybridization::toString() const
{
    return this->toRDKit().toLower();
}

/** Return the SDF-format value for this bond. This returns
    '0' if the stereoscopy is undefined
 */
int Hybridization::toSDF() const
{
    return 0;
}

/** Return a string representation of the RDKit stereo value */
QString Hybridization::toRDKit() const
{
    switch (this->hybrid_type)
    {
    case -1:
        return "UNKNOWN";
    case 0:
        return "UNSPECIFIED";
    case 1:
        return "S";
    case 2:
        return "SP";
    case 3:
        return "SP2";
    case 4:
        return "SP3";
    case 5:
        return "SP2D";
    case 6:
        return "SP3D";
    case 7:
        return "SP3D2";
    case 8:
        return "OTHER";
    default:
        return "UNSPECIFIED";
    }
}

Hybridization Hybridization::s()
{
    return Hybridization::fromRDKit("S");
}

Hybridization Hybridization::sp()
{
    return Hybridization::fromRDKit("SP");
}

Hybridization Hybridization::sp2()
{
    return Hybridization::fromRDKit("SP2");
}

Hybridization Hybridization::sp3()
{
    return Hybridization::fromRDKit("SP3");
}

Hybridization Hybridization::other()
{
    return Hybridization::fromRDKit("OTHER");
}

Hybridization Hybridization::unspecified()
{
    return Hybridization::fromRDKit("UNSPECIFIED");
}

Hybridization Hybridization::unknown()
{
    return Hybridization::fromRDKit("UNKNOWN");
}

bool Hybridization::is_s() const
{
    return this->hybrid_type == 1;
}

bool Hybridization::is_sp() const
{
    return this->hybrid_type == 2;
}

bool Hybridization::is_sp2() const
{
    return this->hybrid_type == 3;
}

bool Hybridization::is_sp3() const
{
    return this->hybrid_type == 4;
}

bool Hybridization::isOther() const
{
    return this->hybrid_type == 8;
}

bool Hybridization::isUnspecified() const
{
    return this->hybrid_type == 0;
}

bool Hybridization::isUnknown() const
{
    return this->hybrid_type == -1;
}
