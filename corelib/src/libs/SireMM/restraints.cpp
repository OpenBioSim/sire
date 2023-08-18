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
  *  at http://sire.openbiosim.org
  *
\*********************************************/

#include "restraints.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Restraints> r_restraints(MAGIC_ONLY,
                                                       Restraints::typeName());

QDataStream &operator<<(QDataStream &ds, const Restraints &restraints)
{
    writeHeader(ds, r_restraints, 1);

    ds << restraints.nme << static_cast<const Property &>(restraints);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Restraints &restraints)
{
    VersionID v = readHeader(ds, r_restraints);

    if (v == 1)
    {
        ds >> restraints.nme >> static_cast<Property &>(restraints);
    }
    else
        throw version_error(v, "1", r_restraints, CODELOC);

    return ds;
}

Restraints::Restraints() : Property(), nme("restraint")
{
}

Restraints::Restraints(const QString &name)
    : Property()
{
    this->setName(name);
}

Restraints::Restraints(const Restraints &other)
    : Property(other), nme(other.nme)
{
}

Restraints::~Restraints()
{
}

/** Return whether or not this is a Restraints object (or derived
    from this object) */
bool Restraints::isRestraints() const
{
    return true;
}

/** Return the name given to this group of restraints */
QString Restraints::name() const
{
    return this->nme;
}

/** Set the name for this group of restraints */
void Restraints::setName(const QString &name)
{
    this->nme = name.simplified();

    if (this->nme.isEmpty())
        this->nme = "restraints";
}

Restraints &Restraints::operator=(const Restraints &other)
{
    nme = other.nme;
    Property::operator=(other);
}

bool Restraints::operator==(const Restraints &other) const
{
    return nme == other.nme;
}

bool Restraints::operator!=(const Restraints &other) const
{
    return not this->operator==(other);
}
