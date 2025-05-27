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

#include "rmsdrestraints.h"

#include "SireUnits/units.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMM;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

////////////
//////////// Implementation of RMSDRestraint
////////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const RMSDRestraint &rmsdrest)
{
    writeHeader(ds, r_rmsdrest, 1);

    SharedDataStream sds(ds);

    sds << rmsdrest.atms << rmsdrest.ref_pos << rmsdrest._k << rmsdrest._r0
        << static_cast<const Property &>(rmsdrest);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, RMSDRestraint &rmsdrest)
{
    VersionID v = readHeader(ds, r_rmsdrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> rmsdrest.atms >> rmsdrest.ref_pos >> rmsdrest._k >> rmsdrest._r0 >> static_cast<Property &>(rmsdrest);
    }
    else
        throw version_error(v, "1", r_rmsdrest, CODELOC);

    return ds;
}

/** Null constructor */
RMSDRestraint::RMSDRestraint()
    : ConcreteProperty<RMSDRestraint, Property>(),
      _k(0), _r0(0)
{
}


/** Construct a restraint that acts on the group of provided atoms,
    restraining these wrt the desired reference positions */
RMSDRestraint::RMSDRestraint(const QList<qint64> &atoms, 
                             const QList<SireMaths::Vector> &ref_positions,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0)
    : ConcreteProperty<RMSDRestraint, Property>(),
      ref_pos(ref_positions.toVector()), _k(k), _r0(r0)
{
    // Make sure that we have a set of distinct atoms
    // QSet<qint64> distinct;
    // distinct.reserve(atoms.size());

    // for (const auto &atom : atoms)
    // {
    //     if (atom >= 0)
    //         distinct.insert(atom);
    // }

    // Convert atom indices to QVector
    atms = atoms.toVector();
}

/** Copy constructor */
RMSDRestraint::RMSDRestraint(const RMSDRestraint &other)
    : ConcreteProperty<RMSDRestraint, Property>(other),
      atms(other.atms), ref_pos(other.ref_pos), _k(other._k), _r0(other._r0)
{
}

/* Destructor */
RMSDRestraint::~RMSDRestraint()
{
}

RMSDRestraint &RMSDRestraint::operator=(const RMSDRestraint &other)
{
    if (this != &other)
    {
        atms = other.atms;
        ref_pos = other.ref_pos;
        _k = other._k;
        _r0 = other._r0;
    }

    return *this;
}

bool RMSDRestraint::operator==(const RMSDRestraint &other) const
{
    return atms == other.atms and ref_pos == other.ref_pos and
           _k == other._k and _r0 == other._r0;
}

bool RMSDRestraint::operator!=(const RMSDRestraint &other) const
{
    return not operator==(other);
}

RMSDRestraints RMSDRestraint::operator+(const RMSDRestraint &other) const
{
    return RMSDRestraints(*this) + other;
}

RMSDRestraints RMSDRestraint::operator+(const RMSDRestraints &other) const
{
    return RMSDRestraints(*this) + other;
}

const char *RMSDRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<RMSDRestraint>());
}

const char *RMSDRestraint::what() const
{
    return RMSDRestraint::typeName();
}

RMSDRestraint *RMSDRestraint::clone() const
{
    return new RMSDRestraint(*this);
}

bool RMSDRestraint::isNull() const
{
    return atms.isEmpty();
}

QString RMSDRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("RMSDRestraint::null");

    else
    {
        return QString("RMSDRestraint( %1 => %2, k=%3 : r0=%4 )")
            .arg(this->atom())
            .arg(ref_pos.toString())
            .arg(_k.toString())
            .arg(_r0.toString());
    }
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicBondConstant RMSDRestraint::k() const
{
    return this->_k;
}

/** Return the width of the flat-bottom well. This is zero for a
 *  pure harmonic restraint
 */
SireUnits::Dimension::Length RMSDRestraint::r0() const
{
    return this->_r0;
}

///////
/////// Implementation of RMSDRestraints
///////

static const RegisterMetaType<RMSDRestraints> r_rmsdrests;

QDataStream &operator<<(QDataStream &ds, const RMSDRestraints &rmsdrests)
{
    writeHeader(ds, r_rmsdrests, 1);

    SharedDataStream sds(ds);

    sds << rmsdrests.r
        << static_cast<const Restraints &>(rmsdrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, RMSDRestraints &rmsdrests)
{
    VersionID v = readHeader(ds, r_rmsdrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> rmsdrests.r >> static_cast<Restraints &>(rmsdrests);
    }
    else
        throw version_error(v, "1", r_rmsdrests, CODELOC);

    return ds;
}

/** Null constructor */
RMSDRestraints::RMSDRestraints()
    : ConcreteProperty<RMSDRestraints, Restraints>()
{
}

RMSDRestraints::RMSDRestraints(const QString &name)
    : ConcreteProperty<RMSDRestraints, Restraints>(name)
{
}

RMSDRestraints::RMSDRestraints(const RMSDRestraint &restraint)
    : ConcreteProperty<RMSDRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

RMSDRestraints::RMSDRestraints(const QList<RMSDRestraint> &restraints)
    : ConcreteProperty<RMSDRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

RMSDRestraints::RMSDRestraints(const QString &name,
                                           const RMSDRestraint &restraint)
    : ConcreteProperty<RMSDRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

RMSDRestraints::RMSDRestraints(const QString &name,
                                           const QList<RMSDRestraint> &restraints)
    : ConcreteProperty<RMSDRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

RMSDRestraints::RMSDRestraints(const RMSDRestraints &other)
    : ConcreteProperty<RMSDRestraints, Restraints>(other), r(other.r)
{
}

RMSDRestraints::~RMSDRestraints()
{
}

RMSDRestraints &RMSDRestraints::operator=(const RMSDRestraints &other)
{
    r = other.r;
    Restraints::operator=(other);
    return *this;
}

bool RMSDRestraints::operator==(const RMSDRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
}

bool RMSDRestraints::operator!=(const RMSDRestraints &other) const
{
    return not operator==(other);
}

const char *RMSDRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<RMSDRestraints>());
}

const char *RMSDRestraints::what() const
{
    return RMSDRestraints::typeName();
}

QString RMSDRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("RMSDRestraints::null");

    QStringList parts;

    const auto n = this->count();

    if (n <= 10)
    {
        for (int i = 0; i < n; ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }
    else
    {
        for (int i = 0; i < 5; ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }

        parts.append("...");

        for (int i = n - 5; i < n; ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }

    return QObject::tr("RMSDRestraints( name=%1, size=%2\n%3\n)").arg(this->name()).arg(n).arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool RMSDRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool RMSDRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int RMSDRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int RMSDRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int RMSDRestraints::size() const
{
    return this->nRestraints();
}

/** Return the ith restraint */
const RMSDRestraint &RMSDRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const RMSDRestraint &RMSDRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<RMSDRestraint> RMSDRestraints::restraints() const
{
    return this->r;
}

/** Add a restraint onto the list */
void RMSDRestraints::add(const RMSDRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void RMSDRestraints::add(const RMSDRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
RMSDRestraints &RMSDRestraints::operator+=(const RMSDRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
RMSDRestraints RMSDRestraints::operator+(const RMSDRestraint &restraint) const
{
    RMSDRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
RMSDRestraints &RMSDRestraints::operator+=(const RMSDRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
RMSDRestraints RMSDRestraints::operator+(const RMSDRestraints &restraints) const
{
    RMSDRestraints ret(*this);
    ret += restraints;
    return *this;
}
