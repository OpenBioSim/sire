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

#include "positionalrestraints.h"

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

///////
/////// Implementation of PositionalRestraint
///////

static const RegisterMetaType<PositionalRestraint> r_posrest;

QDataStream &operator<<(QDataStream &ds, const PositionalRestraint &posrest)
{
    writeHeader(ds, r_posrest, 1);

    SharedDataStream sds(ds);

    sds << posrest.atms << posrest.pos << posrest._k << posrest._r0
        << static_cast<const Property &>(posrest);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PositionalRestraint &posrest)
{
    VersionID v = readHeader(ds, r_posrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> posrest.atms >> posrest.pos >> posrest._k >> posrest._r0 >> static_cast<Property &>(posrest);
    }
    else
        throw version_error(v, "1", r_posrest, CODELOC);

    return ds;
}

/** Null constructor */
PositionalRestraint::PositionalRestraint()
    : ConcreteProperty<PositionalRestraint, Property>(),
      _k(0), _r0(0)
{
}

/** Construct to restrain the atom at index 'atom' to the specified position
 *  using the specified force constant and flat-bottom well-width
 */
PositionalRestraint::PositionalRestraint(qint64 atom, const SireMaths::Vector &position,
                                         const SireUnits::Dimension::HarmonicBondConstant &k,
                                         const SireUnits::Dimension::Length &r0)
    : ConcreteProperty<PositionalRestraint, Property>(),
      pos(position), _k(k), _r0(r0)
{
    atms = QVector<qint64>(1, atom);
    atms.squeeze();
}

/** Construct to restrain the centroid of the atoms whose indicies are
 *  in 'atoms' to the specified position using the specified force constant
 *  and flat-bottom well width
 */
PositionalRestraint::PositionalRestraint(const QList<qint64> &atoms,
                                         const SireMaths::Vector &position,
                                         const SireUnits::Dimension::HarmonicBondConstant &k,
                                         const SireUnits::Dimension::Length &r0)
    : ConcreteProperty<PositionalRestraint, Property>(),
      pos(position), _k(k), _r0(r0)
{
    if (atoms.isEmpty())
        return;

    // remove duplicates
    atms.reserve(atoms.count());

    auto sorted = atoms;
    std::sort(sorted.begin(), sorted.end());

    atms.append(sorted.at(0));

    for (const auto &atom : sorted)
    {
        if (atom != atms.last())
            atms.append(atom);
    }

    atms.squeeze();
}

/** Copy constructor */
PositionalRestraint::PositionalRestraint(const PositionalRestraint &other)
    : ConcreteProperty<PositionalRestraint, Property>(other),
      atms(other.atms), pos(other.pos), _k(other._k), _r0(other._r0)
{
}

PositionalRestraint::~PositionalRestraint()
{
}

PositionalRestraint &PositionalRestraint::operator=(const PositionalRestraint &other)
{
    if (this != &other)
    {
        atms = other.atms;
        pos = other.pos;
        _k = other._k;
        _r0 = other._r0;
    }

    return *this;
}

bool PositionalRestraint::operator==(const PositionalRestraint &other) const
{
    return atms == other.atms and pos == other.pos and
           _k == other._k and _r0 == other._r0;
}

bool PositionalRestraint::operator!=(const PositionalRestraint &other) const
{
    return not operator==(other);
}

PositionalRestraints PositionalRestraint::operator+(const PositionalRestraint &other) const
{
    return PositionalRestraints(*this) + other;
}

PositionalRestraints PositionalRestraint::operator+(const PositionalRestraints &other) const
{
    return PositionalRestraints(*this) + other;
}

const char *PositionalRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PositionalRestraint>());
}

const char *PositionalRestraint::what() const
{
    return PositionalRestraint::typeName();
}

PositionalRestraint *PositionalRestraint::clone() const
{
    return new PositionalRestraint(*this);
}

bool PositionalRestraint::isNull() const
{
    return atms.isEmpty();
}

QString PositionalRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("PositionalRestraint::null");

    else if (this->isAtomRestraint())
    {
        return QString("PositionalRestraint( %1 => %2, k=%3 : r0=%4 )")
            .arg(this->atom())
            .arg(pos.toString())
            .arg(_k.toString())
            .arg(_r0.toString());
    }
    else
    {
        QStringList a;

        for (const auto &atom : atms)
        {
            a.append(QString::number(atom));
        }

        return QString("PositionalRestraint( [%1] => %2, k=%3 : r0=%4 )")
            .arg(a.join(", "))
            .arg(pos.toString())
            .arg(_k.toString())
            .arg(_r0.toString());
    }
}

/** Return whether this is a single-atom restraint */
bool PositionalRestraint::isAtomRestraint() const
{
    return atms.count() == 1;
}

/** Return whether this restraint acts on the centroid of a group
 *  of atoms */
bool PositionalRestraint::isCentroidRestraint() const
{
    return atms.count() > 1;
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 PositionalRestraint::atom() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms.at(0);
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> PositionalRestraint::atoms() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms;
}

/** Return the position in space where the restraint is placed */
SireMaths::Vector PositionalRestraint::position() const
{
    return this->pos;
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicBondConstant PositionalRestraint::k() const
{
    return this->_k;
}

/** Return the width of the flat-bottom well. This is zero for a
 *  pure harmonic restraint
 */
SireUnits::Dimension::Length PositionalRestraint::r0() const
{
    return this->_r0;
}

///////
/////// Implementation of PositionalRestraints
///////

static const RegisterMetaType<PositionalRestraints> r_posrests;

QDataStream &operator<<(QDataStream &ds, const PositionalRestraints &posrests)
{
    writeHeader(ds, r_posrests, 2);

    SharedDataStream sds(ds);

    sds << posrests.r << posrests.use_pbc
        << static_cast<const Restraints &>(posrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PositionalRestraints &posrests)
{
    VersionID v = readHeader(ds, r_posrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> posrests.r >> static_cast<Restraints &>(posrests);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> posrests.r >> posrests.use_pbc
            >> static_cast<Restraints &>(posrests);
    }
    else
        throw version_error(v, "1,2", r_posrests, CODELOC);

    return ds;
}

/** Null constructor */
PositionalRestraints::PositionalRestraints()
    : ConcreteProperty<PositionalRestraints, Restraints>()
{
}

PositionalRestraints::PositionalRestraints(const QString &name)
    : ConcreteProperty<PositionalRestraints, Restraints>(name)
{
}

PositionalRestraints::PositionalRestraints(const PositionalRestraint &restraint)
    : ConcreteProperty<PositionalRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

PositionalRestraints::PositionalRestraints(const QList<PositionalRestraint> &restraints)
    : ConcreteProperty<PositionalRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

PositionalRestraints::PositionalRestraints(const QString &name,
                                           const PositionalRestraint &restraint)
    : ConcreteProperty<PositionalRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

PositionalRestraints::PositionalRestraints(const QString &name,
                                           const QList<PositionalRestraint> &restraints)
    : ConcreteProperty<PositionalRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

PositionalRestraints::PositionalRestraints(const PositionalRestraints &other)
    : ConcreteProperty<PositionalRestraints, Restraints>(other), r(other.r), use_pbc(other.use_pbc)
{
}

PositionalRestraints::~PositionalRestraints()
{
}

PositionalRestraints &PositionalRestraints::operator=(const PositionalRestraints &other)
{
    r = other.r;
    use_pbc = other.use_pbc;
    Restraints::operator=(other);
    return *this;
}

bool PositionalRestraints::operator==(const PositionalRestraints &other) const
{
    return r == other.r and Restraints::operator==(other) and use_pbc == other.use_pbc;
}

bool PositionalRestraints::operator!=(const PositionalRestraints &other) const
{
    return not operator==(other);
}

const char *PositionalRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PositionalRestraints>());
}

const char *PositionalRestraints::what() const
{
    return PositionalRestraints::typeName();
}

PositionalRestraints *PositionalRestraints::clone() const
{
    return new PositionalRestraints(*this);
}

QString PositionalRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("PositionalRestraints::null");

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

    return QObject::tr("PositionalRestraints( name=%1, size=%2, use_pbc=%3\n%4\n)")
            .arg(this->name())
            .arg(n)
            .arg(this->use_pbc ? "true" : "false")
            .arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool PositionalRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool PositionalRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int PositionalRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int PositionalRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int PositionalRestraints::size() const
{
    return this->nRestraints();
}

/** Return the number of atom restraints */
int PositionalRestraints::nAtomRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isAtomRestraint());
    }

    return n;
}

/** Return the number of centroid restraints */
int PositionalRestraints::nCentroidRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isCentroidRestraint());
    }

    return n;
}

/** Return whether or not there are any atom restraints */
bool PositionalRestraints::hasAtomRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isAtomRestraint())
            return true;
    }

    return false;
}

/** Return whether or not there are any centroid restraints */
bool PositionalRestraints::hasCentroidRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isCentroidRestraint())
            return true;
    }

    return false;
}

/** Return the ith restraint */
const PositionalRestraint &PositionalRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const PositionalRestraint &PositionalRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<PositionalRestraint> PositionalRestraints::restraints() const
{
    return this->r;
}

/** Return all of the atom restraints */
QList<PositionalRestraint> PositionalRestraints::atomRestraints() const
{
    if (this->hasCentroidRestraints())
    {
        QList<PositionalRestraint> ar;

        for (const auto &restraint : this->r)
        {
            if (restraint.isAtomRestraint())
                ar.append(restraint);
        }

        return ar;
    }
    else
        return this->restraints();
}

/** Return all of the centroid restraints */
QList<PositionalRestraint> PositionalRestraints::centroidRestraints() const
{
    if (this->hasAtomRestraints())
    {
        QList<PositionalRestraint> cr;

        for (const auto &restraint : this->r)
        {
            if (restraint.isCentroidRestraint())
                cr.append(restraint);
        }

        return cr;
    }
    else
        return this->restraints();
}

/** Add a restraint onto the list */
void PositionalRestraints::add(const PositionalRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void PositionalRestraints::add(const PositionalRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
PositionalRestraints &PositionalRestraints::operator+=(const PositionalRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
PositionalRestraints PositionalRestraints::operator+(const PositionalRestraint &restraint) const
{
    PositionalRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
PositionalRestraints &PositionalRestraints::operator+=(const PositionalRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
PositionalRestraints PositionalRestraints::operator+(const PositionalRestraints &restraints) const
{
    PositionalRestraints ret(*this);
    ret += restraints;
    return *this;
}

/** Set whether or not periodic boundary conditions are to be used */
void PositionalRestraints::setUsesPeriodicBoundaryConditions(bool use_pbc)
{
    this->use_pbc = use_pbc;
}

/** Return whether or not periodic boundary conditions are to be used */
bool PositionalRestraints::getUsesPeriodicBoundaryConditions() const
{
    return this->use_pbc;
}
