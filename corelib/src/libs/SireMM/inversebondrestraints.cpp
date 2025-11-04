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

#include "inversebondrestraints.h"

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
/////// Implementation of InverseBondRestraint
///////

static const RegisterMetaType<InverseBondRestraint> r_bndrest;

QDataStream &operator<<(QDataStream &ds, const InverseBondRestraint &bndrest)
{
    writeHeader(ds, r_bndrest, 1);

    SharedDataStream sds(ds);

    sds << bndrest.atms0 << bndrest.atms1 << bndrest._k << bndrest._r0
        << static_cast<const Property &>(bndrest);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, InverseBondRestraint &bndrest)
{
    VersionID v = readHeader(ds, r_bndrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> bndrest.atms0 >> bndrest.atms1 >> bndrest._k >> bndrest._r0 >> static_cast<Property &>(bndrest);
    }
    else
        throw version_error(v, "1", r_bndrest, CODELOC);

    return ds;
}

/** Null constructor */
InverseBondRestraint::InverseBondRestraint()
    : ConcreteProperty<InverseBondRestraint, Property>(),
      _k(0), _r0(0)
{
}

/** Construct to restrain the atom at index 'atom' to the specified position
 *  using the specified force constant and flat-bottom well-width
 */
InverseBondRestraint::InverseBondRestraint(qint64 atom0, qint64 atom1,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0)
    : ConcreteProperty<InverseBondRestraint, Property>(),
      _k(k), _r0(r0)
{
    if (atom0 == atom1)
        throw SireError::invalid_arg(QObject::tr(
                                         "You cannot create a bond restraint between identical atoms! %1-%2")
                                         .arg(atom0)
                                         .arg(atom1),
                                     CODELOC);

    atms0 = QVector<qint64>(1, atom0);
    atms0.squeeze();

    atms1 = QVector<qint64>(1, atom1);
    atms1.squeeze();
}

/** Construct to restrain the centroid of the atoms whose indicies are
 *  in 'atoms' to the specified position using the specified force constant
 *  and flat-bottom well width
 */
InverseBondRestraint::InverseBondRestraint(const QList<qint64> &atoms0,
                             const QList<qint64> &atoms1,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0)
    : ConcreteProperty<InverseBondRestraint, Property>(),
      _k(k), _r0(r0)
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    // remove duplicates
    atms0.reserve(atoms0.count());

    auto sorted = atoms0;
    std::sort(sorted.begin(), sorted.end());

    atms0.append(sorted.at(0));

    for (const auto &atom : sorted)
    {
        if (atom != atms0.last())
            atms0.append(atom);
    }

    atms0.squeeze();

    // now the same for atoms1
    atms1.reserve(atoms1.count());

    sorted = atoms1;
    std::sort(sorted.begin(), sorted.end());

    if (atms0.indexOf(sorted.at(0)) != -1)
        throw SireError::invalid_arg(QObject::tr(
                                         "You cannot have an overlap in atoms between the two groups. Atom "
                                         "%1 appears in both!")
                                         .arg(sorted.at(0)),
                                     CODELOC);

    atms1.append(sorted.at(0));

    for (const auto &atom : sorted)
    {
        if (atom != atms1.last())
        {
            if (atms0.indexOf(atom) != -1)
                throw SireError::invalid_arg(QObject::tr(
                                                 "You cannot have an overlap in atoms between the two groups. Atom "
                                                 "%1 appears in both!")
                                                 .arg(sorted.at(0)),
                                             CODELOC);

            atms1.append(atom);
        }

        atms1.squeeze();
    }
}

/** Copy constructor */
InverseBondRestraint::InverseBondRestraint(const InverseBondRestraint &other)
    : ConcreteProperty<InverseBondRestraint, Property>(other),
      atms0(other.atms0), atms1(other.atms1), _k(other._k), _r0(other._r0)
{
}

InverseBondRestraint::~InverseBondRestraint()
{
}

InverseBondRestraint &InverseBondRestraint::operator=(const InverseBondRestraint &other)
{
    if (this != &other)
    {
        atms0 = other.atms0;
        atms1 = other.atms1;
        _k = other._k;
        _r0 = other._r0;
    }

    return *this;
}

bool InverseBondRestraint::operator==(const InverseBondRestraint &other) const
{
    return atms0 == other.atms0 and atms1 == other.atms1 and
           _k == other._k and _r0 == other._r0;
}

bool InverseBondRestraint::operator!=(const InverseBondRestraint &other) const
{
    return not operator==(other);
}

InverseBondRestraints InverseBondRestraint::operator+(const InverseBondRestraint &other) const
{
    return InverseBondRestraints(*this) + other;
}

InverseBondRestraints InverseBondRestraint::operator+(const InverseBondRestraints &other) const
{
    return InverseBondRestraints(*this) + other;
}

const char *InverseBondRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<InverseBondRestraint>());
}

const char *InverseBondRestraint::what() const
{
    return InverseBondRestraint::typeName();
}

InverseBondRestraint *InverseBondRestraint::clone() const
{
    return new InverseBondRestraint(*this);
}

bool InverseBondRestraint::isNull() const
{
    return atms0.isEmpty() or atms1.isEmpty();
}

QString InverseBondRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("InverseBondRestraint::null");

    else if (this->isAtomRestraint())
    {
        return QString("InverseBondRestraint( %1 <=> %2, k=%3 : r0=%4 )")
            .arg(this->atom0())
            .arg(this->atom1())
            .arg(_k.toString())
            .arg(_r0.toString());
    }
    else
    {
        QStringList a0, a1;

        for (const auto &atom : atms0)
        {
            a0.append(QString::number(atom));
        }

        for (const auto &atom : atms1)
        {
            a1.append(QString::number(atom));
        }

        return QString("InverseBondRestraint( [%1] <=> [%2], k=%3 : r0=%4 )")
            .arg(a0.join(", "))
            .arg(a1.join(", "))
            .arg(_k.toString())
            .arg(_r0.toString());
    }
}

/** Return whether this is a single-atom restraint */
bool InverseBondRestraint::isAtomRestraint() const
{
    return atms0.count() == 1 and atms1.count() == 1;
}

/** Return whether this restraint acts on the centroid of a group
 *  of atoms */
bool InverseBondRestraint::isCentroidRestraint() const
{
    return atms0.count() > 1 or atms1.count() > 1;
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 InverseBondRestraint::atom0() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms0.at(0);
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 InverseBondRestraint::atom1() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms1.at(0);
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> InverseBondRestraint::atoms0() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms0;
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> InverseBondRestraint::atoms1() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms1;
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicBondConstant InverseBondRestraint::k() const
{
    return this->_k;
}

/** Return the width of the harmonic bond. */
SireUnits::Dimension::Length InverseBondRestraint::r0() const
{
    return this->_r0;
}

///////
/////// Implementation of InverseBondRestraints
///////

static const RegisterMetaType<InverseBondRestraints> r_bndrests;

QDataStream &operator<<(QDataStream &ds, const InverseBondRestraints &bndrests)
{
    writeHeader(ds, r_bndrests, 2);

    SharedDataStream sds(ds);

    sds << bndrests.r << bndrests.use_pbc
        << static_cast<const Restraints &>(bndrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, InverseBondRestraints &bndrests)
{
    VersionID v = readHeader(ds, r_bndrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> bndrests.r >> static_cast<Restraints &>(bndrests);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> bndrests.r >> bndrests.use_pbc
            >> static_cast<Restraints &>(bndrests);
    }
    else
        throw version_error(v, "1", r_bndrests, CODELOC);

    return ds;
}

/** Null constructor */
InverseBondRestraints::InverseBondRestraints()
    : ConcreteProperty<InverseBondRestraints, Restraints>()
{
}

InverseBondRestraints::InverseBondRestraints(const QString &name)
    : ConcreteProperty<InverseBondRestraints, Restraints>(name)
{
}

InverseBondRestraints::InverseBondRestraints(const InverseBondRestraint &restraint)
    : ConcreteProperty<InverseBondRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

InverseBondRestraints::InverseBondRestraints(const QList<InverseBondRestraint> &restraints)
    : ConcreteProperty<InverseBondRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

InverseBondRestraints::InverseBondRestraints(const QString &name,
                               const InverseBondRestraint &restraint)
    : ConcreteProperty<InverseBondRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

InverseBondRestraints::InverseBondRestraints(const QString &name,
                               const QList<InverseBondRestraint> &restraints)
    : ConcreteProperty<InverseBondRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

InverseBondRestraints::InverseBondRestraints(const InverseBondRestraints &other)
    : ConcreteProperty<InverseBondRestraints, Restraints>(other), r(other.r), use_pbc(other.use_pbc)
{
}

InverseBondRestraints::~InverseBondRestraints()
{
}

InverseBondRestraints &InverseBondRestraints::operator=(const InverseBondRestraints &other)
{
    r = other.r;
    use_pbc = other.use_pbc;
    Restraints::operator=(other);
    return *this;
}

bool InverseBondRestraints::operator==(const InverseBondRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
}

bool InverseBondRestraints::operator!=(const InverseBondRestraints &other) const
{
    return not operator==(other);
}

const char *InverseBondRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<InverseBondRestraints>());
}

const char *InverseBondRestraints::what() const
{
    return InverseBondRestraints::typeName();
}

InverseBondRestraints *InverseBondRestraints::clone() const
{
    return new InverseBondRestraints(*this);
}

QString InverseBondRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("InverseBondRestraints::null");

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

    return QObject::tr("InverseBondRestraints( name=%1, size=%2, use_pbc=%3\n%4\n)")
            .arg(this->name())
            .arg(n)
            .arg(this->use_pbc ? "true" : "false")
            .arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool InverseBondRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool InverseBondRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int InverseBondRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int InverseBondRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int InverseBondRestraints::size() const
{
    return this->nRestraints();
}

/** Return the number of atom restraints */
int InverseBondRestraints::nAtomRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isAtomRestraint());
    }

    return n;
}

/** Return the number of centroid restraints */
int InverseBondRestraints::nCentroidRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isCentroidRestraint());
    }

    return n;
}

/** Return whether or not there are any atom restraints */
bool InverseBondRestraints::hasAtomRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isAtomRestraint())
            return true;
    }

    return false;
}

/** Return whether or not there are any centroid restraints */
bool InverseBondRestraints::hasCentroidRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isCentroidRestraint())
            return true;
    }

    return false;
}

/** Return the ith restraint */
const InverseBondRestraint &InverseBondRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const InverseBondRestraint &InverseBondRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<InverseBondRestraint> InverseBondRestraints::restraints() const
{
    return this->r;
}

/** Return all of the atom restraints */
QList<InverseBondRestraint> InverseBondRestraints::atomRestraints() const
{
    if (this->hasCentroidRestraints())
    {
        QList<InverseBondRestraint> ar;

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
QList<InverseBondRestraint> InverseBondRestraints::centroidRestraints() const
{
    if (this->hasAtomRestraints())
    {
        QList<InverseBondRestraint> cr;

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
void InverseBondRestraints::add(const InverseBondRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void InverseBondRestraints::add(const InverseBondRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
InverseBondRestraints &InverseBondRestraints::operator+=(const InverseBondRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
InverseBondRestraints InverseBondRestraints::operator+(const InverseBondRestraint &restraint) const
{
    InverseBondRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
InverseBondRestraints &InverseBondRestraints::operator+=(const InverseBondRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
InverseBondRestraints InverseBondRestraints::operator+(const InverseBondRestraints &restraints) const
{
    InverseBondRestraints ret(*this);
    ret += restraints;
    return *this;
}

/** Set whether or not periodic boundary conditions are to be used */
void InverseBondRestraints::setUsesPeriodicBoundaryConditions(bool use_pbc)
{
    this->use_pbc = use_pbc;
}

/** Return whether or not periodic boundary conditions are to be used */
bool InverseBondRestraints::getUsesPeriodicBoundaryConditions() const
{
    return this->use_pbc;
}
