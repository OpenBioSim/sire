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

#include "bondrestraints.h"

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
/////// Implementation of BondRestraint
///////

static const RegisterMetaType<BondRestraint> r_bndrest;

QDataStream &operator<<(QDataStream &ds, const BondRestraint &bndrest)
{
    writeHeader(ds, r_bndrest, 1);

    SharedDataStream sds(ds);

    sds << bndrest.atms0 << bndrest.atms1 << bndrest._k << bndrest._r0
        << static_cast<const Property &>(bndrest);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BondRestraint &bndrest)
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
BondRestraint::BondRestraint()
    : ConcreteProperty<BondRestraint, Property>(),
      _k(0), _r0(0)
{
}

/** Construct to restrain the atom at index 'atom' to the specified position
 *  using the specified force constant and flat-bottom well-width
 */
BondRestraint::BondRestraint(qint64 atom0, qint64 atom1,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0)
    : ConcreteProperty<BondRestraint, Property>(),
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
BondRestraint::BondRestraint(const QList<qint64> &atoms0,
                             const QList<qint64> &atoms1,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0)
    : ConcreteProperty<BondRestraint, Property>(),
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
BondRestraint::BondRestraint(const BondRestraint &other)
    : ConcreteProperty<BondRestraint, Property>(other),
      atms0(other.atms0), atms1(other.atms1), _k(other._k), _r0(other._r0)
{
}

BondRestraint::~BondRestraint()
{
}

BondRestraint &BondRestraint::operator=(const BondRestraint &other)
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

bool BondRestraint::operator==(const BondRestraint &other) const
{
    return atms0 == other.atms0 and atms1 == other.atms1 and
           _k == other._k and _r0 == other._r0;
}

bool BondRestraint::operator!=(const BondRestraint &other) const
{
    return not operator==(other);
}

BondRestraints BondRestraint::operator+(const BondRestraint &other) const
{
    return BondRestraints(*this) + other;
}

BondRestraints BondRestraint::operator+(const BondRestraints &other) const
{
    return BondRestraints(*this) + other;
}

const char *BondRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<BondRestraint>());
}

const char *BondRestraint::what() const
{
    return BondRestraint::typeName();
}

BondRestraint *BondRestraint::clone() const
{
    return new BondRestraint(*this);
}

bool BondRestraint::isNull() const
{
    return atms0.isEmpty() or atms1.isEmpty();
}

QString BondRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("BondRestraint::null");

    else if (this->isAtomRestraint())
    {
        return QString("BondRestraint( %1 <=> %2, k=%3 : r0=%4 )")
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

        return QString("BondRestraint( [%1] <=> [%2], k=%3 : r0=%4 )")
            .arg(a0.join(", "))
            .arg(a1.join(", "))
            .arg(_k.toString())
            .arg(_r0.toString());
    }
}

/** Return whether this is a single-atom restraint */
bool BondRestraint::isAtomRestraint() const
{
    return atms0.count() == 1 and atms1.count() == 1;
}

/** Return whether this restraint acts on the centroid of a group
 *  of atoms */
bool BondRestraint::isCentroidRestraint() const
{
    return atms0.count() > 1 or atms1.count() > 1;
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 BondRestraint::atom0() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms0.at(0);
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 BondRestraint::atom1() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms1.at(0);
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> BondRestraint::atoms0() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms0;
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> BondRestraint::atoms1() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms1;
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicBondConstant BondRestraint::k() const
{
    return this->_k;
}

/** Return the width of the harmonic bond. */
SireUnits::Dimension::Length BondRestraint::r0() const
{
    return this->_r0;
}

///////
/////// Implementation of BondRestraints
///////

static const RegisterMetaType<BondRestraints> r_bndrests;

QDataStream &operator<<(QDataStream &ds, const BondRestraints &bndrests)
{
    writeHeader(ds, r_bndrests, 2);

    SharedDataStream sds(ds);

    sds << bndrests.r << bndrests.use_pbc
        << static_cast<const Restraints &>(bndrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BondRestraints &bndrests)
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
        throw version_error(v, "1,2", r_bndrests, CODELOC);

    return ds;
}

/** Null constructor */
BondRestraints::BondRestraints()
    : ConcreteProperty<BondRestraints, Restraints>()
{
}

BondRestraints::BondRestraints(const QString &name)
    : ConcreteProperty<BondRestraints, Restraints>(name)
{
}

BondRestraints::BondRestraints(const BondRestraint &restraint)
    : ConcreteProperty<BondRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

BondRestraints::BondRestraints(const QList<BondRestraint> &restraints)
    : ConcreteProperty<BondRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

BondRestraints::BondRestraints(const QString &name,
                               const BondRestraint &restraint)
    : ConcreteProperty<BondRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

BondRestraints::BondRestraints(const QString &name,
                               const QList<BondRestraint> &restraints)
    : ConcreteProperty<BondRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

BondRestraints::BondRestraints(const BondRestraints &other)
    : ConcreteProperty<BondRestraints, Restraints>(other), r(other.r), use_pbc(other.use_pbc)
{
}

BondRestraints::~BondRestraints()
{
}

BondRestraints &BondRestraints::operator=(const BondRestraints &other)
{
    r = other.r;
    use_pbc = other.use_pbc;
    Restraints::operator=(other);
    return *this;
}

bool BondRestraints::operator==(const BondRestraints &other) const
{
    return r == other.r and Restraints::operator==(other) and use_pbc == other.use_pbc;
}

bool BondRestraints::operator!=(const BondRestraints &other) const
{
    return not operator==(other);
}

const char *BondRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<BondRestraints>());
}

const char *BondRestraints::what() const
{
    return BondRestraints::typeName();
}

BondRestraints *BondRestraints::clone() const
{
    return new BondRestraints(*this);
}

QString BondRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("BondRestraints::null");

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

    return QObject::tr("BondRestraints( name=%1, size=%2, use_pbc=$3\n%4\n)")
            .arg(this->name())
            .arg(n)
            .arg(this->use_pbc ? "true" : "false")
            .arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool BondRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool BondRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int BondRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int BondRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int BondRestraints::size() const
{
    return this->nRestraints();
}

/** Return the number of atom restraints */
int BondRestraints::nAtomRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isAtomRestraint());
    }

    return n;
}

/** Return the number of centroid restraints */
int BondRestraints::nCentroidRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isCentroidRestraint());
    }

    return n;
}

/** Return whether or not there are any atom restraints */
bool BondRestraints::hasAtomRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isAtomRestraint())
            return true;
    }

    return false;
}

/** Return whether or not there are any centroid restraints */
bool BondRestraints::hasCentroidRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isCentroidRestraint())
            return true;
    }

    return false;
}

/** Return the ith restraint */
const BondRestraint &BondRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const BondRestraint &BondRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<BondRestraint> BondRestraints::restraints() const
{
    return this->r;
}

/** Return all of the atom restraints */
QList<BondRestraint> BondRestraints::atomRestraints() const
{
    if (this->hasCentroidRestraints())
    {
        QList<BondRestraint> ar;

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
QList<BondRestraint> BondRestraints::centroidRestraints() const
{
    if (this->hasAtomRestraints())
    {
        QList<BondRestraint> cr;

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
void BondRestraints::add(const BondRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void BondRestraints::add(const BondRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
BondRestraints &BondRestraints::operator+=(const BondRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
BondRestraints BondRestraints::operator+(const BondRestraint &restraint) const
{
    BondRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
BondRestraints &BondRestraints::operator+=(const BondRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
BondRestraints BondRestraints::operator+(const BondRestraints &restraints) const
{
    BondRestraints ret(*this);
    ret += restraints;
    return *this;
}

/** Set whether or not periodic boundary conditions are to be used */
void BondRestraints::setUsesPbc(bool use_pbc)
{
    this->use_pbc = use_pbc;
}

/** Return whether or not periodic boundary conditions are to be used */
bool BondRestraints::usesPbc() const
{
    return this->use_pbc;
}
