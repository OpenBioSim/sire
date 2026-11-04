/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2026  Lester Hedges
  *
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published
  *  by the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  *
  *  For full details of the license please see the LICENSE file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the website
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#include "flatbottombondrestraints.h"

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
/////// Implementation of FlatBottomBondRestraint
///////

static const RegisterMetaType<FlatBottomBondRestraint> r_fbbndrest;

QDataStream &operator<<(QDataStream &ds, const FlatBottomBondRestraint &fbbndrest)
{
    writeHeader(ds, r_fbbndrest, 1);

    SharedDataStream sds(ds);

    sds << fbbndrest.atms0 << fbbndrest.atms1 << fbbndrest._k
        << fbbndrest._r0 << fbbndrest._rfb
        << static_cast<const Property &>(fbbndrest);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, FlatBottomBondRestraint &fbbndrest)
{
    VersionID v = readHeader(ds, r_fbbndrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> fbbndrest.atms0 >> fbbndrest.atms1 >> fbbndrest._k >> fbbndrest._r0 >> fbbndrest._rfb >> static_cast<Property &>(fbbndrest);
    }
    else
        throw version_error(v, "1", r_fbbndrest, CODELOC);

    return ds;
}

static void assert_valid_rfb(const SireUnits::Dimension::Length &rfb)
{
    if (rfb.value() < 0)
        throw SireError::invalid_arg(QObject::tr(
                                         "The flat-bottom radius cannot be negative! %1")
                                         .arg(rfb.toString()),
                                     CODELOC);
}

/** Null constructor */
FlatBottomBondRestraint::FlatBottomBondRestraint()
    : ConcreteProperty<FlatBottomBondRestraint, Property>(),
      _k(0), _r0(0), _rfb(0)
{
}

/** Construct to restrain the distance between atoms 'atom0' and 'atom1'
 *  using the specified force constant, equilibrium distance and
 *  flat-bottom radius
 */
FlatBottomBondRestraint::FlatBottomBondRestraint(qint64 atom0, qint64 atom1,
                                                 const SireUnits::Dimension::HarmonicBondConstant &k,
                                                 const SireUnits::Dimension::Length &r0,
                                                 const SireUnits::Dimension::Length &rfb)
    : ConcreteProperty<FlatBottomBondRestraint, Property>(),
      _k(k), _r0(r0), _rfb(rfb)
{
    if (atom0 == atom1)
        throw SireError::invalid_arg(QObject::tr(
                                         "You cannot create a bond restraint between identical atoms! %1-%2")
                                         .arg(atom0)
                                         .arg(atom1),
                                     CODELOC);

    assert_valid_rfb(rfb);

    atms0 = QVector<qint64>(1, atom0);
    atms0.squeeze();

    atms1 = QVector<qint64>(1, atom1);
    atms1.squeeze();
}

/** Construct to restrain the distance between the centroids of the atoms
 *  whose indices are in 'atoms0' and 'atoms1' using the specified force
 *  constant, equilibrium distance and flat-bottom radius
 */
FlatBottomBondRestraint::FlatBottomBondRestraint(const QList<qint64> &atoms0,
                                                 const QList<qint64> &atoms1,
                                                 const SireUnits::Dimension::HarmonicBondConstant &k,
                                                 const SireUnits::Dimension::Length &r0,
                                                 const SireUnits::Dimension::Length &rfb)
    : ConcreteProperty<FlatBottomBondRestraint, Property>(),
      _k(k), _r0(r0), _rfb(rfb)
{
    assert_valid_rfb(rfb);

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

    for (const auto &atom : sorted)
    {
        if (atms1.isEmpty() or atom != atms1.last())
        {
            if (atms0.indexOf(atom) != -1)
                throw SireError::invalid_arg(QObject::tr(
                                                 "You cannot have an overlap in atoms between the two groups. Atom "
                                                 "%1 appears in both!")
                                                 .arg(atom),
                                             CODELOC);

            atms1.append(atom);
        }
    }

    atms1.squeeze();
}

/** Copy constructor */
FlatBottomBondRestraint::FlatBottomBondRestraint(const FlatBottomBondRestraint &other)
    : ConcreteProperty<FlatBottomBondRestraint, Property>(other),
      atms0(other.atms0), atms1(other.atms1),
      _k(other._k), _r0(other._r0), _rfb(other._rfb)
{
}

FlatBottomBondRestraint::~FlatBottomBondRestraint()
{
}

FlatBottomBondRestraint &FlatBottomBondRestraint::operator=(const FlatBottomBondRestraint &other)
{
    if (this != &other)
    {
        Property::operator=(other);
        atms0 = other.atms0;
        atms1 = other.atms1;
        _k = other._k;
        _r0 = other._r0;
        _rfb = other._rfb;
    }

    return *this;
}

bool FlatBottomBondRestraint::operator==(const FlatBottomBondRestraint &other) const
{
    return atms0 == other.atms0 and atms1 == other.atms1 and
           _k == other._k and _r0 == other._r0 and _rfb == other._rfb;
}

bool FlatBottomBondRestraint::operator!=(const FlatBottomBondRestraint &other) const
{
    return not operator==(other);
}

FlatBottomBondRestraints FlatBottomBondRestraint::operator+(const FlatBottomBondRestraint &other) const
{
    return FlatBottomBondRestraints(*this) + other;
}

FlatBottomBondRestraints FlatBottomBondRestraint::operator+(const FlatBottomBondRestraints &other) const
{
    return FlatBottomBondRestraints(*this) + other;
}

const char *FlatBottomBondRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<FlatBottomBondRestraint>());
}

const char *FlatBottomBondRestraint::what() const
{
    return FlatBottomBondRestraint::typeName();
}

FlatBottomBondRestraint *FlatBottomBondRestraint::clone() const
{
    return new FlatBottomBondRestraint(*this);
}

bool FlatBottomBondRestraint::isNull() const
{
    return atms0.isEmpty() or atms1.isEmpty();
}

QString FlatBottomBondRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("FlatBottomBondRestraint::null");

    else if (this->isAtomRestraint())
    {
        return QString("FlatBottomBondRestraint( %1 <=> %2, k=%3 : r0=%4 : rfb=%5 )")
            .arg(this->atom0())
            .arg(this->atom1())
            .arg(_k.toString())
            .arg(_r0.toString())
            .arg(_rfb.toString());
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

        return QString("FlatBottomBondRestraint( [%1] <=> [%2], k=%3 : r0=%4 : rfb=%5 )")
            .arg(a0.join(", "))
            .arg(a1.join(", "))
            .arg(_k.toString())
            .arg(_r0.toString())
            .arg(_rfb.toString());
    }
}

/** Return whether this is a single-atom restraint */
bool FlatBottomBondRestraint::isAtomRestraint() const
{
    return atms0.count() == 1 and atms1.count() == 1;
}

/** Return whether this restraint acts on the centroid of a group
 *  of atoms */
bool FlatBottomBondRestraint::isCentroidRestraint() const
{
    return atms0.count() > 1 or atms1.count() > 1;
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 FlatBottomBondRestraint::atom0() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms0.at(0);
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 FlatBottomBondRestraint::atom1() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms1.at(0);
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> FlatBottomBondRestraint::atoms0() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms0;
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> FlatBottomBondRestraint::atoms1() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms1;
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicBondConstant FlatBottomBondRestraint::k() const
{
    return this->_k;
}

/** Return the equilibrium distance for the restraint */
SireUnits::Dimension::Length FlatBottomBondRestraint::r0() const
{
    return this->_r0;
}

/** Return the flat-bottom radius about the equilibrium distance */
SireUnits::Dimension::Length FlatBottomBondRestraint::rfb() const
{
    return this->_rfb;
}

///////
/////// Implementation of FlatBottomBondRestraints
///////

static const RegisterMetaType<FlatBottomBondRestraints> r_fbbndrests;

QDataStream &operator<<(QDataStream &ds, const FlatBottomBondRestraints &fbbndrests)
{
    writeHeader(ds, r_fbbndrests, 1);

    SharedDataStream sds(ds);

    sds << fbbndrests.r << fbbndrests.use_pbc
        << static_cast<const Restraints &>(fbbndrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, FlatBottomBondRestraints &fbbndrests)
{
    VersionID v = readHeader(ds, r_fbbndrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> fbbndrests.r >> fbbndrests.use_pbc >> static_cast<Restraints &>(fbbndrests);
    }
    else
        throw version_error(v, "1", r_fbbndrests, CODELOC);

    return ds;
}

/** Null constructor */
FlatBottomBondRestraints::FlatBottomBondRestraints()
    : ConcreteProperty<FlatBottomBondRestraints, Restraints>()
{
}

FlatBottomBondRestraints::FlatBottomBondRestraints(const QString &name)
    : ConcreteProperty<FlatBottomBondRestraints, Restraints>(name)
{
}

FlatBottomBondRestraints::FlatBottomBondRestraints(const FlatBottomBondRestraint &restraint)
    : ConcreteProperty<FlatBottomBondRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

FlatBottomBondRestraints::FlatBottomBondRestraints(const QList<FlatBottomBondRestraint> &restraints)
    : ConcreteProperty<FlatBottomBondRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

FlatBottomBondRestraints::FlatBottomBondRestraints(const QString &name,
                                                   const FlatBottomBondRestraint &restraint)
    : ConcreteProperty<FlatBottomBondRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

FlatBottomBondRestraints::FlatBottomBondRestraints(const QString &name,
                                                   const QList<FlatBottomBondRestraint> &restraints)
    : ConcreteProperty<FlatBottomBondRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

FlatBottomBondRestraints::FlatBottomBondRestraints(const FlatBottomBondRestraints &other)
    : ConcreteProperty<FlatBottomBondRestraints, Restraints>(other), r(other.r), use_pbc(other.use_pbc)
{
}

FlatBottomBondRestraints::~FlatBottomBondRestraints()
{
}

FlatBottomBondRestraints &FlatBottomBondRestraints::operator=(const FlatBottomBondRestraints &other)
{
    r = other.r;
    use_pbc = other.use_pbc;
    Restraints::operator=(other);
    return *this;
}

bool FlatBottomBondRestraints::operator==(const FlatBottomBondRestraints &other) const
{
    return r == other.r and Restraints::operator==(other) and use_pbc == other.use_pbc;
}

bool FlatBottomBondRestraints::operator!=(const FlatBottomBondRestraints &other) const
{
    return not operator==(other);
}

const char *FlatBottomBondRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<FlatBottomBondRestraints>());
}

const char *FlatBottomBondRestraints::what() const
{
    return FlatBottomBondRestraints::typeName();
}

FlatBottomBondRestraints *FlatBottomBondRestraints::clone() const
{
    return new FlatBottomBondRestraints(*this);
}

QString FlatBottomBondRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("FlatBottomBondRestraints::null");

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

    return QObject::tr("FlatBottomBondRestraints( name=%1, size=%2, use_pbc=%3\n%4\n)")
        .arg(this->name())
        .arg(n)
        .arg(this->use_pbc ? "true" : "false")
        .arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool FlatBottomBondRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool FlatBottomBondRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int FlatBottomBondRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int FlatBottomBondRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int FlatBottomBondRestraints::size() const
{
    return this->nRestraints();
}

/** Return the number of atom restraints */
int FlatBottomBondRestraints::nAtomRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isAtomRestraint());
    }

    return n;
}

/** Return the number of centroid restraints */
int FlatBottomBondRestraints::nCentroidRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isCentroidRestraint());
    }

    return n;
}

/** Return whether or not there are any atom restraints */
bool FlatBottomBondRestraints::hasAtomRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isAtomRestraint())
            return true;
    }

    return false;
}

/** Return whether or not there are any centroid restraints */
bool FlatBottomBondRestraints::hasCentroidRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isCentroidRestraint())
            return true;
    }

    return false;
}

/** Return the ith restraint */
const FlatBottomBondRestraint &FlatBottomBondRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const FlatBottomBondRestraint &FlatBottomBondRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<FlatBottomBondRestraint> FlatBottomBondRestraints::restraints() const
{
    return this->r;
}

/** Return all of the atom restraints */
QList<FlatBottomBondRestraint> FlatBottomBondRestraints::atomRestraints() const
{
    if (this->hasCentroidRestraints())
    {
        QList<FlatBottomBondRestraint> ar;

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
QList<FlatBottomBondRestraint> FlatBottomBondRestraints::centroidRestraints() const
{
    if (this->hasAtomRestraints())
    {
        QList<FlatBottomBondRestraint> cr;

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
void FlatBottomBondRestraints::add(const FlatBottomBondRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void FlatBottomBondRestraints::add(const FlatBottomBondRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
FlatBottomBondRestraints &FlatBottomBondRestraints::operator+=(const FlatBottomBondRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
FlatBottomBondRestraints FlatBottomBondRestraints::operator+(const FlatBottomBondRestraint &restraint) const
{
    FlatBottomBondRestraints ret(*this);
    ret += restraint;
    return ret;
}

/** Add restraints onto the list */
FlatBottomBondRestraints &FlatBottomBondRestraints::operator+=(const FlatBottomBondRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
FlatBottomBondRestraints FlatBottomBondRestraints::operator+(const FlatBottomBondRestraints &restraints) const
{
    FlatBottomBondRestraints ret(*this);
    ret += restraints;
    return ret;
}

/** Set whether or not periodic boundary conditions are to be used */
void FlatBottomBondRestraints::setUsesPbc(bool use_pbc)
{
    this->use_pbc = use_pbc;
}

/** Return whether or not periodic boundary conditions are to be used */
bool FlatBottomBondRestraints::usesPbc() const
{
    return this->use_pbc;
}
