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

#include "morsepotentialrestraints.h"

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
/////// Implementation of MorsePotentialRestraint
///////

static const RegisterMetaType<MorsePotentialRestraint> r_morsepotentialrest;

QDataStream &operator<<(QDataStream &ds, const MorsePotentialRestraint &morsepotentialrest)
{
    writeHeader(ds, r_morsepotentialrest, 1);

    SharedDataStream sds(ds);

    sds << morsepotentialrest.atms0 << morsepotentialrest.atms1 << morsepotentialrest._k << morsepotentialrest._r0 << morsepotentialrest._de
        << static_cast<const Property &>(morsepotentialrest);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MorsePotentialRestraint &morsepotentialrest)
{
    VersionID v = readHeader(ds, r_morsepotentialrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> morsepotentialrest.atms0 >> morsepotentialrest.atms1 >> morsepotentialrest._k >> morsepotentialrest._r0 >> morsepotentialrest._de >> static_cast<Property &>(morsepotentialrest);
    }
    else
        throw version_error(v, "1", r_morsepotentialrest, CODELOC);

    return ds;
}

/** Null constructor */
MorsePotentialRestraint::MorsePotentialRestraint()
    : ConcreteProperty<MorsePotentialRestraint, Property>(),
      _k(0), _r0(0), _de(0)
{
}

/** Construct to restrain the atom at index 'atom' to the specified position
 *  using the specified force constant and flat-bottom well-width
 */
MorsePotentialRestraint::MorsePotentialRestraint(qint64 atom0, qint64 atom1,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0,
                             const SireUnits::Dimension::MolarEnergy &de)
    : ConcreteProperty<MorsePotentialRestraint, Property>(),
      _k(k), _r0(r0), _de(de)
{
    if (atom0 == atom1)
        throw SireError::invalid_arg(QObject::tr(
                                         "You cannot create a MorsePotential restraint between identical atoms! %1-%2")
                                         .arg(atom0)
                                         .arg(atom1),
                                     CODELOC);

    atms0 = QVector<qint64>(1, atom0);
    atms0.squeeze();

    atms1 = QVector<qint64>(1, atom1);
    atms1.squeeze();
}

/** Construct to restrain the atoms whose indicies are
 *  in 'atoms' to the specified position using the specified force constant
 *  and dissociation energy
 */
MorsePotentialRestraint::MorsePotentialRestraint(const QList<qint64> &atoms0,
    const QList<qint64> &atoms1,
    const SireUnits::Dimension::HarmonicBondConstant &k,
    const SireUnits::Dimension::Length &r0,
    const SireUnits::Dimension::MolarEnergy &de)
: ConcreteProperty<MorsePotentialRestraint, Property>(),
_k(k), _r0(r0), _de(de)
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
MorsePotentialRestraint::MorsePotentialRestraint(const MorsePotentialRestraint &other)
    : ConcreteProperty<MorsePotentialRestraint, Property>(other),
      atms0(other.atms0), atms1(other.atms1), _k(other._k), _r0(other._r0), _de(other._de)
{
}

MorsePotentialRestraint::~MorsePotentialRestraint()
{
}

MorsePotentialRestraint &MorsePotentialRestraint::operator=(const MorsePotentialRestraint &other)
{
    if (this != &other)
    {
        atms0 = other.atms0;
        atms1 = other.atms1;
        _k = other._k;
        _r0 = other._r0;
        _de = other._de;
    }

    return *this;
}

bool MorsePotentialRestraint::operator==(const MorsePotentialRestraint &other) const
{
    return atms0 == other.atms0 and atms1 == other.atms1 and
           _k == other._k and _r0 == other._r0 and _de == other._de;
}

bool MorsePotentialRestraint::operator!=(const MorsePotentialRestraint &other) const
{
    return not operator==(other);
}

MorsePotentialRestraints MorsePotentialRestraint::operator+(const MorsePotentialRestraint &other) const
{
    return MorsePotentialRestraints(*this) + other;
}

MorsePotentialRestraints MorsePotentialRestraint::operator+(const MorsePotentialRestraints &other) const
{
    return MorsePotentialRestraints(*this) + other;
}

const char *MorsePotentialRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<MorsePotentialRestraint>());
}

const char *MorsePotentialRestraint::what() const
{
    return MorsePotentialRestraint::typeName();
}

MorsePotentialRestraint *MorsePotentialRestraint::clone() const
{
    return new MorsePotentialRestraint(*this);
}

bool MorsePotentialRestraint::isNull() const
{
    return atms0.isEmpty() or atms1.isEmpty();
}

QString MorsePotentialRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("MorsePotentialRestraint::null");

    else if (this->isAtomRestraint())
    {
        return QString("MorsePotentialRestraint( %1 <=> %2, k=%3 : r0=%4 : de=%5 )")
            .arg(this->atom0())
            .arg(this->atom1())
            .arg(_k.toString())
            .arg(_r0.toString())
            .arg(_de.toString());
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

        return QString("MorsePotentialRestraint( [%1] <=> [%2], k=%3 : r0=%4 : de=%5 )")
            .arg(a0.join(", "))
            .arg(a1.join(", "))
            .arg(_k.toString())
            .arg(_r0.toString())
            .arg(_de.toString());
    }
}

/** Return whether this is a single-atom restraint */
bool MorsePotentialRestraint::isAtomRestraint() const
{
    return atms0.count() == 1 and atms1.count() == 1;
}

/** Return whether this restraint acts on the centroid of a group
 *  of atoms */
bool MorsePotentialRestraint::isCentroidRestraint() const
{
    return atms0.count() > 1 or atms1.count() > 1;
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 MorsePotentialRestraint::atom0() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms0.at(0);
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 MorsePotentialRestraint::atom1() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms1.at(0);
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> MorsePotentialRestraint::atoms0() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms0;
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> MorsePotentialRestraint::atoms1() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms1;
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicBondConstant MorsePotentialRestraint::k() const
{
    return this->_k;
}

/** Return the equilibrium position of the harmonic MovingHarmonic. */
SireUnits::Dimension::Length MorsePotentialRestraint::r0() const
{
    return this->_r0;
}

/** Return the dissociation energy (depth) of MovingHarmonic bond */
// SireUnits::Dimension::MolarEnergy MorsePotentialRestraint::r1() const
SireUnits::Dimension::MolarEnergy MorsePotentialRestraint::de() const
{
    return this->_de;
}

///////
/////// Implementation of MorsePotentialRestraints
///////

static const RegisterMetaType<MorsePotentialRestraints> r_morsepotentialrests;

QDataStream &operator<<(QDataStream &ds, const MorsePotentialRestraints &morser_morsepotentialrests)
{
    writeHeader(ds, r_morsepotentialrests, 1);

    SharedDataStream sds(ds);

    sds << morser_morsepotentialrests.r
        << static_cast<const Restraints &>(morser_morsepotentialrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MorsePotentialRestraints &morser_morsepotentialrests)
{
    VersionID v = readHeader(ds, r_morsepotentialrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> morser_morsepotentialrests.r >> static_cast<Restraints &>(morser_morsepotentialrests);
    }
    else
        throw version_error(v, "1", r_morsepotentialrests, CODELOC);

    return ds;
}

/** Null constructor */
MorsePotentialRestraints::MorsePotentialRestraints()
    : ConcreteProperty<MorsePotentialRestraints, Restraints>()
{
}

MorsePotentialRestraints::MorsePotentialRestraints(const QString &name)
    : ConcreteProperty<MorsePotentialRestraints, Restraints>(name)
{
}

MorsePotentialRestraints::MorsePotentialRestraints(const MorsePotentialRestraint &restraint)
    : ConcreteProperty<MorsePotentialRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

MorsePotentialRestraints::MorsePotentialRestraints(const QList<MorsePotentialRestraint> &restraints)
    : ConcreteProperty<MorsePotentialRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

MorsePotentialRestraints::MorsePotentialRestraints(const QString &name,
                               const MorsePotentialRestraint &restraint)
    : ConcreteProperty<MorsePotentialRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

MorsePotentialRestraints::MorsePotentialRestraints(const QString &name,
                               const QList<MorsePotentialRestraint> &restraints)
    : ConcreteProperty<MorsePotentialRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

MorsePotentialRestraints::MorsePotentialRestraints(const MorsePotentialRestraints &other)
    : ConcreteProperty<MorsePotentialRestraints, Restraints>(other), r(other.r)
{
}

MorsePotentialRestraints::~MorsePotentialRestraints()
{
}

MorsePotentialRestraints &MorsePotentialRestraints::operator=(const MorsePotentialRestraints &other)
{
    r = other.r;
    Restraints::operator=(other);
    return *this;
}

bool MorsePotentialRestraints::operator==(const MorsePotentialRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
}

bool MorsePotentialRestraints::operator!=(const MorsePotentialRestraints &other) const
{
    return not operator==(other);
}

const char *MorsePotentialRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<MorsePotentialRestraints>());
}

const char *MorsePotentialRestraints::what() const
{
    return MorsePotentialRestraints::typeName();
}

MorsePotentialRestraints *MorsePotentialRestraints::clone() const
{
    return new MorsePotentialRestraints(*this);
}

QString MorsePotentialRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("MorsePotentialRestraints::null");

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

    return QObject::tr("MorsePotentialRestraints( name=%1, size=%2\n%3\n)").arg(this->name()).arg(n).arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool MorsePotentialRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool MorsePotentialRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int MorsePotentialRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int MorsePotentialRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int MorsePotentialRestraints::size() const
{
    return this->nRestraints();
}

/** Return the number of atom restraints */
int MorsePotentialRestraints::nAtomRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isAtomRestraint());
    }

    return n;
}

/** Return the number of centroid restraints */
int MorsePotentialRestraints::nCentroidRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isCentroidRestraint());
    }

    return n;
}

/** Return whether or not there are any atom restraints */
bool MorsePotentialRestraints::hasAtomRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isAtomRestraint())
            return true;
    }

    return false;
}

/** Return whether or not there are any centroid restraints */
bool MorsePotentialRestraints::hasCentroidRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isCentroidRestraint())
            return true;
    }

    return false;
}

/** Return the ith restraint */
const MorsePotentialRestraint &MorsePotentialRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const MorsePotentialRestraint &MorsePotentialRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<MorsePotentialRestraint> MorsePotentialRestraints::restraints() const
{
    return this->r;
}

/** Return all of the atom restraints */
QList<MorsePotentialRestraint> MorsePotentialRestraints::atomRestraints() const
{
    if (this->hasCentroidRestraints())
    {
        QList<MorsePotentialRestraint> ar;

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
QList<MorsePotentialRestraint> MorsePotentialRestraints::centroidRestraints() const
{
    if (this->hasAtomRestraints())
    {
        QList<MorsePotentialRestraint> cr;

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
void MorsePotentialRestraints::add(const MorsePotentialRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void MorsePotentialRestraints::add(const MorsePotentialRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
MorsePotentialRestraints &MorsePotentialRestraints::operator+=(const MorsePotentialRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
MorsePotentialRestraints MorsePotentialRestraints::operator+(const MorsePotentialRestraint &restraint) const
{
    MorsePotentialRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
MorsePotentialRestraints &MorsePotentialRestraints::operator+=(const MorsePotentialRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
MorsePotentialRestraints MorsePotentialRestraints::operator+(const MorsePotentialRestraints &restraints) const
{
    MorsePotentialRestraints ret(*this);
    ret += restraints;
    return *this;
}