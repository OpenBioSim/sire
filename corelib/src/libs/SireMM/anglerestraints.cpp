/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2025  Christopher Woods
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

#include "anglerestraints.h"

#include "SireID/index.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireCAS/errors.h"

#include <QDebug>

using namespace SireMM;
using namespace SireID;
using namespace SireBase;
using namespace SireMaths;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;

////////////
//////////// Implementation of AngleRestraint
////////////

static const RegisterMetaType<AngleRestraint> r_angrest;

/** Serialise to a binary datastream */

QDataStream &operator<<(QDataStream &ds, const AngleRestraint &angrest)
{
    writeHeader(ds, r_angrest, 1);

    SharedDataStream sds(ds);

    sds << angrest.atms << angrest._theta0 << angrest._ktheta;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AngleRestraint &angrest)
{
    VersionID v = readHeader(ds, r_angrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> angrest.atms >> angrest._theta0 >> angrest._ktheta;
    }
    else
        throw version_error(v, "1", r_angrest, CODELOC);

    return ds;
}

/** Null constructor */
AngleRestraint::AngleRestraint()
    : ConcreteProperty<AngleRestraint, Property>(),
      _ktheta(0), _theta0(0)
{
}

/** Construct a restraint that acts on the angle within the
    three atoms 'atom0', 'atom1' and 'atom2' (theta == a(012)),
    restraining the angle within these atoms */
AngleRestraint::AngleRestraint(const QList<qint64> &atoms,
                               const SireUnits::Dimension::Angle &theta0,
                               const SireUnits::Dimension::HarmonicAngleConstant &ktheta)
    : ConcreteProperty<AngleRestraint, Property>(),
      _theta0(theta0), _ktheta(ktheta)
{

    // Make sure that we have 3 distinct atoms
    QSet<qint64> distinct;
    distinct.reserve(3);

    for (const auto &atom : atoms)
    {
        if (atom >= 0)
            distinct.insert(atom);
    }

    atms = atoms.toVector();
}

/* Copy constructor*/
AngleRestraint::AngleRestraint(const AngleRestraint &other)
    : ConcreteProperty<AngleRestraint, Property>(other),
      atms(other.atms), _theta0(other._theta0), _ktheta(other._ktheta)

{
}

/* Destructor */
AngleRestraint::~AngleRestraint()
{
}

AngleRestraint &AngleRestraint::operator=(const AngleRestraint &other)
{
    if (this != &other)
    {
        Property::operator=(other);
        atms = other.atms;
        _theta0 = other._theta0;
        _ktheta = other._ktheta;
    }

    return *this;
}

bool AngleRestraint::operator==(const AngleRestraint &other) const
{
    return atms == other.atms and
           _theta0 == other._theta0 and
           _ktheta == other._ktheta;
}

bool AngleRestraint::operator!=(const AngleRestraint &other) const
{
    return not operator==(other);
}

AngleRestraints AngleRestraint::operator+(const AngleRestraint &other) const
{
    return AngleRestraints(*this) + other;
}

AngleRestraints AngleRestraint::operator+(const AngleRestraints &other) const
{
    return AngleRestraints(*this) + other;
}

const char *AngleRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AngleRestraint>());
}

const char *AngleRestraint::what() const
{
    return AngleRestraint::typeName();
}

AngleRestraint *AngleRestraint::clone() const
{
    return new AngleRestraint(*this);
}

bool AngleRestraint::isNull() const
{
    return atms.isEmpty();
}

QString AngleRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("AngleRestraint::null");
    else
    {
        QStringList a;

        for (const auto &atom : atms)
        {
            a.append(QString::number(atom));
        }
        return QString("AngleRestraint( [%1], theta0=%2, ktheta=%3 )")
            .arg(a.join(", "))
            .arg(_theta0.toString())
            .arg(_ktheta.toString());
    }
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicAngleConstant AngleRestraint::ktheta() const
{
    return this->_ktheta;
}

/** Return the equilibrium angle for the restraint */
SireUnits::Dimension::Angle AngleRestraint::theta0() const
{
    return this->_theta0;
}

/** Return the atoms involved in the restraint */
QVector<qint64> AngleRestraint::atoms() const
{
    return this->atms;
}

///////
/////// Implementation of AngleRestraints
///////

/** Serialise to a binary datastream */

static const RegisterMetaType<AngleRestraints> r_angrests;

QDataStream &operator<<(QDataStream &ds, const AngleRestraints &angrests)
{
    writeHeader(ds, r_angrests, 2);

    SharedDataStream sds(ds);

    sds << angrests.r << angrests.use_pbc
        << static_cast<const Restraints &>(angrests);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AngleRestraints &angrests)
{
    VersionID v = readHeader(ds, r_angrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> angrests.r >>
            static_cast<Restraints &>(angrests);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> angrests.r >> angrests.use_pbc
            >> static_cast<Restraints &>(angrests);
    }
    else
        throw version_error(v, "1,2", r_angrests, CODELOC);

    return ds;
}

/** Null constructor */
AngleRestraints::AngleRestraints()
    : ConcreteProperty<AngleRestraints, Restraints>()
{
}

AngleRestraints::AngleRestraints(const QString &name)
    : ConcreteProperty<AngleRestraints, Restraints>(name)
{
}

AngleRestraints::AngleRestraints(const AngleRestraint &restraint)
    : ConcreteProperty<AngleRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

AngleRestraints::AngleRestraints(const QList<AngleRestraint> &restraints)
    : ConcreteProperty<AngleRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

AngleRestraints::AngleRestraints(const QString &name,
                                 const AngleRestraint &restraint)
    : ConcreteProperty<AngleRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

AngleRestraints::AngleRestraints(const QString &name,
                                 const QList<AngleRestraint> &restraints)
    : ConcreteProperty<AngleRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

AngleRestraints::AngleRestraints(const AngleRestraints &other)
    : ConcreteProperty<AngleRestraints, Restraints>(other), r(other.r), use_pbc(other.use_pbc)
{
}

/* Desctructor */
AngleRestraints::~AngleRestraints()
{
}

AngleRestraints &AngleRestraints::operator=(const AngleRestraints &other)
{
    if (this != &other)
    {
        Restraints::operator=(other);
        r = other.r;
        use_pbc = other.use_pbc;
    }

    return *this;
}

bool AngleRestraints::operator==(const AngleRestraints &other) const
{
    return r == other.r and Restraints::operator==(other) and use_pbc == other.use_pbc;
}

bool AngleRestraints::operator!=(const AngleRestraints &other) const
{
    return not operator==(other);
}

const char *AngleRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AngleRestraints>());
}

const char *AngleRestraints::what() const
{
    return AngleRestraints::typeName();
}

AngleRestraints *AngleRestraints::clone() const
{
    return new AngleRestraints(*this);
}

QString AngleRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("AngleRestraints::null");

    QStringList parts;

    const auto n = this->count();

    if (n <= 10)
    {
        for (int i = 0; i < n; i++)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }
    else
    {
        for (int i = 0; i < 5; i++)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }

        parts.append("...");

        for (int i = n - 5; i < n; i++)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }

    return QObject::tr("AngleRestraints( name=%1, size=%2, use_pbc=%3\n%4\n )")
        .arg(this->name())
        .arg(n)
        .arg(this->use_pbc ? "true" : "false")
        .arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool AngleRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool AngleRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int AngleRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int AngleRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int AngleRestraints::size() const
{
    return this->nRestraints();
}

/** Return the ith restraint */
const AngleRestraint &AngleRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const AngleRestraint &AngleRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<AngleRestraint> AngleRestraints::restraints() const
{
    return this->r;
}

/** Add a restraints onto the list */
void AngleRestraints::add(const AngleRestraint &restraint)
{
    if (not restraint.isNull())
        r.append(restraint);
}

/** Add a restraint onto the list */
void AngleRestraints::add(const AngleRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
AngleRestraints &AngleRestraints::operator+=(const AngleRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
AngleRestraints AngleRestraints::operator+(const AngleRestraint &restraint) const
{
    AngleRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
AngleRestraints &AngleRestraints::operator+=(const AngleRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
AngleRestraints AngleRestraints::operator+(const AngleRestraints &restraints) const
{
    AngleRestraints ret(*this);
    ret += restraints;
    return *this;
}

/** Set whether or not periodic boundary conditions are to be used */
void AngleRestraints::setUsesPbc(bool use_pbc)
{
    this->use_pbc = use_pbc;
}

/** Return whether or not periodic boundary conditions are to be used */
bool AngleRestraints::usesPbc() const
{
    return this->use_pbc;
}
