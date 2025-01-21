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

#include "dihedralrestraints.h"

#include "SireID/index.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireUnits/units.h"

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
//////////// Implementation of DihedralRestraint
////////////

static const RegisterMetaType<DihedralRestraint> r_dihrest;

/** Serialise to a binary datastream */

QDataStream &operator<<(QDataStream &ds, const DihedralRestraint &dihrest)
{
    writeHeader(ds, r_dihrest, 1);

    SharedDataStream sds(ds);

    sds << dihrest.atms << dihrest._phi0 << dihrest._kphi;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, DihedralRestraint &dihrest)
{
    VersionID v = readHeader(ds, r_dihrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> dihrest.atms >> dihrest._phi0 >> dihrest._kphi;
    }
    else
        throw version_error(v, "1", r_dihrest, CODELOC);

    return ds;
}

/** Null constructor */
DihedralRestraint::DihedralRestraint()
    : ConcreteProperty<DihedralRestraint, Property>(),
      _phi0(0), _kphi(0)
{
}

/** Construct a restraint that acts on the angle within the
    four atoms 'atom0', 'atom1', 'atom2' 'atom3' (phi == a(0123)),
    restraining the angle within these atoms */
DihedralRestraint::DihedralRestraint(const QList<qint64> &atoms,
                                     const SireUnits::Dimension::Angle &phi0,
                                     const SireUnits::Dimension::HarmonicAngleConstant &kphi)
    : ConcreteProperty<DihedralRestraint, Property>(),
      _phi0(phi0), _kphi(kphi)
{

    // Make sure that we have 4 distinct atoms
    QSet<qint64> distinct;
    distinct.reserve(4);

    for (const auto &atom : atoms)
    {
        if (atom >= 0)
            distinct.insert(atom);
    }

    atms = atoms.toVector();
}

/* Copy constructor*/
DihedralRestraint::DihedralRestraint(const DihedralRestraint &other)
    : ConcreteProperty<DihedralRestraint, Property>(other),
      atms(other.atms), _phi0(other._phi0), _kphi(other._kphi)

{
}

/* Destructor */
DihedralRestraint::~DihedralRestraint()
{
}

DihedralRestraint &DihedralRestraint::operator=(const DihedralRestraint &other)
{
    if (this != &other)
    {
        Property::operator=(other);
        atms = other.atms;
        _phi0 = other._phi0;
        _kphi = other._kphi;
    }

    return *this;
}

bool DihedralRestraint::operator==(const DihedralRestraint &other) const
{
    return atms == other.atms and
           _phi0 == other._phi0 and
           _kphi == other._kphi;
}

bool DihedralRestraint::operator!=(const DihedralRestraint &other) const
{
    return not operator==(other);
}

DihedralRestraints DihedralRestraint::operator+(const DihedralRestraint &other) const
{
    return DihedralRestraints(*this) + other;
}

DihedralRestraints DihedralRestraint::operator+(const DihedralRestraints &other) const
{
    return DihedralRestraints(*this) + other;
}

const char *DihedralRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<DihedralRestraint>());
}

const char *DihedralRestraint::what() const
{
    return DihedralRestraint::typeName();
}

DihedralRestraint *DihedralRestraint::clone() const
{
    return new DihedralRestraint(*this);
}

bool DihedralRestraint::isNull() const
{
    return atms.isEmpty();
}

QString DihedralRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("DihedralRestraint::null");
    else
    {
        QStringList a;

        for (const auto &atom : atms)
        {
            a.append(QString::number(atom));
        }
        return QString("DihedralRestraint( [%1], phi0=%2, kphi=%3 )")
            .arg(a.join(", "))
            .arg(_phi0.toString())
            .arg(_kphi.toString());
    }
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicAngleConstant DihedralRestraint::kphi() const
{
    return this->_kphi;
}

/** Return the equilibrium angle for the restraint */
SireUnits::Dimension::Angle DihedralRestraint::phi0() const
{
    return this->_phi0;
}

/** Return the atoms involved in the restraint */
QVector<qint64> DihedralRestraint::atoms() const
{
    return this->atms;
}

///////
/////// Implementation of DihedralRestraints
///////

/** Serialise to a binary datastream */

static const RegisterMetaType<DihedralRestraints> r_dihrests;

QDataStream &operator<<(QDataStream &ds, const DihedralRestraints &dihrests)
{
    writeHeader(ds, r_dihrests, 1);

    SharedDataStream sds(ds);

    sds << dihrests.r
        << static_cast<const Restraints &>(dihrests);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, DihedralRestraints &dihrests)
{
    VersionID v = readHeader(ds, r_dihrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> dihrests.r >>
            static_cast<Restraints &>(dihrests);
    }
    else
        throw version_error(v, "1", r_dihrests, CODELOC);

    return ds;
}

/** Null constructor */
DihedralRestraints::DihedralRestraints()
    : ConcreteProperty<DihedralRestraints, Restraints>()
{
}

DihedralRestraints::DihedralRestraints(const QString &name)
    : ConcreteProperty<DihedralRestraints, Restraints>(name)
{
}

DihedralRestraints::DihedralRestraints(const DihedralRestraint &restraint)
    : ConcreteProperty<DihedralRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

DihedralRestraints::DihedralRestraints(const QList<DihedralRestraint> &restraints)
    : ConcreteProperty<DihedralRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

DihedralRestraints::DihedralRestraints(const QString &name,
                                       const DihedralRestraint &restraint)
    : ConcreteProperty<DihedralRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

DihedralRestraints::DihedralRestraints(const QString &name,
                                       const QList<DihedralRestraint> &restraints)
    : ConcreteProperty<DihedralRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

DihedralRestraints::DihedralRestraints(const DihedralRestraints &other)
    : ConcreteProperty<DihedralRestraints, Restraints>(other), r(other.r)
{
}

/* Desctructor */
DihedralRestraints::~DihedralRestraints()
{
}

DihedralRestraints &DihedralRestraints::operator=(const DihedralRestraints &other)
{
    if (this != &other)
    {
        Restraints::operator=(other);
        r = other.r;
    }

    return *this;
}

bool DihedralRestraints::operator==(const DihedralRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
}

bool DihedralRestraints::operator!=(const DihedralRestraints &other) const
{
    return not operator==(other);
}

const char *DihedralRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<DihedralRestraints>());
}

const char *DihedralRestraints::what() const
{
    return DihedralRestraints::typeName();
}

DihedralRestraints *DihedralRestraints::clone() const
{
    return new DihedralRestraints(*this);
}

QString DihedralRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("DihedralRestraints::null");

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

    return QObject::tr("DihedralRestraints( name=%1, size=%2\n%3\n )")
        .arg(this->name())
        .arg(n)
        .arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool DihedralRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool DihedralRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int DihedralRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int DihedralRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int DihedralRestraints::size() const
{
    return this->nRestraints();
}

/** Return the ith restraint */
const DihedralRestraint &DihedralRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const DihedralRestraint &DihedralRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<DihedralRestraint> DihedralRestraints::restraints() const
{
    return this->r;
}

/** Add a restraints onto the list */
void DihedralRestraints::add(const DihedralRestraint &restraint)
{
    if (not restraint.isNull())
        r.append(restraint);
}

/** Add a restraint onto the list */
void DihedralRestraints::add(const DihedralRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
DihedralRestraints &DihedralRestraints::operator+=(const DihedralRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
DihedralRestraints DihedralRestraints::operator+(const DihedralRestraint &restraint) const
{
    DihedralRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
DihedralRestraints &DihedralRestraints::operator+=(const DihedralRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
DihedralRestraints DihedralRestraints::operator+(const DihedralRestraints &restraints) const
{
    DihedralRestraints ret(*this);
    ret += restraints;
    return *this;
}