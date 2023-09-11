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

#include "boreschrestraints.h"

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
/////// Implementation of BoreschRestraint
///////

static const RegisterMetaType<BoreschRestraint> r_borrest;

QDataStream &operator<<(QDataStream &ds, const BoreschRestraint &borrest)
{
    writeHeader(ds, r_borrest, 1);

    SharedDataStream sds(ds);

    sds << borrest.receptor_atms << borrest.ligand_atms
        << borrest._r0 << borrest._theta0 << borrest._phi0
        << borrest._kr << borrest._ktheta << borrest._kphi
        << static_cast<const Property &>(borrest);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BoreschRestraint &borrest)
{
    VersionID v = readHeader(ds, r_borrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> borrest.receptor_atms >> borrest.ligand_atms >> borrest._r0 >> borrest._theta0 >> borrest._phi0 >> borrest._kr >> borrest._ktheta >> borrest._kphi >> static_cast<Property &>(borrest);
    }
    else
        throw version_error(v, "1", r_borrest, CODELOC);

    return ds;
}

/** Null constructor */
BoreschRestraint::BoreschRestraint()
    : ConcreteProperty<BoreschRestraint, Property>(),
      _r0(0), _kr(0)
{
}

/** Construct to restrain the atom at index 'atom' to the specified position
 *  using the specified force constant and flat-bottom well-width
 */
BoreschRestraint::BoreschRestraint(const QList<qint64> &receptor_atoms,
                                   const QList<qint64> &ligand_atoms,
                                   const SireUnits::Dimension::Length &r0,
                                   const QVector<SireUnits::Dimension::Angle> &theta0,
                                   const QVector<SireUnits::Dimension::Angle> &phi0,
                                   const SireUnits::Dimension::HarmonicBondConstant &kr,
                                   const QVector<SireUnits::Dimension::HarmonicAngleConstant> &ktheta,
                                   const QVector<SireUnits::Dimension::HarmonicAngleConstant> &kphi)
    : ConcreteProperty<BoreschRestraint, Property>(),
      _r0(r0), _theta0(theta0), _phi0(phi0),
      _kr(kr), _ktheta(ktheta), _kphi(kphi)
{
    if (receptor_atoms.count() != 3 or ligand_atoms.count() != 3 or
        theta0.count() != 2 or phi0.count() != 3 or
        ktheta.count() != 2 or kphi.count() != 3)
    {
        throw SireError::invalid_arg(QObject::tr(
                                         "Wrong number of inputs for a Boresch restraint. You need to "
                                         "provide 3 receptor atoms (%1), 3 ligand atoms (%2), 2 theta0 values (%3) "
                                         "3 phi0 values (%4), 2 ktheta values (%5) and 3 kphi values (%6).")
                                         .arg(receptor_atoms.count())
                                         .arg(ligand_atoms.count())
                                         .arg(theta0.count())
                                         .arg(phi0.count())
                                         .arg(ktheta.count())
                                         .arg(kphi.count()),
                                     CODELOC);
    }

    // make sure that we have 6 distinct atoms
    QSet<qint64> distinct;
    distinct.reserve(6);

    for (const auto &atom : receptor_atoms)
    {
        if (atom >= 0)
            distinct.insert(atom);
    }

    for (const auto &atom : ligand_atoms)
    {
        if (atom >= 0)
            distinct.insert(atom);
    }

    if (distinct.count() != 6)
        throw SireError::invalid_arg(QObject::tr(
                                         "There is something wrong with the receptor or ligand atoms. "
                                         "They should all be unique and all greater than or equal to 0."),
                                     CODELOC);

    receptor_atms = receptor_atoms.toVector();
    receptor_atms.squeeze();

    ligand_atms = ligand_atoms.toVector();
    ligand_atms.squeeze();
}

/** Copy constructor */
BoreschRestraint::BoreschRestraint(const BoreschRestraint &other)
    : ConcreteProperty<BoreschRestraint, Property>(other),
      receptor_atms(other.receptor_atms), ligand_atms(other.ligand_atms),
      _r0(other._r0), _theta0(other._theta0), _phi0(other._phi0),
      _kr(other._kr), _ktheta(other._ktheta), _kphi(other._kphi)
{
}

BoreschRestraint::~BoreschRestraint()
{
}

BoreschRestraint &BoreschRestraint::operator=(const BoreschRestraint &other)
{
    if (this != &other)
    {
        receptor_atms = other.receptor_atms;
        ligand_atms = other.ligand_atms;
        _r0 = other._r0;
        _theta0 = other._theta0;
        _phi0 = other._phi0;
        _kr = other._kr;
        _ktheta = other._ktheta;
        _kphi = other._kphi;
        Property::operator=(other);
    }

    return *this;
}

bool BoreschRestraint::operator==(const BoreschRestraint &other) const
{
    return receptor_atms == other.receptor_atms and
           ligand_atms == other.ligand_atms and
           _kr == other._kr and _ktheta == other._ktheta and
           _kphi == other._kphi and _r0 == other._r0 and
           _theta0 == other._theta0 and _phi0 == other._phi0;
}

bool BoreschRestraint::operator!=(const BoreschRestraint &other) const
{
    return not operator==(other);
}

BoreschRestraints BoreschRestraint::operator+(const BoreschRestraint &other) const
{
    return BoreschRestraints(*this) + other;
}

BoreschRestraints BoreschRestraint::operator+(const BoreschRestraints &other) const
{
    return BoreschRestraints(*this) + other;
}

const char *BoreschRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<BoreschRestraint>());
}

const char *BoreschRestraint::what() const
{
    return BoreschRestraint::typeName();
}

BoreschRestraint *BoreschRestraint::clone() const
{
    return new BoreschRestraint(*this);
}

bool BoreschRestraint::isNull() const
{
    return receptor_atms.isEmpty() and ligand_atms.isEmpty();
}

QString BoreschRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("BoreschRestraint::null");

    else
    {
        QStringList r;

        for (const auto &atom : receptor_atms)
        {
            r.append(QString::number(atom));
        }

        QStringList l;

        for (const auto &atom : ligand_atms)
        {
            l.append(QString::number(atom));
        }

        QStringList k;

        k.append(_kr.toString());

        for (const auto &kt : _ktheta)
        {
            k.append(kt.toString());
        }

        for (const auto &kp : _kphi)
        {
            k.append(kp.toString());
        }

        QStringList t;

        for (const auto &t0 : _theta0)
        {
            t.append(t0.toString());
        }

        QStringList p;

        for (const auto &p0 : _phi0)
        {
            p.append(p0.toString());
        }

        return QString("BoreschRestraint( [%1] => [%2],\n"
                       "                  k=[%3]\n"
                       "                  r0=%4, theta0=[%5],\n"
                       "                  phi0=[%6] )")
            .arg(r.join(", "))
            .arg(l.join(", "))
            .arg(k.join(", "))
            .arg(_r0.toString())
            .arg(t.join(", "))
            .arg(p.join(', '));
    }
}

/** Return the indexes of the three receptor atoms */
QVector<qint64> BoreschRestraint::receptorAtoms() const
{
    return this->receptor_atms;
}

/** Return the indexes of the three ligand atoms */
QVector<qint64> BoreschRestraint::ligandAtoms() const
{
    return this->ligand_atms;
}

/** Return the equilibrium length of the bond restraint */
SireUnits::Dimension::Length BoreschRestraint::r0() const
{
    return _r0;
}

/** Return the equilibrium size of the two angle restraints */
QVector<SireUnits::Dimension::Angle> BoreschRestraint::theta0() const
{
    return _theta0;
}

/** Return the equilibium size of the three dihedral restraints */
QVector<SireUnits::Dimension::Angle> BoreschRestraint::phi0() const
{
    return _phi0;
}

/** Return the force constant for the bond restraint */
SireUnits::Dimension::HarmonicBondConstant BoreschRestraint::kr() const
{
    return _kr;
}

/** Return the force constant for the two angle restraints */
QVector<SireUnits::Dimension::HarmonicAngleConstant> BoreschRestraint::ktheta() const
{
    return _ktheta;
}

/** Return the force constant for the three dihedral restraints */
QVector<SireUnits::Dimension::HarmonicAngleConstant> BoreschRestraint::kphi() const
{
    return _kphi;
}

///////
/////// Implementation of BoreschRestraints
///////

static const RegisterMetaType<BoreschRestraints> r_borrests;

QDataStream &operator<<(QDataStream &ds, const BoreschRestraints &borrests)
{
    writeHeader(ds, r_borrests, 1);

    SharedDataStream sds(ds);

    sds << borrests.r
        << static_cast<const Restraints &>(borrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BoreschRestraints &borrests)
{
    VersionID v = readHeader(ds, r_borrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> borrests.r >> static_cast<Restraints &>(borrests);
    }
    else
        throw version_error(v, "1", r_borrests, CODELOC);

    return ds;
}

/** Null constructor */
BoreschRestraints::BoreschRestraints()
    : ConcreteProperty<BoreschRestraints, Restraints>()
{
}

BoreschRestraints::BoreschRestraints(const QString &name)
    : ConcreteProperty<BoreschRestraints, Restraints>(name)
{
}

BoreschRestraints::BoreschRestraints(const BoreschRestraint &restraint)
    : ConcreteProperty<BoreschRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

BoreschRestraints::BoreschRestraints(const QList<BoreschRestraint> &restraints)
    : ConcreteProperty<BoreschRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

BoreschRestraints::BoreschRestraints(const QString &name,
                                     const BoreschRestraint &restraint)
    : ConcreteProperty<BoreschRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

BoreschRestraints::BoreschRestraints(const QString &name,
                                     const QList<BoreschRestraint> &restraints)
    : ConcreteProperty<BoreschRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

BoreschRestraints::BoreschRestraints(const BoreschRestraints &other)
    : ConcreteProperty<BoreschRestraints, Restraints>(other), r(other.r)
{
}

BoreschRestraints::~BoreschRestraints()
{
}

BoreschRestraints &BoreschRestraints::operator=(const BoreschRestraints &other)
{
    r = other.r;
    Restraints::operator=(other);
    return *this;
}

bool BoreschRestraints::operator==(const BoreschRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
}

bool BoreschRestraints::operator!=(const BoreschRestraints &other) const
{
    return not operator==(other);
}

const char *BoreschRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<BoreschRestraints>());
}

const char *BoreschRestraints::what() const
{
    return BoreschRestraints::typeName();
}

BoreschRestraints *BoreschRestraints::clone() const
{
    return new BoreschRestraints(*this);
}

QString BoreschRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("BoreschRestraints::null");

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

    return QObject::tr("BoreschRestraints( name=%1, size=%2\n%3\n)").arg(this->name()).arg(n).arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool BoreschRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool BoreschRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int BoreschRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int BoreschRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int BoreschRestraints::size() const
{
    return this->nRestraints();
}

/** Return the ith restraint */
const BoreschRestraint &BoreschRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const BoreschRestraint &BoreschRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<BoreschRestraint> BoreschRestraints::restraints() const
{
    return this->r;
}

/** Add a restraint onto the list */
void BoreschRestraints::add(const BoreschRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void BoreschRestraints::add(const BoreschRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
BoreschRestraints &BoreschRestraints::operator+=(const BoreschRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
BoreschRestraints BoreschRestraints::operator+(const BoreschRestraint &restraint) const
{
    BoreschRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
BoreschRestraints &BoreschRestraints::operator+=(const BoreschRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
BoreschRestraints BoreschRestraints::operator+(const BoreschRestraints &restraints) const
{
    BoreschRestraints ret(*this);
    ret += restraints;
    return *this;
}
