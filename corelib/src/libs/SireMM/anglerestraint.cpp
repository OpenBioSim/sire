/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "anglerestraint.h"

// #include "SireFF/forcetable.h"

// #include "SireCAS/conditional.h"
// #include "SireCAS/power.h"
// #include "SireCAS/symbols.h"
// #include "SireCAS/values.h"

#include "SireID/index.h"

// #include "SireUnits/angle.h"
#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireCAS/errors.h"

#include <QDebug>

using namespace SireMM;
// using namespace SireMol;
// using namespace SireFF;
using namespace SireID;
using namespace SireBase;
// using namespace SireCAS;
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
    // Need to think here about validating the angle and force constant values
    // if (atoms.count() != 3)
    // {
    //     throw SireError::invalid_arg(QObject::tr(
    //                                      "Wrong number of inputs for an Angle restraint. You need to "
    //                                      "provide 3 atoms (%1).")
    //                                      .arg(atoms.count()),
    //                                  //  .arg(theta0.count())
    //                                  //  .arg(ktheta.count()),
    //                                  CODELOC);
    // }

    // make sure that we have 3 distinct atoms
    QSet<qint64> distinct;
    distinct.reserve(3);

    for (const auto &atom : atoms)
    {
        if (atom >= 0)
            distinct.insert(atom);
    }

    // if (distinct.count() != 3)
    //     throw SireError::invalid_arg(QObject::tr(
    //                                      "There is something wrong with the atoms provided. "
    //                                      "They should all be unique and all greater than or equal to 0."),
    //                                  CODELOC);

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
    writeHeader(ds, r_angrests, 1);

    SharedDataStream sds(ds);

    sds << angrests.r
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
    else
        throw version_error(v, "1", r_angrests, CODELOC);

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
    : ConcreteProperty<AngleRestraints, Restraints>(other), r(other.r)
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
    }

    return *this;
}

bool AngleRestraints::operator==(const AngleRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
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

    return QObject::tr("AngleRestraints( name=%1, size=%2\n%3\n )")
        .arg(this->name())
        .arg(n)
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

// QDataStream &operator<<(QDataStream &ds, const AngleRestraint &angrest)
// {
//     writeHeader(ds, r_angrest, 1);

//     SharedDataStream sds(ds);

//     sds << angrest.p[0] << angrest.p[1] << angrest.p[2] << angrest.force_expression
//         << static_cast<const ExpressionRestraint3D &>(angrest);

//     return ds;
// }

// /** Extract from a binary datastream */
// QDataStream &operator>>(QDataStream &ds, AngleRestraint &angrest)
// {
//     VersionID v = readHeader(ds, r_angrest);

//     if (v == 1)
//     {
//         SharedDataStream sds(ds);

//         sds >> angrest.p[0] >> angrest.p[1] >> angrest.p[2] >> angrest.force_expression >>
//             static_cast<ExpressionRestraint3D &>(angrest);

//         angrest.intra_molecule_points = Point::areIntraMoleculePoints(angrest.p[0], angrest.p[1]) and
//                                         Point::areIntraMoleculePoints(angrest.p[0], angrest.p[2]);
//     }
//     else
//         throw version_error(v, "1", r_angrest, CODELOC);

//     return ds;
// }

// Q_GLOBAL_STATIC_WITH_ARGS(Symbol, getThetaSymbol, ("theta"));

/** Return the symbol that represents the angle between the points (theta) */
// const Symbol &AngleRestraint::theta()
// {
//     return *(getThetaSymbol());
// }

// /** Constructor */
// AngleRestraint::AngleRestraint() : ConcreteProperty<AngleRestraint, ExpressionRestraint3D>()
// {
// }

// void AngleRestraint::calculateTheta()
// {
//     if (this->restraintFunction().isFunction(theta()))
//     {
//         SireUnits::Dimension::Angle angle;

//         if (intra_molecule_points)
//             // we don't use the space when calculating intra-molecular angles
//             angle = Vector::angle(p[0].read().point(), p[1].read().point(), p[2].read().point());
//         else
//             angle = this->space().calcAngle(p[0].read().point(), p[1].read().point(), p[2].read().point());

//         ExpressionRestraint3D::_pvt_setValue(theta(), angle);
//     }
// }

// AngleRestraint::AngleRestraint(const PointRef &point0, const PointRef &point1, const PointRef &point2,
//                                const Expression &restraint)
//     : ConcreteProperty<AngleRestraint, ExpressionRestraint3D>(restraint)
// {
//     p[0] = point0;
//     p[1] = point1;
//     p[2] = point2;

//     force_expression = this->restraintFunction().differentiate(theta());

//     if (force_expression.isConstant())
//         force_expression = force_expression.evaluate(Values());

//     intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]) and Point::areIntraMoleculePoints(p[0], p[2]);

//     this->calculateTheta();
// }

// /** Construct a restraint that acts on the angle within the
//     three points 'point0', 'point1' and 'point2' (theta == a(012)),
//     restraining the angle within these points using the expression
//     'restraint' */
// AngleRestraint::AngleRestraint(const PointRef &point0, const PointRef &point1, const PointRef &point2,
//                                const Expression &restraint, const Values &values)
//     : ConcreteProperty<AngleRestraint, ExpressionRestraint3D>(restraint, values)
// {
//     p[0] = point0;
//     p[1] = point1;
//     p[2] = point2;

//     force_expression = this->restraintFunction().differentiate(theta());

//     if (force_expression.isConstant())
//         force_expression = force_expression.evaluate(Values());

//     intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]) and Point::areIntraMoleculePoints(p[0], p[2]);

//     this->calculateTheta();
// }

// /** Internal constructor used to construct a restraint using the specified
//     points, energy expression and force expression */
// AngleRestraint::AngleRestraint(const PointRef &point0, const PointRef &point1, const PointRef &point2,
//                                const Expression &nrg_restraint, const Expression &force_restraint)
//     : ConcreteProperty<AngleRestraint, ExpressionRestraint3D>(nrg_restraint), force_expression(force_restraint)
// {
//     p[0] = point0;
//     p[1] = point1;
//     p[2] = point2;

//     if (force_expression.isConstant())
//     {
//         force_expression = force_expression.evaluate(Values());
//     }
//     else
//     {
//         if (not this->restraintFunction().symbols().contains(force_expression.symbols()))
//             throw SireError::incompatible_error(
//                 QObject::tr("You cannot use a force function which uses more symbols "
//                             "(%1) than the energy function (%2).")
//                     .arg(Sire::toString(force_expression.symbols()), Sire::toString(restraintFunction().symbols())),
//                 CODELOC);
//     }

//     intra_molecule_points = Point::areIntraMoleculePoints(p[0], p[1]) and Point::areIntraMoleculePoints(p[0], p[2]);

//     this->calculateTheta();
// }

// /** Copy constructor */
// AngleRestraint::AngleRestraint(const AngleRestraint &other)
//     : ConcreteProperty<AngleRestraint, ExpressionRestraint3D>(other), force_expression(other.force_expression),
//       intra_molecule_points(other.intra_molecule_points)
// {
//     for (int i = 0; i < this->nPoints(); ++i)
//     {
//         p[i] = other.p[i];
//     }
// }

// /** Destructor */
// AngleRestraint::~AngleRestraint()
// {
// }

// /** Copy assignment operator */
// AngleRestraint &AngleRestraint::operator=(const AngleRestraint &other)
// {
//     if (this != &other)
//     {
//         ExpressionRestraint3D::operator=(other);

//         for (int i = 0; i < this->nPoints(); ++i)
//         {
//             p[i] = other.p[i];
//         }

//         force_expression = other.force_expression;
//         intra_molecule_points = other.intra_molecule_points;
//     }

//     return *this;
// }

// /** Comparison operator */
// bool AngleRestraint::operator==(const AngleRestraint &other) const
// {
//     return this == &other or (ExpressionRestraint3D::operator==(other) and p[0] == other.p[0] and p[1] == other.p[1] and
//                               p[2] == other.p[2] and force_expression == other.force_expression);
// }

// /** Comparison operator */
// bool AngleRestraint::operator!=(const AngleRestraint &other) const
// {
//     return not AngleRestraint::operator==(other);
// }

// const char *AngleRestraint::typeName()
// {
//     return QMetaType::typeName(qMetaTypeId<AngleRestraint>());
// }

// /** This restraint involves three points */
// int AngleRestraint::nPoints() const
// {
//     return 3;
// }

// /** Return the ith point */
// const Point &AngleRestraint::point(int i) const
// {
//     i = Index(i).map(this->nPoints());

//     return p[i].read();
// }

// /** Return the first point */
// const Point &AngleRestraint::point0() const
// {
//     return p[0].read();
// }

// /** Return the second point */
// const Point &AngleRestraint::point1() const
// {
//     return p[1].read();
// }

// /** Return the third point */
// const Point &AngleRestraint::point2() const
// {
//     return p[2].read();
// }

// /** Return the built-in symbols of this restraint */
// Symbols AngleRestraint::builtinSymbols() const
// {
//     if (this->restraintFunction().isFunction(theta()))
//         return theta();
//     else
//         return Symbols();
// }

// /** Return the values of the built-in symbols of this restraint */
// Values AngleRestraint::builtinValues() const
// {
//     if (this->restraintFunction().isFunction(theta()))
//         return theta() == this->values()[theta()];
//     else
//         return Values();
// }

// /** Return the differential of this restraint with respect to
//     the symbol 'symbol'

//     \throw SireCAS::unavailable_differential
// */
// RestraintPtr AngleRestraint::differentiate(const Symbol &symbol) const
// {
//     if (this->restraintFunction().isFunction(symbol))
//         return AngleRestraint(p[0], p[1], p[2], restraintFunction().differentiate(symbol), this->values());
//     else
//         return NullRestraint();
// }

// /** Set the space used to evaluate the energy of this restraint

//     \throw SireVol::incompatible_space
// */
// void AngleRestraint::setSpace(const Space &new_space)
// {
//     if (not this->space().equals(new_space))
//     {
//         AngleRestraint old_state(*this);

//         try
//         {
//             for (int i = 0; i < this->nPoints(); ++i)
//             {
//                 p[i].edit().setSpace(new_space);
//             }

//             Restraint::setSpace(new_space);

//             this->calculateTheta();
//         }
//         catch (...)
//         {
//             AngleRestraint::operator=(old_state);
//             throw;
//         }
//     }
// }

// /** Return the function used to calculate the restraint force */
// const Expression &AngleRestraint::differentialRestraintFunction() const
// {
//     return force_expression;
// }

// /** Calculate the force acting on the molecule in the forcetable 'forcetable'
//     caused by this restraint, and add it on to the forcetable scaled by
//     'scale_force' */
// void AngleRestraint::force(MolForceTable &forcetable, double scale_force) const
// {
//     bool in_p0 = p[0].read().contains(forcetable.molNum());
//     bool in_p1 = p[1].read().contains(forcetable.molNum());
//     bool in_p2 = p[2].read().contains(forcetable.molNum());

//     if (not(in_p0 or in_p1 or in_p2))
//         // this molecule is not affected by the restraint
//         return;

//     throw SireError::incomplete_code(QObject::tr("Haven't yet written the code to calculate forces caused "
//                                                  "by an angle restraint."),
//                                      CODELOC);
// }

// /** Calculate the force acting on the molecules in the forcetable 'forcetable'
//     caused by this restraint, and add it on to the forcetable scaled by
//     'scale_force' */
// void AngleRestraint::force(ForceTable &forcetable, double scale_force) const
// {
//     bool in_p0 = p[0].read().usesMoleculesIn(forcetable);
//     bool in_p1 = p[1].read().usesMoleculesIn(forcetable);
//     bool in_p2 = p[2].read().usesMoleculesIn(forcetable);

//     if (not(in_p0 or in_p1 or in_p2))
//         // this molecule is not affected by the restraint
//         return;

//     throw SireError::incomplete_code(QObject::tr("Haven't yet written the code to calculate forces caused "
//                                                  "by an angle restraint."),
//                                      CODELOC);
// }

// /** Update the points of this restraint using new molecule data from 'moldata'

//     \throw SireBase::missing_property
//     \throw SireError::invalid_cast
//     \throw SireError::incompatible_error
// */
// void AngleRestraint::update(const MoleculeData &moldata)
// {
//     if (this->contains(moldata.number()))
//     {
//         AngleRestraint old_state(*this);

//         try
//         {
//             for (int i = 0; i < this->nPoints(); ++i)
//             {
//                 p[i].edit().update(moldata);
//             }

//             this->calculateTheta();
//         }
//         catch (...)
//         {
//             AngleRestraint::operator=(old_state);
//             throw;
//         }
//     }
// }

// /** Update the points of this restraint using new molecule data from 'molecules'

//     \throw SireBase::missing_property
//     \throw SireError::invalid_cast
//     \throw SireError::incompatible_error
// */
// void AngleRestraint::update(const Molecules &molecules)
// {
//     if (this->usesMoleculesIn(molecules))
//     {
//         AngleRestraint old_state(*this);

//         try
//         {
//             for (int i = 0; i < this->nPoints(); ++i)
//             {
//                 p[i].edit().update(molecules);
//             }

//             this->calculateTheta();
//         }
//         catch (...)
//         {
//             AngleRestraint::operator=(old_state);
//             throw;
//         }
//     }
// }

// /** Return the molecules used in this restraint */
// Molecules AngleRestraint::molecules() const
// {
//     Molecules mols;

//     for (int i = 0; i < this->nPoints(); ++i)
//     {
//         mols += p[i].read().molecules();
//     }

//     return mols;
// }

// /** Return whether or not this restraint affects the molecule
//     with number 'molnum' */
// bool AngleRestraint::contains(MolNum molnum) const
// {
//     return p[0].read().contains(molnum) or p[1].read().contains(molnum) or p[2].read().contains(molnum);
// }

// /** Return whether or not this restraint affects the molecule
//     with ID 'molid' */
// bool AngleRestraint::contains(const MolID &molid) const
// {
//     return p[0].read().contains(molid) or p[1].read().contains(molid) or p[2].read().contains(molid);
// }

// /** Return whether or not this restraint involves any of the molecules
//     that are in the forcetable 'forcetable' */
// bool AngleRestraint::usesMoleculesIn(const ForceTable &forcetable) const
// {
//     return p[0].read().usesMoleculesIn(forcetable) or p[1].read().usesMoleculesIn(forcetable) or
//            p[2].read().usesMoleculesIn(forcetable);
// }

// /** Return whether or not this restraint involves any of the molecules
//     in 'molecules' */
// bool AngleRestraint::usesMoleculesIn(const Molecules &molecules) const
// {
//     return p[0].read().usesMoleculesIn(molecules) or p[1].read().usesMoleculesIn(molecules) or
//            p[2].read().usesMoleculesIn(molecules);
// }

// static Expression harmonicFunction(double force_constant)
// {
//     if (SireMaths::isZero(force_constant))
//         return 0;
//     else
//         return force_constant * pow(AngleRestraint::theta(), 2);
// }

// static Expression diffHarmonicFunction(double force_constant)
// {
//     if (SireMaths::isZero(force_constant))
//         return 0;
//     else
//         return (2 * force_constant) * AngleRestraint::theta();
// }

// /** Return a distance restraint that applies a harmonic potential between
//     the points 'point0' and 'point1' using a force constant 'force_constant' */
// AngleRestraint AngleRestraint::harmonic(const PointRef &point0, const PointRef &point1, const PointRef &point2,
//                                         const HarmonicAngleForceConstant &force_constant)
// {
//     return AngleRestraint(point0, point1, point2, ::harmonicFunction(force_constant),
//                           ::diffHarmonicFunction(force_constant));
// }

// static Expression halfHarmonicFunction(double force_constant, double angle)
// {
//     if (SireMaths::isZero(force_constant))
//         return 0;

//     else if (angle <= 0)
//         // this is just a harmonic function
//         return ::harmonicFunction(force_constant);

//     else
//     {
//         const Symbol &theta = AngleRestraint::theta();
//         return Conditional(GreaterThan(theta, angle), force_constant * pow(theta - angle, 2), 0);
//     }
// }

// static Expression diffHalfHarmonicFunction(double force_constant, double angle)
// {
//     if (SireMaths::isZero(force_constant))
//         return 0;

//     else if (angle <= 0)
//         // this is just a harmonic function
//         return ::diffHarmonicFunction(force_constant);

//     else
//     {
//         const Symbol &theta = AngleRestraint::theta();
//         return Conditional(GreaterThan(theta, angle), (2 * force_constant) * (theta - angle), 0);
//     }
// }

// /** Return a distance restraint that applied a half-harmonic potential
//     between the points 'point0' and 'point1' above a distance 'distance'
//     using a force constant 'force_constant' */
// AngleRestraint AngleRestraint::halfHarmonic(const PointRef &point0, const PointRef &point1, const PointRef &point2,
//                                             const Angle &angle, const HarmonicAngleForceConstant &force_constant)
// {
//     double acute_angle = acute(angle).to(radians);

//     return AngleRestraint(point0, point1, point2, ::halfHarmonicFunction(force_constant, acute_angle),
//                           ::diffHalfHarmonicFunction(force_constant, acute_angle));
// }
