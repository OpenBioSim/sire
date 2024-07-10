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

#ifndef SIREMM_ANGLERESTRAINT_H
#define SIREMM_ANGLERESTRAINT_H

// #include "SireFF/point.h"

#include "restraints.h"

// #include "SireCAS/expression.h"
// #include "SireCAS/symbol.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class AngleRestraint;
    class AngleRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::AngleRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::AngleRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::AngleRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::AngleRestraints &);

namespace SireMM
{

    /** This class represents a single angle restraint between any three
     *  atoms in a system
     * @author Christopher Woods
     */
    class SIREMM_EXPORT AngleRestraint
        : public SireBase::ConcreteProperty<AngleRestraint, SireBase::Property>
    {

        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const SireMM::AngleRestraint &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, SireMM::AngleRestraint &);

    public:
        AngleRestraint();
        AngleRestraint(const QList<qint64> &atoms,
                       const SireUnits::Dimension::Angle &theta0,
                       const SireUnits::Dimension::HarmonicAngleConstant &ktheta);

        AngleRestraint(const AngleRestraint &other);

        ~AngleRestraint();

        AngleRestraint &operator=(const AngleRestraint &other);

        bool operator==(const AngleRestraint &other) const;
        bool operator!=(const AngleRestraint &other) const;

        AngleRestraints operator+(const AngleRestraint &other) const;
        AngleRestraints operator+(const AngleRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        AngleRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        QVector<qint64> atoms() const;

        SireUnits::Dimension::Angle theta0() const;
        SireUnits::Dimension::HarmonicAngleConstant ktheta() const;

    private:
        /** Atoms involved in the angle restraint */
        QVector<qint64> atms;

        /** Equilibrium angle */
        SireUnits::Dimension::Angle _theta0;

        /** Harmonic angle constant */
        SireUnits::Dimension::HarmonicAngleConstant _ktheta;
    };

    /** This class represents a collection of angle restraints */
    class SIREMM_EXPORT AngleRestraints
        : public SireBase::ConcreteProperty<AngleRestraints, Restraints>
    {
        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const SireMM::AngleRestraints &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, SireMM::AngleRestraints &);

    public:
        AngleRestraints();
        AngleRestraints(const QString &name);

        AngleRestraints(const AngleRestraint &restraint);
        AngleRestraints(const QList<AngleRestraint> &restraints);

        AngleRestraints(const QString &name,
                        const AngleRestraint &restraint);

        AngleRestraints(const QString &name,
                        const QList<AngleRestraint> &restraints);

        AngleRestraints(const AngleRestraints &other);

        ~AngleRestraints();

        AngleRestraints &operator=(const AngleRestraints &other);

        bool operator==(const AngleRestraints &other) const;
        bool operator!=(const AngleRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        AngleRestraints *clone() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;
        int nRestraints() const;

        const AngleRestraint &at(int i) const;
        const AngleRestraint &operator[](int i) const;

        QList<AngleRestraint> restraints() const;

        void add(const AngleRestraint &restraint);
        void add(const AngleRestraints &restraints);

        AngleRestraints &operator+=(const AngleRestraint &restraint);
        AngleRestraints &operator+=(const AngleRestraints &restraints);

        AngleRestraints operator+(const AngleRestraint &restraint) const;
        AngleRestraints operator+(const AngleRestraints &restraints) const;

    private:
        /** List of restraints */
        QList<AngleRestraint> r;
    };
}

Q_DECLARE_METATYPE(SireMM::AngleRestraint)
Q_DECLARE_METATYPE(SireMM::AngleRestraints)

SIRE_EXPOSE_CLASS(SireMM::AngleRestraint)
SIRE_EXPOSE_CLASS(SireMM::AngleRestraints)
SIRE_END_HEADER

#endif
