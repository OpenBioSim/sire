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

#ifndef SIREMM_DIHEDRALRESTRAINT_H
#define SIREMM_DIHEDRALRESTRAINT_H

// #include "SireFF/point.h"

#include "restraints.h"

// #include "SireCAS/expression.h"
// #include "SireCAS/symbol.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class DihedralRestraint;
    class DihedralRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::DihedralRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::DihedralRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::DihedralRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::DihedralRestraints &);

namespace SireMM
{

    /** This class represents a single angle restraint between any three
     *  atoms in a system
     * @author Christopher Woods
     */
    class SIREMM_EXPORT DihedralRestraint
        : public SireBase::ConcreteProperty<DihedralRestraint, SireBase::Property>
    {

        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const SireMM::DihedralRestraint &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, SireMM::DihedralRestraint &);

    public:
        DihedralRestraint();
        DihedralRestraint(const QList<qint64> &atoms,
                       const SireUnits::Dimension::Angle &phi0,
                       const SireUnits::Dimension::HarmonicAngleConstant &kphi);

        DihedralRestraint(const DihedralRestraint &other);

        ~DihedralRestraint();

        DihedralRestraint &operator=(const DihedralRestraint &other);

        bool operator==(const DihedralRestraint &other) const;
        bool operator!=(const DihedralRestraint &other) const;

        DihedralRestraints operator+(const DihedralRestraint &other) const;
        DihedralRestraints operator+(const DihedralRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        DihedralRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        QVector<qint64> atoms() const;

        SireUnits::Dimension::Angle phi0() const;
        SireUnits::Dimension::HarmonicAngleConstant kphi() const;

    private:
        /** Atoms involved in the angle restraint */
        QVector<qint64> atms;

        /** Equilibrium angle */
        SireUnits::Dimension::Angle _phi0;

        /** Harmonic angle constant */
        SireUnits::Dimension::HarmonicAngleConstant _kphi;
    };

    /** This class represents a collection of angle restraints */
    class SIREMM_EXPORT DihedralRestraints
        : public SireBase::ConcreteProperty<DihedralRestraints, Restraints>
    {
        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const SireMM::DihedralRestraints &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, SireMM::DihedralRestraints &);

    public:
        DihedralRestraints();
        DihedralRestraints(const QString &name);

        DihedralRestraints(const DihedralRestraint &restraint);
        DihedralRestraints(const QList<DihedralRestraint> &restraints);

        DihedralRestraints(const QString &name,
                        const DihedralRestraint &restraint);

        DihedralRestraints(const QString &name,
                        const QList<DihedralRestraint> &restraints);

        DihedralRestraints(const DihedralRestraints &other);

        ~DihedralRestraints();

        DihedralRestraints &operator=(const DihedralRestraints &other);

        bool operator==(const DihedralRestraints &other) const;
        bool operator!=(const DihedralRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        DihedralRestraints *clone() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;
        int nRestraints() const;

        const DihedralRestraint &at(int i) const;
        const DihedralRestraint &operator[](int i) const;

        QList<DihedralRestraint> restraints() const;

        void add(const DihedralRestraint &restraint);
        void add(const DihedralRestraints &restraints);

        DihedralRestraints &operator+=(const DihedralRestraint &restraint);
        DihedralRestraints &operator+=(const DihedralRestraints &restraints);

        DihedralRestraints operator+(const DihedralRestraint &restraint) const;
        DihedralRestraints operator+(const DihedralRestraints &restraints) const;

    private:
        /** List of restraints */
        QList<DihedralRestraint> r;
    };
}

Q_DECLARE_METATYPE(SireMM::DihedralRestraint)
Q_DECLARE_METATYPE(SireMM::DihedralRestraints)

SIRE_EXPOSE_CLASS(SireMM::DihedralRestraint)
SIRE_EXPOSE_CLASS(SireMM::DihedralRestraints)
SIRE_END_HEADER

#endif