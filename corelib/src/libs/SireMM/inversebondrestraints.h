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

#ifndef SIREMM_INVERSEBONDRESTRAINTS_H
#define SIREMM_INVERSEBONDRESTRAINTS_H

#include "restraints.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class InverseBondRestraint;
    class InverseBondRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::InverseBondRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::InverseBondRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::InverseBondRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::InverseBondRestraints &);

namespace SireMM
{

    /** This class represents a single bond restraint between any two
     *  atoms in a system (or between the centroids of any two groups
     *  of atoms in a system)
     */
    class SIREMM_EXPORT InverseBondRestraint
        : public SireBase::ConcreteProperty<InverseBondRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::InverseBondRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::InverseBondRestraint &);

    public:
        InverseBondRestraint();
        InverseBondRestraint(qint64 atom0, qint64 atom1,
                      const SireUnits::Dimension::HarmonicBondConstant &k,
                      const SireUnits::Dimension::Length &r0);

        InverseBondRestraint(const QList<qint64> &atoms0,
                      const QList<qint64> &atoms1,
                      const SireUnits::Dimension::HarmonicBondConstant &k,
                      const SireUnits::Dimension::Length &r0);

        InverseBondRestraint(const InverseBondRestraint &other);

        ~InverseBondRestraint();

        InverseBondRestraint &operator=(const InverseBondRestraint &other);

        bool operator==(const InverseBondRestraint &other) const;
        bool operator!=(const InverseBondRestraint &other) const;

        InverseBondRestraints operator+(const InverseBondRestraint &other) const;
        InverseBondRestraints operator+(const InverseBondRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        InverseBondRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        bool isAtomRestraint() const;
        bool isCentroidRestraint() const;

        qint64 atom0() const;
        qint64 atom1() const;

        QVector<qint64> atoms0() const;
        QVector<qint64> atoms1() const;

        SireUnits::Dimension::HarmonicBondConstant k() const;
        SireUnits::Dimension::Length r0() const;

    private:
        /** The first set of atoms involved in the restraint */
        QVector<qint64> atms0;

        /** The second set of atoms involved in the restraint */
        QVector<qint64> atms1;

        /** The force constant */
        SireUnits::Dimension::HarmonicBondConstant _k;

        /** The equilibrium distance for the restraint */
        SireUnits::Dimension::Length _r0;
    };

    /** This class provides the information for a collection of bond
     *  restraints that can be added to a collection of molecues. Each
     *  restraint can act on a pair of particles or a pair of the
     *  centroids of two collections of particles.
     *  The restaints are spherically symmetric, and
     *  are simple harmonic potentials
     */
    class SIREMM_EXPORT InverseBondRestraints
        : public SireBase::ConcreteProperty<InverseBondRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::InverseBondRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::InverseBondRestraints &);

    public:
        InverseBondRestraints();
        InverseBondRestraints(const QString &name);

        InverseBondRestraints(const InverseBondRestraint &restraint);
        InverseBondRestraints(const QList<InverseBondRestraint> &restraints);

        InverseBondRestraints(const QString &name,
                       const InverseBondRestraint &restraint);

        InverseBondRestraints(const QString &name,
                       const QList<InverseBondRestraint> &restraints);

        InverseBondRestraints(const InverseBondRestraints &other);

        ~InverseBondRestraints();

        InverseBondRestraints &operator=(const InverseBondRestraints &other);

        bool operator==(const InverseBondRestraints &other) const;
        bool operator!=(const InverseBondRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        InverseBondRestraints *clone() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;
        int nRestraints() const;

        int nAtomRestraints() const;
        int nCentroidRestraints() const;

        bool hasAtomRestraints() const;
        bool hasCentroidRestraints() const;

        const InverseBondRestraint &at(int i) const;
        const InverseBondRestraint &operator[](int i) const;

        QList<InverseBondRestraint> restraints() const;

        QList<InverseBondRestraint> atomRestraints() const;
        QList<InverseBondRestraint> centroidRestraints() const;

        void add(const InverseBondRestraint &restraint);
        void add(const InverseBondRestraints &restraints);

        InverseBondRestraints &operator+=(const InverseBondRestraint &restraint);
        InverseBondRestraints &operator+=(const InverseBondRestraints &restraints);

        InverseBondRestraints operator+(const InverseBondRestraint &restraint) const;
        InverseBondRestraints operator+(const InverseBondRestraints &restraints) const;

        void setUsesPeriodicBoundaryConditions(bool use_pbc);
        bool getUsesPeriodicBoundaryConditions() const;

    private:
        /** The actual list of restraints*/
        QList<InverseBondRestraint> r;

        /** Whether the restraints use periodic boundary conditions */
        bool use_pbc = false;
    };

}

Q_DECLARE_METATYPE(SireMM::InverseBondRestraint)
Q_DECLARE_METATYPE(SireMM::InverseBondRestraints)

SIRE_EXPOSE_CLASS(SireMM::InverseBondRestraint)
SIRE_EXPOSE_CLASS(SireMM::InverseBondRestraints)

SIRE_END_HEADER

#endif
