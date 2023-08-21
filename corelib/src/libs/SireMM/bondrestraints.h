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
  *  at http://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREMM_BONDRESTRAINTS_H
#define SIREMM_BONDRESTRAINTS_H

#include "restraints.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class BondRestraint;
    class BondRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::BondRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::BondRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::BondRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::BondRestraints &);

namespace SireMM
{

    /** This class represents a single bond restraint between any two
     *  atoms in a system (or between the centroids of any two groups
     *  of atoms in a system)
     */
    class SIREMM_EXPORT BondRestraint
        : public SireBase::ConcreteProperty<BondRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::BondRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::BondRestraint &);

    public:
        BondRestraint();
        BondRestraint(qint64 atom0, qint64 atom1,
                      const SireUnits::Dimension::GeneralUnit &k,
                      const SireUnits::Dimension::Length &r0);

        BondRestraint(const QList<qint64> &atoms0,
                      const QList<qint64> &atoms1,
                      const SireUnits::Dimension::GeneralUnit &k,
                      const SireUnits::Dimension::Length &r0);

        BondRestraint(const BondRestraint &other);

        ~BondRestraint();

        BondRestraint &operator=(const BondRestraint &other);

        bool operator==(const BondRestraint &other) const;
        bool operator!=(const BondRestraint &other) const;

        BondRestraints operator+(const BondRestraint &other) const;
        BondRestraints operator+(const BondRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        BondRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        bool isAtomRestraint() const;
        bool isCentroidRestraint() const;

        qint64 atom0() const;
        qint64 atom1() const;

        QVector<qint64> atoms0() const;
        QVector<qint64> atoms1() const;

        SireUnits::Dimension::GeneralUnit k() const;
        SireUnits::Dimension::Length r0() const;

    private:
        /** The first set of atoms involved in the restraint */
        QVector<qint64> atms0;

        /** The second set of atoms involved in the restraint */
        QVector<qint64> atms1;

        /** The force constant */
        double _k;

        /** The equilibrium distance for the restraint */
        double _r0;
    };

    /** This class provides the information for a collection of bond
     *  restraints that can be added to a collection of molecues. Each
     *  restraint can act on a pair of particles or a pair of the
     *  centroids of two collections of particles.
     *  The restaints are spherically symmetric, and
     *  are simple harmonic potentials
     */
    class SIREMM_EXPORT BondRestraints
        : public SireBase::ConcreteProperty<BondRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::BondRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::BondRestraints &);

    public:
        BondRestraints();

        BondRestraints(const BondRestraint &restraint);
        BondRestraints(const QList<BondRestraint> &restraints);

        BondRestraints(const QString &name,
                       const BondRestraint &restraint);

        BondRestraints(const QString &name,
                       const QList<BondRestraint> &restraints);

        BondRestraints(const BondRestraints &other);

        ~BondRestraints();

        BondRestraints &operator=(const BondRestraints &other);

        bool operator==(const BondRestraints &other) const;
        bool operator!=(const BondRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        BondRestraints *clone() const;

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

        const BondRestraint &at(int i) const;
        const BondRestraint &operator[](int i) const;

        QList<BondRestraint> restraints() const;

        QList<BondRestraint> atomRestraints() const;
        QList<BondRestraint> centroidRestraints() const;

        void add(const BondRestraint &restraint);
        void add(const BondRestraints &restraints);

        BondRestraints &operator+=(const BondRestraint &restraint);
        BondRestraints &operator+=(const BondRestraints &restraints);

        BondRestraints operator+(const BondRestraint &restraint) const;
        BondRestraints operator+(const BondRestraints &restraints) const;

    private:
        /** The actual list of restraints*/
        QList<BondRestraint> r;
    };

}

Q_DECLARE_METATYPE(SireMM::BondRestraint)
Q_DECLARE_METATYPE(SireMM::BondRestraints)

SIRE_EXPOSE_CLASS(SireMM::BondRestraint)
SIRE_EXPOSE_CLASS(SireMM::BondRestraints)

SIRE_END_HEADER

#endif
