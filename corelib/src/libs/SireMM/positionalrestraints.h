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

#ifndef SIREMM_POSITIONALRESTRAINTS_H
#define SIREMM_POSITIONALRESTRAINTS_H

#include "restraints.h"

#include "SireMaths/vector.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class PositionalRestraint;
    class PositionalRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::PositionalRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::PositionalRestraints &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::PositionalRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::PositionalRestraint &);

namespace SireMM
{
    /** This class provides information about a single positional restraint.
     *  This is spherically symmetric and can act on either a single particle,
     *  or the centroid of a group of particles. The restraints are either
     *  flat-bottom harmonics or harmonic potentials
     */
    class SIREMM_EXPORT PositionalRestraint
        : public SireBase::ConcreteProperty<PositionalRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::PositionalRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::PositionalRestraint &);

    public:
        PositionalRestraint();
        PositionalRestraint(qint64 atom, const SireMaths::Vector &position,
                            const SireUnits::Dimension::HarmonicBondConstant &k,
                            const SireUnits::Dimension::Length &r0);

        PositionalRestraint(const QList<qint64> &atoms,
                            const SireMaths::Vector &position,
                            const SireUnits::Dimension::HarmonicBondConstant &k,
                            const SireUnits::Dimension::Length &r0);

        PositionalRestraint(const PositionalRestraint &other);

        ~PositionalRestraint();

        PositionalRestraint &operator=(const PositionalRestraint &other);

        bool operator==(const PositionalRestraint &other) const;
        bool operator!=(const PositionalRestraint &other) const;

        PositionalRestraints operator+(const PositionalRestraint &other) const;
        PositionalRestraints operator+(const PositionalRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        PositionalRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        bool isAtomRestraint() const;
        bool isCentroidRestraint() const;

        qint64 atom() const;
        QVector<qint64> atoms() const;

        SireMaths::Vector position() const;

        SireUnits::Dimension::HarmonicBondConstant k() const;
        SireUnits::Dimension::Length r0() const;

    private:
        /** The atoms involved in the restraint */
        QVector<qint64> atms;

        /** The location of the restraint */
        SireMaths::Vector pos;

        /** The force constant */
        SireUnits::Dimension::HarmonicBondConstant _k;

        /** The flat-bottom width (0 of harmonic restraints) */
        SireUnits::Dimension::Length _r0;
    };

    /** This class provides the information for a collection of positional
     *  restraints that can be added to a collection of molecues. Each
     *  restraint can act on a particle or the centroid of a collection
     *  of particles. The restaints are spherically symmetric, and
     *  are either flat-bottom harmonics or harmonic potentials
     */
    class SIREMM_EXPORT PositionalRestraints
        : public SireBase::ConcreteProperty<PositionalRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::PositionalRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::PositionalRestraints &);

    public:
        PositionalRestraints();
        PositionalRestraints(const QString &name);

        PositionalRestraints(const PositionalRestraint &restraint);
        PositionalRestraints(const QList<PositionalRestraint> &restraints);

        PositionalRestraints(const QString &name,
                             const PositionalRestraint &restraint);

        PositionalRestraints(const QString &name,
                             const QList<PositionalRestraint> &restraints);

        PositionalRestraints(const PositionalRestraints &other);

        ~PositionalRestraints();

        PositionalRestraints &operator=(const PositionalRestraints &other);

        bool operator==(const PositionalRestraints &other) const;
        bool operator!=(const PositionalRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        PositionalRestraints *clone() const;

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

        const PositionalRestraint &at(int i) const;
        const PositionalRestraint &operator[](int i) const;

        QList<PositionalRestraint> restraints() const;

        QList<PositionalRestraint> atomRestraints() const;
        QList<PositionalRestraint> centroidRestraints() const;

        void add(const PositionalRestraint &restraint);
        void add(const PositionalRestraints &restraints);

        PositionalRestraints &operator+=(const PositionalRestraint &restraint);
        PositionalRestraints &operator+=(const PositionalRestraints &restraints);

        PositionalRestraints operator+(const PositionalRestraint &restraint) const;
        PositionalRestraints operator+(const PositionalRestraints &restraints) const;

        void setUsesPeriodicBoundaryConditions(bool use_pbc);
        bool getUsesPeriodicBoundaryConditions() const;

    private:
        /** The actual list of restraints*/
        QList<PositionalRestraint> r;

        /** Whether the restraints use periodic boundary conditions */
        bool use_pbc = true;
    };

}

Q_DECLARE_METATYPE(SireMM::PositionalRestraint)
Q_DECLARE_METATYPE(SireMM::PositionalRestraints)

SIRE_EXPOSE_CLASS(SireMM::PositionalRestraint)
SIRE_EXPOSE_CLASS(SireMM::PositionalRestraints)

SIRE_END_HEADER

#endif
