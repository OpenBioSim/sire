/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2026  Lester Hedges
  *
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published
  *  by the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  *
  *  For full details of the license please see the LICENSE file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the website
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREMM_FLATBOTTOMBONDRESTRAINTS_H
#define SIREMM_FLATBOTTOMBONDRESTRAINTS_H

#include "restraints.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class FlatBottomBondRestraint;
    class FlatBottomBondRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::FlatBottomBondRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::FlatBottomBondRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::FlatBottomBondRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::FlatBottomBondRestraints &);

namespace SireMM
{

    /** This class represents a single flat-bottomed bond restraint
     *  between any two atoms in a system (or between the centroids
     *  of any two groups of atoms in a system). The energy is
     *  k*(|r - r0| - rfb)^2 when |r - r0| > rfb, and zero otherwise.
     */
    class SIREMM_EXPORT FlatBottomBondRestraint
        : public SireBase::ConcreteProperty<FlatBottomBondRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::FlatBottomBondRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::FlatBottomBondRestraint &);

    public:
        FlatBottomBondRestraint();
        FlatBottomBondRestraint(qint64 atom0, qint64 atom1,
                                const SireUnits::Dimension::HarmonicBondConstant &k,
                                const SireUnits::Dimension::Length &r0,
                                const SireUnits::Dimension::Length &rfb);

        FlatBottomBondRestraint(const QList<qint64> &atoms0,
                                const QList<qint64> &atoms1,
                                const SireUnits::Dimension::HarmonicBondConstant &k,
                                const SireUnits::Dimension::Length &r0,
                                const SireUnits::Dimension::Length &rfb);

        FlatBottomBondRestraint(const FlatBottomBondRestraint &other);

        ~FlatBottomBondRestraint();

        FlatBottomBondRestraint &operator=(const FlatBottomBondRestraint &other);

        bool operator==(const FlatBottomBondRestraint &other) const;
        bool operator!=(const FlatBottomBondRestraint &other) const;

        FlatBottomBondRestraints operator+(const FlatBottomBondRestraint &other) const;
        FlatBottomBondRestraints operator+(const FlatBottomBondRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        FlatBottomBondRestraint *clone() const;

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
        SireUnits::Dimension::Length rfb() const;

    private:
        /** The first set of atoms involved in the restraint */
        QVector<qint64> atms0;

        /** The second set of atoms involved in the restraint */
        QVector<qint64> atms1;

        /** The force constant */
        SireUnits::Dimension::HarmonicBondConstant _k;

        /** The equilibrium distance for the restraint */
        SireUnits::Dimension::Length _r0;

        /** The flat-bottom radius about the equilibrium distance */
        SireUnits::Dimension::Length _rfb;
    };

    /** This class provides the information for a collection of
     *  flat-bottomed bond restraints that can be added to a collection
     *  of molecules. Each restraint can act on a pair of particles or
     *  a pair of the centroids of two collections of particles.
     */
    class SIREMM_EXPORT FlatBottomBondRestraints
        : public SireBase::ConcreteProperty<FlatBottomBondRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::FlatBottomBondRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::FlatBottomBondRestraints &);

    public:
        FlatBottomBondRestraints();
        FlatBottomBondRestraints(const QString &name);

        FlatBottomBondRestraints(const FlatBottomBondRestraint &restraint);
        FlatBottomBondRestraints(const QList<FlatBottomBondRestraint> &restraints);

        FlatBottomBondRestraints(const QString &name,
                                 const FlatBottomBondRestraint &restraint);

        FlatBottomBondRestraints(const QString &name,
                                 const QList<FlatBottomBondRestraint> &restraints);

        FlatBottomBondRestraints(const FlatBottomBondRestraints &other);

        ~FlatBottomBondRestraints();

        FlatBottomBondRestraints &operator=(const FlatBottomBondRestraints &other);

        bool operator==(const FlatBottomBondRestraints &other) const;
        bool operator!=(const FlatBottomBondRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        FlatBottomBondRestraints *clone() const;

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

        const FlatBottomBondRestraint &at(int i) const;
        const FlatBottomBondRestraint &operator[](int i) const;

        QList<FlatBottomBondRestraint> restraints() const;

        QList<FlatBottomBondRestraint> atomRestraints() const;
        QList<FlatBottomBondRestraint> centroidRestraints() const;

        void add(const FlatBottomBondRestraint &restraint);
        void add(const FlatBottomBondRestraints &restraints);

        FlatBottomBondRestraints &operator+=(const FlatBottomBondRestraint &restraint);
        FlatBottomBondRestraints &operator+=(const FlatBottomBondRestraints &restraints);

        FlatBottomBondRestraints operator+(const FlatBottomBondRestraint &restraint) const;
        FlatBottomBondRestraints operator+(const FlatBottomBondRestraints &restraints) const;

        void setUsesPbc(bool use_pbc);
        bool usesPbc() const;

    private:
        /** The actual list of restraints*/
        QList<FlatBottomBondRestraint> r;

        /** Whether the restraints use periodic boundary conditions */
        bool use_pbc = false;
    };

}

Q_DECLARE_METATYPE(SireMM::FlatBottomBondRestraint)
Q_DECLARE_METATYPE(SireMM::FlatBottomBondRestraints)

SIRE_EXPOSE_CLASS(SireMM::FlatBottomBondRestraint)
SIRE_EXPOSE_CLASS(SireMM::FlatBottomBondRestraints)

SIRE_END_HEADER

#endif
