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
  *  You can contact the authors via the website
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREMM_MORSEPOTENTIALRESTRAINTS_H
#define SIREMM_MORSEPOTENTIALRESTRAINTS_H

#include "restraints.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class MorsePotentialRestraint;
    class MorsePotentialRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::MorsePotentialRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::MorsePotentialRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::MorsePotentialRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::MorsePotentialRestraints &);

namespace SireMM
{

    /** This class represents a single Morse restraint between any two
     *  atoms in a system
     */
    class SIREMM_EXPORT MorsePotentialRestraint
        : public SireBase::ConcreteProperty<MorsePotentialRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::MorsePotentialRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::MorsePotentialRestraint &);

    public:
        MorsePotentialRestraint();
        MorsePotentialRestraint(qint64 atom0, qint64 atom1,
                       const SireUnits::Dimension::HarmonicBondConstant &k,
                       const SireUnits::Dimension::Length &r0,
                       const SireUnits::Dimension::MolarEnergy &de);

        MorsePotentialRestraint(const QList<qint64> &atoms0,
                       const QList<qint64> &atoms1,
                       const SireUnits::Dimension::HarmonicBondConstant &k,
                       const SireUnits::Dimension::Length &r0,
                       const SireUnits::Dimension::MolarEnergy &de);

        MorsePotentialRestraint(const MorsePotentialRestraint &other);

        ~MorsePotentialRestraint();

        MorsePotentialRestraint &operator=(const MorsePotentialRestraint &other);

        bool operator==(const MorsePotentialRestraint &other) const;
        bool operator!=(const MorsePotentialRestraint &other) const;

        MorsePotentialRestraints operator+(const MorsePotentialRestraint &other) const;
        MorsePotentialRestraints operator+(const MorsePotentialRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        MorsePotentialRestraint *clone() const;

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
        SireUnits::Dimension::MolarEnergy de() const;

    private:
        /** The first set of atoms involved in the restraint */
        QVector<qint64> atms0;

        /** The second set of atoms involved in the restraint */
        QVector<qint64> atms1;

        /** The force constant */
        SireUnits::Dimension::HarmonicBondConstant _k;

        /** The equilibrium distance for the restraint */
        SireUnits::Dimension::Length _r0;

        /** The bond dissociation energy */
        SireUnits::Dimension::MolarEnergy _de;
    };

    /** This class provides the information for a collection of Morse potential
     *  restraints that can be added to a collection of molecules. Each
     *  restraint can act on a pair of particles.
     */
    class SIREMM_EXPORT MorsePotentialRestraints
        : public SireBase::ConcreteProperty<MorsePotentialRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::MorsePotentialRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::MorsePotentialRestraints &);

    public:
        MorsePotentialRestraints();
        MorsePotentialRestraints(const QString &name);

        MorsePotentialRestraints(const MorsePotentialRestraint &restraint);
        MorsePotentialRestraints(const QList<MorsePotentialRestraint> &restraints);

        MorsePotentialRestraints(const QString &name,
                        const MorsePotentialRestraint &restraint);

        MorsePotentialRestraints(const QString &name,
                        const QList<MorsePotentialRestraint> &restraints);

        MorsePotentialRestraints(const MorsePotentialRestraints &other);

        ~MorsePotentialRestraints();

        MorsePotentialRestraints &operator=(const MorsePotentialRestraints &other);

        bool operator==(const MorsePotentialRestraints &other) const;
        bool operator!=(const MorsePotentialRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        MorsePotentialRestraints *clone() const;

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

        const MorsePotentialRestraint &at(int i) const;
        const MorsePotentialRestraint &operator[](int i) const;

        QList<MorsePotentialRestraint> restraints() const;

        QList<MorsePotentialRestraint> atomRestraints() const;
        QList<MorsePotentialRestraint> centroidRestraints() const;

        void add(const MorsePotentialRestraint &restraint);
        void add(const MorsePotentialRestraints &restraints);

        MorsePotentialRestraints &operator+=(const MorsePotentialRestraint &restraint);
        MorsePotentialRestraints &operator+=(const MorsePotentialRestraints &restraints);

        MorsePotentialRestraints operator+(const MorsePotentialRestraint &restraint) const;
        MorsePotentialRestraints operator+(const MorsePotentialRestraints &restraints) const;

        void setUsesPbc(bool use_pbc);
        bool usesPbc() const;

    private:
        /** The actual list of restraints*/
        QList<MorsePotentialRestraint> r;

        /** Whether the restraints use periodic boundary conditions */
        bool use_pbc = false;
    };

}

Q_DECLARE_METATYPE(SireMM::MorsePotentialRestraint)
Q_DECLARE_METATYPE(SireMM::MorsePotentialRestraints)

SIRE_EXPOSE_CLASS(SireMM::MorsePotentialRestraint)
SIRE_EXPOSE_CLASS(SireMM::MorsePotentialRestraints)

SIRE_END_HEADER

#endif
