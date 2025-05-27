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

#ifndef SIREMM_RMSDRESTRAINTS_H
#define SIREMM_RMSDRESTRAINTS_H

#include "restraints.h"

#include "SireMaths/vector.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class RMSDRestraint;
    class RMSDRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::RMSDRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::RMSDRestraints &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::RMSDRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::RMSDRestraint &);


namespace SireMM
{
    /** This class provides information about a single RMSD restraint.
     *  The restraints are either flat-bottom harmonics or harmonic potentials
     */
    class SIREMM_EXPORT RMSDRestraint
        : public SireBase::ConcreteProperty<RMSDRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::RMSDRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::RMSDRestraint &);

    public:
        RMSDRestraint();
        RMSDRestraint(const QList<qint64> &atoms, 
                        const QList<SireMaths::Vector> &ref_positions,
                        const SireUnits::Dimension::HarmonicBondConstant &k,
                        const SireUnits::Dimension::Length &r0);

        RMSDRestraint(const RMSDRestraint &other);

        ~RMSDRestraint();

        RMSDRestraint &operator=(const RMSDRestraint &other);

        bool operator==(const RMSDRestraint &other) const;
        bool operator!=(const RMSDRestraint &other) const;

        RMSDRestraints operator+(const RMSDRestraint &other) const;
        RMSDRestraints operator+(const RMSDRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        RMSDRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        QList<qint64> atoms() const;
        QList<SireMaths::Vector> ref_positions() const;

        SireUnits::Dimension::HarmonicBondConstant k() const;
        SireUnits::Dimension::Length r0() const;

    private:
        /** The atoms involved in the restraint */
        QVector<qint64> atms;

        /** The reference positions */
        QVector<SireMaths::Vector> ref_pos;

        /** The force constant */
        SireUnits::Dimension::HarmonicBondConstant _k;

        /** The flat-bottom width (0 of harmonic restraints) */
        SireUnits::Dimension::Length _r0;
    };

    /** This class provides the information for a collection of RMSD
     *  restraints that can be added to a collection of molecues. Each
     *  restraint can act on a particle or the centroid of a collection
     *  of particles. The restaints are spherically symmetric, and
     *  are either flat-bottom harmonics or harmonic potentials
     */
    class SIREMM_EXPORT RMSDRestraints
        : public SireBase::ConcreteProperty<RMSDRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::RMSDRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::RMSDRestraints &);

    public:
        RMSDRestraints();
        RMSDRestraints(const QString &name);

        RMSDRestraints(const RMSDRestraint &restraint);
        RMSDRestraints(const QList<RMSDRestraint> &restraints);

        RMSDRestraints(const QString &name,
                             const RMSDRestraint &restraint);

        RMSDRestraints(const QString &name,
                             const QList<RMSDRestraint> &restraints);

        RMSDRestraints(const RMSDRestraints &other);

        ~RMSDRestraints();

        RMSDRestraints &operator=(const RMSDRestraints &other);

        bool operator==(const RMSDRestraints &other) const;
        bool operator!=(const RMSDRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        RMSDRestraints *clone() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;
        int nRestraints() const;

        const RMSDRestraint &at(int i) const;
        const RMSDRestraint &operator[](int i) const;

        QList<RMSDRestraint> restraints() const;

        void add(const RMSDRestraint &restraint);
        void add(const RMSDRestraints &restraints);

        RMSDRestraints &operator+=(const RMSDRestraint &restraint);
        RMSDRestraints &operator+=(const RMSDRestraints &restraints);

        RMSDRestraints operator+(const RMSDRestraint &restraint) const;
        RMSDRestraints operator+(const RMSDRestraints &restraints) const;

    private:
        /** The actual list of restraints*/
        QList<RMSDRestraint> r;
    };

}

Q_DECLARE_METATYPE(SireMM::RMSDRestraint)
Q_DECLARE_METATYPE(SireMM::RMSDRestraints)

SIRE_EXPOSE_CLASS(SireMM::RMSDRestraint)
SIRE_EXPOSE_CLASS(SireMM::RMSDRestraints)

SIRE_END_HEADER

#endif