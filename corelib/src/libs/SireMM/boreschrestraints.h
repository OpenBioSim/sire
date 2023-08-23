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

#ifndef SIREMM_BORESCH_RESTRAINTS_H
#define SIREMM_BORESCH_RESTRAINTS_H

#include "restraints.h"

#include "SireMaths/vector.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class BoreschRestraint;
    class BoreschRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::BoreschRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::BoreschRestraints &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::BoreschRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::BoreschRestraint &);

namespace SireMM
{
    /** This class provides information about a single Boresch restaint.
     *  This is a collection of distance, angle and torsion restraints
     *  that hold a ligand in a binding pose relative to a receptor
     */
    class SIREMM_EXPORT BoreschRestraint
        : public SireBase::ConcreteProperty<BoreschRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::BoreschRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::BoreschRestraint &);

    public:
        BoreschRestraint();

        BoreschRestraint(const QList<qint64> &receptor_atoms,
                         const QList<qint64> &ligand_atoms,
                         const SireUnits::Dimension::Length &r0,
                         const QVector<SireUnits::Dimension::Angle> &theta0,
                         const QVector<SireUnits::Dimension::Angle> &phi0,
                         const SireUnits::Dimension::HarmonicBondConstant &kr,
                         const QVector<SireUnits::Dimension::HarmonicAngleConstant> &ktheta,
                         const QVector<SireUnits::Dimension::HarmonicAngleConstant> &kphi);

        BoreschRestraint(const BoreschRestraint &other);

        ~BoreschRestraint();

        BoreschRestraint &operator=(const BoreschRestraint &other);

        bool operator==(const BoreschRestraint &other) const;
        bool operator!=(const BoreschRestraint &other) const;

        BoreschRestraints operator+(const BoreschRestraint &other) const;
        BoreschRestraints operator+(const BoreschRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        BoreschRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        QVector<qint64> receptorAtoms() const;
        QVector<qint64> ligandAtoms() const;

        SireUnits::Dimension::Length r0() const;
        QVector<SireUnits::Dimension::Angle> theta0() const;
        QVector<SireUnits::Dimension::Angle> phi0() const;

        SireUnits::Dimension::HarmonicBondConstant kr() const;
        QVector<SireUnits::Dimension::HarmonicAngleConstant> ktheta() const;
        QVector<SireUnits::Dimension::HarmonicAngleConstant> kphi() const;

    private:
        /** The three receptor atoms involved in the restraint */
        QVector<qint64> receptor_atms;

        /** The three ligand atoms involved in the restraint */
        QVector<qint64> ligand_atms;

        /** The equilibium bond length */
        SireUnits::Dimension::Length _r0;

        /** The equilibrium angle sizes */
        QVector<SireUnits::Dimension::Angle> _theta0;

        /** The equilibrium torsion sizes */
        QVector<SireUnits::Dimension::Angle> _phi0;

        /** The bond force constant */
        SireUnits::Dimension::HarmonicBondConstant _kr;

        /** The angle force constants */
        QVector<SireUnits::Dimension::HarmonicAngleConstant> _ktheta;

        /** The dihedral force constants */
        QVector<SireUnits::Dimension::HarmonicAngleConstant> _kphi;
    };

    /** This class provides the information for a collection of positional
     *  restraints that can be added to a collection of molecues. Each
     *  restraint can act on a particle or the centroid of a collection
     *  of particles. The restaints are spherically symmetric, and
     *  are either flat-bottom harmonics or harmonic potentials
     */
    class SIREMM_EXPORT BoreschRestraints
        : public SireBase::ConcreteProperty<BoreschRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::BoreschRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::BoreschRestraints &);

    public:
        BoreschRestraints();

        BoreschRestraints(const BoreschRestraint &restraint);
        BoreschRestraints(const QList<BoreschRestraint> &restraints);

        BoreschRestraints(const QString &name,
                          const BoreschRestraint &restraint);

        BoreschRestraints(const QString &name,
                          const QList<BoreschRestraint> &restraints);

        BoreschRestraints(const BoreschRestraints &other);

        ~BoreschRestraints();

        BoreschRestraints &operator=(const BoreschRestraints &other);

        bool operator==(const BoreschRestraints &other) const;
        bool operator!=(const BoreschRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        BoreschRestraints *clone() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;
        int nRestraints() const;

        const BoreschRestraint &at(int i) const;
        const BoreschRestraint &operator[](int i) const;

        QList<BoreschRestraint> restraints() const;

        void add(const BoreschRestraint &restraint);
        void add(const BoreschRestraints &restraints);

        BoreschRestraints &operator+=(const BoreschRestraint &restraint);
        BoreschRestraints &operator+=(const BoreschRestraints &restraints);

        BoreschRestraints operator+(const BoreschRestraint &restraint) const;
        BoreschRestraints operator+(const BoreschRestraints &restraints) const;

    private:
        /** The actual list of restraints*/
        QList<BoreschRestraint> r;
    };

}

Q_DECLARE_METATYPE(SireMM::BoreschRestraint)
Q_DECLARE_METATYPE(SireMM::BoreschRestraints)

SIRE_EXPOSE_CLASS(SireMM::BoreschRestraint)
SIRE_EXPOSE_CLASS(SireMM::BoreschRestraints)

SIRE_END_HEADER

#endif
