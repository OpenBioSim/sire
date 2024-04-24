/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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

#ifndef SIRESYSTEM_SYSTEMTRAJECTORY_H
#define SIRESYSTEM_SYSTEMTRAJECTORY_H

#include "SireMol/trajectory.h"
#include "SireMol/molecules.h"

#include "SireVol/space.h"
#include "SireUnits/dimensions.h"

#include <memory>

SIRE_BEGIN_HEADER

namespace SireSystem
{
    class SystemTrajectory;
    class MolSystemTrajectory;
}

SIRESYSTEM_EXPORT QDataStream &operator<<(QDataStream &, const SireSystem::SystemTrajectory &);
SIRESYSTEM_EXPORT QDataStream &operator>>(QDataStream &, SireSystem::SystemTrajectory &);

SIRESYSTEM_EXPORT QDataStream &operator<<(QDataStream &, const SireSystem::MolSystemTrajectory &);
SIRESYSTEM_EXPORT QDataStream &operator>>(QDataStream &, SireSystem::MolSystemTrajectory &);

namespace SireSystem
{
    class SystemFrames;

    /** This is a TrajectoryData object that is used to hold the trajectory
     *  data for all molecules in a system. This is used to both speed
     *  up the reading and writing of trajectory data for entire
     *  systems (e.g. during a dynamics simulation), plus also to
     *  support offloading of trajectory data from memory to disk
     *  as the simulation progresses
     */
    class SIRESYSTEM_EXPORT SystemTrajectory : public SireMol::TrajectoryData
    {
        friend QDataStream & ::operator<<(QDataStream &, const SystemTrajectory &);
        friend QDataStream & ::operator>>(QDataStream &, SystemTrajectory &);

        friend class MolSystemTrajectory;

    public:
        SystemTrajectory();
        SystemTrajectory(const QList<SireMol::MolNum> &molnums,
                         const SireMol::Molecules &mols,
                         const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SystemTrajectory(const SystemTrajectory &other);

        ~SystemTrajectory();

        SystemTrajectory &operator=(const SystemTrajectory &other);

        bool operator==(const SystemTrajectory &other) const;
        bool operator!=(const SystemTrajectory &other) const;

        static const char *typeName();

        const char *what() const;

        SystemTrajectory *clone() const;

        bool isLive() const;

        void clear();

        bool isCompatibleWith(const QList<SireMol::MolNum> &molnums,
                              const SireMol::Molecules &mols,
                              const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        void saveFrame(const SireMol::Molecules &mols,
                       const SireVol::Space &space,
                       SireUnits::Dimension::Time time,
                       const SireBase::Properties &props,
                       const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SireMol::TrajectoryDataPtr getTrajectory(SireMol::MolNum molnum) const;

        int nFrames() const;
        int nAtoms() const;

        QStringList filenames() const;

        SireMol::Frame getFrame(int i) const;
        SireMol::Frame getFrame(int i, const SireBase::LazyEvaluator &evaluator) const;

        bool isEditable() const;

    protected:
        bool _equals(const TrajectoryData &other) const;

    private:
        /** Shared pointer to the underlying trajectory frames */
        std::shared_ptr<SystemFrames> d;
    };

    /** This is the view of a SystemTrajectory that is used to
     *  access the trajectory data for a single molecule in the system
     */
    class SIRESYSTEM_EXPORT MolSystemTrajectory : public SireMol::TrajectoryData
    {
        friend QDataStream & ::operator<<(QDataStream &, const MolSystemTrajectory &);
        friend QDataStream & ::operator>>(QDataStream &, MolSystemTrajectory &);

    public:
        MolSystemTrajectory();

        MolSystemTrajectory(const SystemTrajectory &trajectory, SireMol::MolNum molnum);

        MolSystemTrajectory(const MolSystemTrajectory &other);

        ~MolSystemTrajectory();

        MolSystemTrajectory &operator=(const MolSystemTrajectory &other);

        bool operator==(const MolSystemTrajectory &other) const;
        bool operator!=(const MolSystemTrajectory &other) const;

        static const char *typeName();

        const char *what() const;

        MolSystemTrajectory *clone() const;

        void clear();

        int nFrames() const;
        int nAtoms() const;

        QStringList filenames() const;

        SireMol::Frame getFrame(int i) const;
        SireMol::Frame getFrame(int i, const SireBase::LazyEvaluator &evaluator) const;

        bool isEditable() const;

    protected:
        bool _equals(const TrajectoryData &other) const;

    private:
        /** Shared pointer to the underlying trajectory frames */
        std::shared_ptr<SystemFrames> d;

        /** The Molecule number for the molecule that owns this trajectory */
        SireMol::MolNum molnum;
    };
}

Q_DECLARE_METATYPE(SireSystem::SystemTrajectory)
Q_DECLARE_METATYPE(SireSystem::MolSystemTrajectory)

SIRE_END_HEADER

#endif
