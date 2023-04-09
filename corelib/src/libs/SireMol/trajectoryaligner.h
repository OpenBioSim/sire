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

#ifndef SIREMOL_TRAJECTORYALIGNER_H
#define SIREMOL_TRAJECTORYALIGNER_H

#include "trajectory.h"
#include "selectorm.hpp"

#include "SireMaths/align.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class FrameTransform;
    class TrajectoryAligner;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::FrameTransform &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::FrameTransform &);

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::TrajectoryAligner &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::TrajectoryAligner &);

namespace SireMol
{
    /** This class represents a transformation that is needed to be
     *  performed for a specific frame of a trajectory
     */
    class SIREMOL_EXPORT FrameTransform : public SireBase::ConcreteProperty<FrameTransform, SireBase::Property>
    {

        friend QDataStream & ::operator<<(QDataStream &, const FrameTransform &);
        friend QDataStream & ::operator>>(QDataStream &, FrameTransform &);

    public:
        FrameTransform();
        FrameTransform(const SireBase::PropertyMap &map);
        FrameTransform(const SireMaths::Transform &transform,
                       const SireMaths::Vector &center,
                       const SireVol::Space &space,
                       int nsmooth,
                       bool wrap);

        FrameTransform(const FrameTransform &other);

        ~FrameTransform();

        FrameTransform &operator=(const FrameTransform &other);

        bool operator==(const FrameTransform &other) const;
        bool operator!=(const FrameTransform &other) const;

        virtual FrameTransform *clone() const;

        const char *what() const;

        static const char *typeName();

        QString toString() const;

        SireVol::SpacePtr apply(const SireVol::Space &space) const;

        SireMaths::Vector apply(const SireMaths::Vector &coords) const;
        QVector<SireMaths::Vector> apply(const QVector<SireMaths::Vector> &coords) const;

        SireMaths::Vector apply(const SireMaths::Vector &coords,
                                const SireVol::Space &space) const;
        QVector<SireMaths::Vector> apply(const QVector<SireMaths::Vector> &coords,
                                         const SireVol::Space &space) const;

        Frame apply(const Frame &frame) const;

        SireVol::SpacePtr reverse(const SireVol::Space &space) const;

        SireMaths::Vector reverse(const SireMaths::Vector &coords) const;
        QVector<SireMaths::Vector> reverse(const QVector<SireMaths::Vector> &coords) const;

        Frame reverse(const Frame &frame) const;

        const SireMaths::Transform &transform() const;
        const SireMaths::Vector &center() const;

        const SireVol::Space &space() const;

        int nSmooth() const;

        bool wrap() const;

        bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    private:
        SireMaths::Transform tform;
        SireMaths::Vector cent;
        SireVol::SpacePtr old_space;
        SireVol::SpacePtr new_space;
        qint32 smooth;
        bool autowrap;
    };

    /** This class can be used to generate the FrameTransform object
     *  needed to align each frame of a trajectory.
     */
    class SIREMOL_EXPORT TrajectoryAligner : public SireBase::ConcreteProperty<TrajectoryAligner, SireBase::Property>
    {

        friend QDataStream & ::operator<<(QDataStream &, const TrajectoryAligner &);
        friend QDataStream & ::operator>>(QDataStream &, TrajectoryAligner &);

    public:
        TrajectoryAligner();
        TrajectoryAligner(const SireMaths::Vector &center,
                          const SireBase::PropertyMap &map = SireBase::PropertyMap());
        TrajectoryAligner(const SelectorM<Atom> &atoms,
                          const SireBase::PropertyMap &map = SireBase::PropertyMap());

        TrajectoryAligner(const TrajectoryAligner &other);

        ~TrajectoryAligner();

        TrajectoryAligner &operator=(const TrajectoryAligner &other);

        bool operator==(const TrajectoryAligner &other) const;
        bool operator!=(const TrajectoryAligner &other) const;

        virtual TrajectoryAligner *clone() const;

        const char *what() const;

        static const char *typeName();

        QString toString() const;

        const SelectorM<Atom> &atoms() const;

        FrameTransform operator[](int i) const;

        QList<FrameTransform> operator[](const QList<qint64> &idxs) const;
        QList<FrameTransform> operator[](const SireBase::Slice &slice) const;

        int count() const;
        int size() const;

        int nFrames() const;

        Vector center() const;

        int nSmooth() const;

        bool wrap() const;

        bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    private:
        void _populate(const SireBase::PropertyMap &map);

        /** The atoms that will be aligned */
        SelectorM<Atom> atms;

        /** The coordinates of these atoms in the reference frame */
        QVector<SireMaths::Vector> refcoords;

        /** The geometric center of the coordinates */
        SireMaths::Vector cent;

        /** The number of frames to smooth over */
        qint32 nsmooth;

        /** Property map used to find properties */
        SireBase::PropertyMap map;

        /** Whether or not we should autowrap after alignment */
        bool autowrap;
    };

}

Q_DECLARE_METATYPE(SireMol::TrajectoryAligner)
Q_DECLARE_METATYPE(SireMol::FrameTransform)

SIRE_EXPOSE_CLASS(SireMol::TrajectoryAligner)
SIRE_EXPOSE_CLASS(SireMol::FrameTransform)

SIRE_END_HEADER

#endif
