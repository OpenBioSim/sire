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
    class TrajectoryAligner;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::TrajectoryAligner &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::TrajectoryAligner &);

namespace SireMol
{

    /** This class can be used to generate the SireMaths::Transform object
     *  needed to align each frame of a trajectory.
     */
    class SIREMOL_EXPORT TrajectoryAligner : public SireBase::ConcreteProperty<TrajectoryAligner, SireBase::Property>
    {

        friend QDataStream & ::operator<<(QDataStream &, const TrajectoryAligner &);
        friend QDataStream & ::operator>>(QDataStream &, TrajectoryAligner &);

    public:
        TrajectoryAligner();
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

        SireMaths::Transform operator[](int i) const;

        QList<SireMaths::Transform> operator[](const QList<qint64> &idxs) const;
        QList<SireMaths::Transform> operator[](const SireBase::Slice &slice) const;

        int count() const;
        int size() const;

        int nFrames() const;

        bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    private:
        /** The atoms that will be aligned */
        SelectorM<Atom> atms;

        /** The coordinates of these atoms in the reference frame */
        QVector<SireMaths::Vector> refcoords;

        /** Property map used to find properties */
        SireBase::PropertyMap map;
    };

}

Q_DECLARE_METATYPE(SireMol::TrajectoryAligner)

SIRE_EXPOSE_CLASS(SireMol::TrajectoryAligner)

SIRE_END_HEADER

#endif
