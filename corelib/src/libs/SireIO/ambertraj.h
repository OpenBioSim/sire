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

#ifndef SIREIO_AMBERTRAJ_H
#define SIREIO_AMBERTRAJ_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

#include "SireMol/trajectory.h"

#include "SireVol/space.h"

#include <memory>

SIRE_BEGIN_HEADER

namespace SireIO
{
    class AmberTraj;
}

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::AmberTraj &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::AmberTraj &);

class AmberTrajFile;

namespace SireIO
{
    /** This class represents an Amber-format trajectory file (ascii)
        currently supporting these files from Amber7 to Amber16.

        This file holds either coordinates OR velocities (but not both).
        It can hold box information, but not time information.

        It holds the coordinates or velocities as a dense block of
        ascii numbers. It is described in the amber manual under
        "AMBER trajectory (coordinate or velocity) file specification"

        https://ambermd.org/FileFormats.php#trajectory
    */
    class SIREIO_EXPORT AmberTraj : public SireBase::ConcreteProperty<AmberTraj, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const AmberTraj &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, AmberTraj &);

    public:
        AmberTraj();
        AmberTraj(const QString &filename, const PropertyMap &map = PropertyMap());
        AmberTraj(const QStringList &lines, const PropertyMap &map = PropertyMap());
        AmberTraj(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        AmberTraj(const AmberTraj &other);

        ~AmberTraj();

        AmberTraj &operator=(const AmberTraj &other);

        bool operator==(const AmberTraj &other) const;
        bool operator!=(const AmberTraj &other) const;

        static const char *typeName();

        const char *what() const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;

        QString formatName() const;
        QString formatDescription() const;
        QStringList formatSuffix() const;

        AmberTraj operator[](int i) const;

        static AmberTraj parse(const QString &filename);

        bool isFrame() const;

        int nFrames() const;

        SireMol::Frame getFrame(int i) const;

        QString title() const;

        int nAtoms() const;

        bool isTextFile() const;

        QStringList writeToFile(const QString &filename) const;

    protected:
        void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

    private:
        void parse();

        /** The title of the file */
        QString ttle;

        /** The current frame */
        SireMol::Frame current_frame;

        /** The number of frames */
        qint64 nframes;

        /** The current frame index */
        qint64 frame_idx;

        /** Pointer to the underlying file */
        std::shared_ptr<::AmberTrajFile> f;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::AmberTraj)

SIRE_EXPOSE_CLASS(SireIO::AmberTraj)

SIRE_END_HEADER

#endif
