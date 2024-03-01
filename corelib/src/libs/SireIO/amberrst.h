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

#ifndef SIREIO_AMBERRST_H
#define SIREIO_AMBERRST_H

#include "moleculeparser.h"

#include "SireMol/trajectory.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
    class AmberRst;
}

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::AmberRst &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::AmberRst &);

class AmberRstFile;

namespace SireIO
{
    /** This class represents an Amber-format restart/coordinate file (binary),
        currently supporting these files from Amber7 to Amber16.

        This is a netcdf file format, which is described here;

        http://ambermd.org/netcdf/nctraj.xhtml

        @author Christopher Woods
    */
    class SIREIO_EXPORT AmberRst : public SireBase::ConcreteProperty<AmberRst, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const AmberRst &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, AmberRst &);

    public:
        AmberRst();
        AmberRst(const QString &filename, const PropertyMap &map = PropertyMap());
        AmberRst(const QStringList &lines, const PropertyMap &map = PropertyMap());
        AmberRst(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        AmberRst(const AmberRst &other);

        ~AmberRst();

        AmberRst &operator=(const AmberRst &other);

        bool operator==(const AmberRst &other) const;
        bool operator!=(const AmberRst &other) const;

        static const char *typeName();

        const char *what() const;

        AmberRst operator[](int i) const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;

        QString formatName() const;
        QString formatDescription() const;
        QStringList formatSuffix() const;

        static AmberRst parse(const QString &filename);

        bool isFrame() const;

        int nFrames() const;
        int count() const;
        int size() const;

        SireMol::Frame getFrame(int i) const;

        int nAtoms() const;

        bool isTextFile() const;

        QStringList writeToFile(const QString &filename) const;

    protected:
        void addToSystem(SireSystem::System &system, const PropertyMap &map) const;
        void reorderLoadedFrame();

    private:
        void parse();

        /** The current frame */
        SireMol::Frame current_frame;

        /** Any warnings that were raised when reading the file */
        QStringList parse_warnings;

        /** The number of frames in this file */
        qint64 nframes;

        /** The current frame index */
        qint64 frame_idx;

        /** Pointer to the underlying file */
        std::shared_ptr<AmberRstFile> f;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::AmberRst)

SIRE_EXPOSE_CLASS(SireIO::AmberRst)

SIRE_END_HEADER

#endif
