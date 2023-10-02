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

#ifndef SIREIO_FILETRAJECTORYPARSER_H
#define SIREIO_FILETRAJECTORYPARSER_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

#include "SireMol/trajectory.h"

#include "SireVol/space.h"

#include <memory>

SIRE_BEGIN_HEADER

namespace SireIO
{
    class FileTrajectoryParser;
}

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::FileTrajectoryParser &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::FileTrajectoryParser &);

class FileTrajectoryParserFile;

namespace SireIO
{
    /** This class represents a parser that parses multiple files of the
     *  same format and assembles them into a single trajectory.
     *
     *  This parser is special as it is one-way (read) only, and is not
     *  registered to be externally visible
     */
    class SIREIO_EXPORT FileTrajectoryParser : public SireBase::ConcreteProperty<FileTrajectoryParser, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const FileTrajectoryParser &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, FileTrajectoryParser &);

    public:
        FileTrajectoryParser();
        FileTrajectoryParser(const QString &filename, const PropertyMap &map = PropertyMap());
        FileTrajectoryParser(const QStringList &lines, const PropertyMap &map = PropertyMap());
        FileTrajectoryParser(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        FileTrajectoryParser(const FileTrajectoryParser &other);

        ~FileTrajectoryParser();

        FileTrajectoryParser &operator=(const FileTrajectoryParser &other);

        bool operator==(const FileTrajectoryParser &other) const;
        bool operator!=(const FileTrajectoryParser &other) const;

        static const char *typeName();

        const char *what() const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;

        QString formatName() const;
        QString formatDescription() const;
        QStringList formatSuffix() const;

        FileTrajectoryParser operator[](int i) const;

        static FileTrajectoryParser parse(const QString &filename);

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

        /** The underlying parser used to read individual frames */
        MoleculeParserPtr frame_parser;

        /** The names of all of the files in the trajectory
            arranged in frame index order */
        QStringList filenames;

        /** The current frame index */
        qint64 current_idx;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::FileTrajectoryParser)

SIRE_EXPOSE_CLASS(SireIO::FileTrajectoryParser)

SIRE_END_HEADER

#endif
