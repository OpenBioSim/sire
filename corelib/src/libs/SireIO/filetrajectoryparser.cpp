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

#include "SireIO/filetrajectoryparser.h"

#include "SireSystem/system.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/core.h"
#include "SireMol/mgname.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/molidx.h"
#include "SireMol/trajectory.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/generalunitproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/timeproperty.h"
#include "SireBase/releasegil.h"
#include "SireBase/progressbar.h"

#include "SireIO/errors.h"
#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"

#include <QDir>
#include <QFileInfo>
#include <QRegularExpression>

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMaths;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<FileTrajectoryParser> r_traj;

QDataStream &operator<<(QDataStream &ds, const FileTrajectoryParser &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);

    sds << traj.frame_parser << traj.filenames << traj.current_idx
        << static_cast<const MoleculeParser &>(traj);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, FileTrajectoryParser &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> traj.frame_parser >> traj.filenames >> traj.current_idx >> static_cast<MoleculeParser &>(traj);
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

/** Constructor */
FileTrajectoryParser::FileTrajectoryParser()
    : ConcreteProperty<FileTrajectoryParser, MoleculeParser>(),
      current_idx(-1)
{
}

/** Return the format name that is used to identify this file format within Sire */
QString FileTrajectoryParser::formatName() const
{
    return "TRAJECTORY_PARSER";
}

/** Return the suffixes that RST7 files will typically have */
QStringList FileTrajectoryParser::formatSuffix() const
{
    return QStringList();
}

/** Return a description of the file format */
QString FileTrajectoryParser::formatDescription() const
{
    return QObject::tr("Parser to load a trajectory from multiple files.");
}

std::tuple<bool, bool, qint64, SireUnits::Dimension::Time> get_info_from_name(const QString &filename)
{
    // search for the frame match first
    {
        static const auto frame_regexp = QRegularExpression("frame_(\\d+)_([\\d-]+)");

        auto it = frame_regexp.globalMatch(filename);

        if (it.hasNext())
        {
            auto last_match = it.next();

            while (it.hasNext())
            {
                last_match = it.next();
            }

            qint64 frame_index = last_match.capturedView(1).toInt();

            double time_in_ps = last_match.capturedView(2).toString().replace("-", ".").toDouble();

            return std::tuple<bool, bool, qint64, SireUnits::Dimension::Time>(true, true, frame_index, time_in_ps * SireUnits::picosecond);
        }
    }

    // we haven't got the frame match, so see if we can get the basic match
    {
        static const auto num_regexp = QRegularExpression("(\\d+)");

        auto it = num_regexp.globalMatch(filename);

        if (it.hasNext())
        {
            auto last_match = it.next();

            while (it.hasNext())
            {
                last_match = it.next();
            }

            qint64 frame_index = last_match.capturedView(1).toInt();

            return std::tuple<bool, bool, qint64, SireUnits::Dimension::Time>(true, false, frame_index, SireUnits::Dimension::Time(0));
        }
    }

    // we couldn't find any match
    return std::tuple<bool, bool, qint64, SireUnits::Dimension::Time>(false, false, 0, SireUnits::Dimension::Time(0));
}

/** Scan the files to work out how many frames there are,
 *  to find out how many files match, and also to parse
 *  the first frame so that we can create the frame_parser
 *  of the right type
 */
void FileTrajectoryParser::parse()
{
    QDir dir(this->filename());

    // find all of the files that have the same extension and
    // follow a clear numbering pattern
    auto files = dir.entryInfoList(QDir::Readable | QDir::Files);

    QHash<QString, QList<QFileInfo>> files_by_suffix;

    for (const auto &file : files)
    {
        files_by_suffix[file.suffix()].append(file);
    }

    QHash<qint64, QStringList> popular_suffixs;

    for (const auto &suffix : files_by_suffix.keys())
    {
        popular_suffixs[files_by_suffix[suffix].count()].append(suffix);
    }

    // go through the suffixes in the order of greatest number
    // of files first
    auto counts = popular_suffixs.keys();
    std::sort(counts.begin(), counts.end(), std::greater());

    QHash<qint64, QStringList> files_by_frameidx;

    for (const auto &count : counts)
    {
        bool all_match = false;

        for (const auto &suffix : popular_suffixs[count])
        {
            all_match = true;
            files_by_frameidx.clear();

            // go through all the filenames and see if there is a pattern,
            // i.e. there is a clear numbering system
            for (const auto &file : files_by_suffix[suffix])
            {
                auto frame_info = get_info_from_name(file.filePath());

                if (std::get<0>(frame_info))
                {
                    files_by_frameidx[std::get<2>(frame_info)].append(file.absoluteFilePath());
                }
                else
                {
                    // nope - not a numerically named file
                    all_match = false;
                    break;
                }
            }
        }

        // do we have a match?
        if (all_match)
        {
            break;
        }
    }

    // now reconstruct the list of filenames in this frame index order
    auto idxs = files_by_frameidx.keys();
    std::sort(idxs.begin(), idxs.end(), std::less());

    for (const auto &idx : idxs)
    {
        for (const auto &file : files_by_frameidx[idx])
        {
            filenames.append(file);
        }
    }

    if (filenames.isEmpty())
    {
        throw SireIO::parse_error(QObject::tr(
                                      "Could not find any files in directory '%1' that followed a "
                                      "predictable numeric naming pattern, e.g. frame_000.rst, "
                                      "frame_001.rst, frame_002.rst etc. Please add an integer "
                                      "frame number to your files so that they can be detected "
                                      "to be part of a trajectory.")
                                      .arg(this->filename()),
                                  CODELOC);
    }

    // let's now try to parse the first file...
    frame_parser = MoleculeParser::_pvt_parse(filenames.at(0), this->propertyMap());

    if (frame_parser.read().isEmpty())
    {
        throw SireIO::parse_error(QObject::tr(
                                      "Despite the files being named sequentially, the first file '%1' "
                                      "did not contain molecular information that could be read.")
                                      .arg(filenames.at(0)),
                                  CODELOC);
    }

    if (not frame_parser.read().isFrame())
    {
        throw SireIO::parse_error(QObject::tr(
                                      "The first file in the trajectory in directory '%1' (%2) does "
                                      "not contain frame-style molecular data. It was recognised as a "
                                      "%3 file, which contains only topology data. The trajectory "
                                      "parser can only read files that contain frame data.")
                                      .arg(this->filename())
                                      .arg(filenames.at(0))
                                      .arg(frame_parser.read().formatName()),
                                  CODELOC);
    }

    if (frame_parser.read().nFrames() != 1)
    {
        throw SireIO::parse_error(QObject::tr(
                                      "The first file of the trajectory in directory '%1' (%2) does not "
                                      "contain a single trajectory frame. The number of contained "
                                      "frames is %3. The file trajectory parser can only be used to "
                                      "create a trajectory from input files that contain just a single "
                                      "frame of molecular data.")
                                      .arg(this->filename())
                                      .arg(filenames.at(0))
                                      .arg(frame_parser.read().nFrames()),
                                  CODELOC);
    }

    current_idx = 0;
}

/** Construct by parsing the passed file */
FileTrajectoryParser::FileTrajectoryParser(const QString &f, const PropertyMap &map)
    : ConcreteProperty<FileTrajectoryParser, MoleculeParser>(map),
      current_idx(-1)
{
    this->setFilename(f);
    this->parse();
}

/** Construct by parsing the data in the passed text lines */
FileTrajectoryParser::FileTrajectoryParser(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<FileTrajectoryParser, MoleculeParser>(map),
      current_idx(-1)
{
    throw SireIO::parse_error(QObject::tr("You cannot create a FileTrajectoryParser from a set of text lines!"),
                              CODELOC);
}

/** Construct by extracting the necessary data from the passed System */
FileTrajectoryParser::FileTrajectoryParser(const System &system, const PropertyMap &map)
    : ConcreteProperty<FileTrajectoryParser, MoleculeParser>(system, map),
      current_idx(-1)
{
    throw SireIO::parse_error(QObject::tr("The FrameTrajectoryParser is read-only!"),
                              CODELOC);
}

/** Copy constructor */
FileTrajectoryParser::FileTrajectoryParser(const FileTrajectoryParser &other)
    : ConcreteProperty<FileTrajectoryParser, MoleculeParser>(other),
      frame_parser(other.frame_parser), filenames(other.filenames),
      current_idx(other.current_idx)
{
}

/** Destructor */
FileTrajectoryParser::~FileTrajectoryParser()
{
}

FileTrajectoryParser &FileTrajectoryParser::operator=(const FileTrajectoryParser &other)
{
    if (this != &other)
    {
        frame_parser = other.frame_parser;
        filenames = other.filenames;
        current_idx = other.current_idx;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool FileTrajectoryParser::operator==(const FileTrajectoryParser &other) const
{
    return frame_parser.read().equals(other.frame_parser.read()) and
           filenames == other.filenames and
           current_idx == other.current_idx and
           MoleculeParser::operator==(other);
}

bool FileTrajectoryParser::operator!=(const FileTrajectoryParser &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *FileTrajectoryParser::typeName()
{
    return QMetaType::typeName(qMetaTypeId<FileTrajectoryParser>());
}

const char *FileTrajectoryParser::what() const
{
    return FileTrajectoryParser::typeName();
}

bool FileTrajectoryParser::isFrame() const
{
    return true;
}

int FileTrajectoryParser::nFrames() const
{
    return filenames.count();
}

Frame FileTrajectoryParser::getFrame(int frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    if (frame == current_idx)
    {
        return frame_parser.read().getFrame(0);
    }

    auto this_frame_parser = frame_parser.read().construct(filenames.at(frame),
                                                           this->propertyMap());

    auto frame_info = get_info_from_name(filenames.at(frame));

    auto f = this_frame_parser.read().getFrame(0);

    if (std::get<0>(frame_info) and std::get<1>(frame_info))
    {
        // update the frame time
        f = Frame(f.coordinates(), f.velocities(), f.forces(),
                  f.space(), std::get<3>(frame_info), f.properties());
    }

    return f;
}

FileTrajectoryParser FileTrajectoryParser::operator[](int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    if (i == current_idx)
    {
        return *this;
    }

    FileTrajectoryParser ret(*this);

    ret.frame_parser = frame_parser.read().construct(filenames.at(i),
                                                     this->propertyMap());

    ret.current_idx = i;

    return ret;
}

QString FileTrajectoryParser::toString() const
{
    if (this->nAtoms() == 0)
    {
        return QObject::tr("FileTrajectoryParser::null");
    }
    else
    {
        return QObject::tr("FileTrajectoryParser( nAtoms() = %1, nFrames() = %2 )")
            .arg(this->nAtoms())
            .arg(this->nFrames());
    }
}

/** Parse from the passed file */
FileTrajectoryParser FileTrajectoryParser::parse(const QString &filename)
{
    return FileTrajectoryParser(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void FileTrajectoryParser::addToSystem(System &system, const PropertyMap &map) const
{
    frame_parser.read().addToSystem(system, map);
}

/** Return the number of atoms whose coordinates are contained in this restart file */
int FileTrajectoryParser::nAtoms() const
{
    return frame_parser.read().nAtoms();
}

/** Return the title */
QString FileTrajectoryParser::title() const
{
    return QString();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr FileTrajectoryParser::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(FileTrajectoryParser(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr FileTrajectoryParser::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(FileTrajectoryParser(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr FileTrajectoryParser::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(FileTrajectoryParser(system, map));
}

/** This is not a text file that should be cached
 *  (it is potentially massive)
 */
bool FileTrajectoryParser::isTextFile() const
{
    return false;
}

/** Write this to the file 'filename' */
QStringList FileTrajectoryParser::writeToFile(const QString &filename) const
{
    // this is not written to a file...
    return QStringList();
}
