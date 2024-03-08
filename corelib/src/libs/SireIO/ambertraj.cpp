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

#include "SireIO/ambertraj.h"
#include "SireIO/amberformat.h"
#include "SireIO/textfile.h"

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

#include <QFile>
#include <QTextStream>

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

static QStringList toLines(const QVector<Vector> &all_coords,
                           bool uses_parallel)
{
    const qint64 nats = all_coords.count();
    QVector<double> coords(nats * 3, 0.0);

    auto coords_data = coords.data();
    const auto all_coords_data = all_coords.constData();

    if (uses_parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nats), [&](const tbb::blocked_range<int> &r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                const Vector &atomcoords = all_coords_data[i];
                coords_data[(3 * i) + 0] = atomcoords.x();
                coords_data[(3 * i) + 1] = atomcoords.y();
                coords_data[(3 * i) + 2] = atomcoords.z();
            } });
    }
    else
    {
        for (int i = 0; i < nats; ++i)
        {
            const Vector &atomcoords = all_coords_data[i];
            coords_data[(3 * i) + 0] = atomcoords.x();
            coords_data[(3 * i) + 1] = atomcoords.y();
            coords_data[(3 * i) + 2] = atomcoords.z();
        }
    }

    return writeFloatData(coords, AmberFormat(AmberPrm::FLOAT, 10, 8, 3), 0, false, 'f');
}

/** This is a specialisation of TextFile for AmberTraj files */
class AmberTrajFile : public TextFile
{
public:
    AmberTrajFile() : TextFile(),
                      start_frame_pos(0), bytes_per_frame(0),
                      lines_per_frame(0),
                      nvalues(0), natoms(0),
                      nframes(0), has_box_dims(false)
    {
    }

    AmberTrajFile(const QString &filename)
        : TextFile(filename),
          start_frame_pos(0), bytes_per_frame(0),
          lines_per_frame(0), nvalues(0), natoms(0),
          nframes(0), has_box_dims(false)
    {
    }

    ~AmberTrajFile()
    {
    }

    bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly)
    {
        QMutexLocker lkr(&mutex);
        this->_lkr_reset();

        if (not this->_lkr_open(mode))
            return false;

        if (f == 0 or ts == 0)
            throw SireError::program_bug(QObject::tr(
                "The file handle should not be null if the file opened correctly..."));

        if ((mode & QIODevice::ReadOnly) == 0)
        {
            // we only want to write, and so this should be an empty file
            return true;
        }

        // go to the start of the file...
        if (not ts->seek(0))
            throw SireIO::parse_error(QObject::tr(
                                          "Can't seek to the start of the file?"),
                                      CODELOC);

        if (ts->atEnd())
            throw SireIO::parse_error(QObject::tr(
                                          "An empty file is not a valid AmberTraj! %1")
                                          .arg(this->filename()),
                                      CODELOC);

        // read in the title - this should not be more than 1024 characters
        // (using a large number to protect against issues, while also
        //  checking against deliberate attempts to break this parser)
        ttle = ts->readLine(1024).simplified();

        if (ts->atEnd())
            throw SireIO::parse_error(QObject::tr(
                                          "There was only the title line in %1. This is not a valid TRAJ file.")
                                          .arg(this->filename()),
                                      CODELOC);

        // this is the position of the first frame
        start_frame_pos = ts->pos();

        // now read in the first frame, trying to work out how big it is...
        nvalues = 0;
        natoms = 0;
        nframes = 0;
        bytes_per_frame = 0;
        has_box_dims = false;

        int nvals_per_line = 0;
        int nvals_last_line = -1;
        int nvals_per_frame = 0;

        auto end_of_frame = [&]()
        {
            // we must have read in the last line of what can only be
            // coordinate data for this frame
            nvals_per_frame += nvals_per_line;

            if (nvals_per_frame % 3 != 0)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "This does not look like a valid Amber Traj file as the "
                                              "number of values in the frame (%1) is not divisible "
                                              "by three.")
                                              .arg(nvals_per_frame),
                                          CODELOC);
            }

            if (natoms <= 0)
            {
                natoms = nvals_per_frame / 3;
            }
            else if (nvals_per_frame / 3 != natoms)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "This does not look like a valid Amber Traj file as the number "
                                              "of atoms in the first frame (%1) does not equal the number "
                                              "of atoms in a subsequent frame (%2)")
                                              .arg(natoms)
                                              .arg(nvals_per_frame / natoms),
                                          CODELOC);
            }

            nvals_per_frame = 0;

            nframes += 1;
        };

        // we will read the data, using the fact it is laid out in
        // blocks of 10 F8.3 values, with different blocks for the
        // atom coordinates (which must be 3 x atoms long) and
        // then box dimensions (a line of 3 x values long), repeating
        // for the number of frames, to work out the number of frame,
        // whether or not there are box dimensions, and the number
        // of atoms. This is not foolproof - there are edge cases
        // where we will read this as one big frame of a large number
        // of atoms, or we will count an atom as being a periodic
        // box (as we will lean towards interpreting a line with
        // three values as the box dimensions). In either of these
        // cases, we can fix the issue when we convert to a System

        lines_per_frame = 0;

        int nread_bytes = 0;
        int end_frame_pos = -1;

        bool check_one_more_line = false;

        int newline_length = -1;

        // lines are a maximum of 80 characters - we will use
        // a buffer that is a power of 2 above that
        const int max_line_length = 128;
        QString line;
        line.reserve(max_line_length);

        while (not ts->atEnd())
        {
            ts->readLineInto(&line, max_line_length);

            if (newline_length == -1)
            {
                // how many bytes for a newline (\n or \r\n?)
                newline_length = ts->pos() - start_frame_pos - line.count();
            }

            const auto length = line.count() + newline_length;
            nread_bytes += length;

            lines_per_frame += 1;

            // read each line - this is a 10F8.3 format
            nvals_per_line = 0;
            for (int j = 0; j < 10; ++j)
            {
                int pos = j * 8;

                bool ok = true;

                if (pos + 8 > length)
                {
                    ok = false;
                }
                else
                {
                    line.midRef(pos, 8).toDouble(&ok);
                }

                if (not ok)
                {
                    break;
                }

                nvals_per_line += 1;
            }

            nvalues += nvals_per_line;

            if (nvals_per_line == 10)
            {
                // we have read in a full set of what can only be coordinate data
                if (nframes > 0)
                {
                    // but this is the start of the next frame
                    nread_bytes -= length;
                    lines_per_frame -= 1;
                    end_frame_pos = ts->pos() - length;
                    break;
                }
                else if (check_one_more_line)
                {
                    // this was either coordinates followed by space,
                    // or just a block of coordinates. We will assume there
                    // is a space
                    nvals_per_line = -3;
                    end_of_frame();
                    has_box_dims = true;
                    nread_bytes -= length;
                    lines_per_frame -= 1;
                    end_frame_pos = ts->pos() - length;
                    break;
                }

                nvals_per_frame += 10;
            }
            else if (nvals_per_line == 3)
            {
                // we must have either reached the end of a frame, or are
                // reading in the periodic box dimensions
                if (nvals_last_line == 3)
                {
                    // this is the edge case where the number of coordinates
                    // happens to have 3 space on the last line of the frame.
                    // This line is the periodic box, and the previous line
                    // is the end of the frame
                    nvals_per_line = 0;
                    end_of_frame();
                    nvals_per_line = 3;
                    has_box_dims = true;

                    if (natoms == 1)
                        throw SireIO::parse_error(QObject::tr(
                                                      "This parser does not support the reading of single-atom "
                                                      "trajectories in Amber Traj format."),
                                                  CODELOC);

                    end_frame_pos = ts->pos();
                    break;
                }
                else if (nvals_last_line == 10)
                {
                    // this is the edge case where the number of coordinates
                    // happens to fill the whole line (10 values) and this
                    // line is the periodic box, OR the end of the line.

                    // We will only be able to tell if the next line
                    // has 3 or 10 values...
                    check_one_more_line = true;
                    nvals_per_frame += 3;
                }
                else
                {
                    // we must have already finished a frame
                    if (nvals_per_frame != 0)
                        throw SireIO::parse_error(QObject::tr(
                                                      "This does not look like a valid Amber Traj file as the "
                                                      "number of consecutive lines without 10 values looks wrong."),
                                                  CODELOC);

                    has_box_dims = true;
                    end_frame_pos = ts->pos();
                    break;
                }
            }
            else
            {
                if (nframes > 0)
                {
                    // this is the next frame for a small system!
                    nread_bytes -= length;
                    lines_per_frame -= 1;
                    break;
                }
                else if (check_one_more_line)
                {
                    nvals_per_frame = -3;
                    nread_bytes -= length;
                    lines_per_frame -= 1;
                    end_of_frame();
                    end_frame_pos = ts->pos() - length;
                    break;
                }

                end_of_frame();
                // we may need to read the space from the next line...
            }

            nvals_last_line = nvals_per_line;
        }

        if (nframes == 0)
        {
            // we haven't finished the number of frames yet
            if (check_one_more_line)
            {
                // this was either coordinates followed by space,
                // or just a block of coordinates. We will assume there
                // is a space
                nvals_per_line = -3;
                end_of_frame();
                has_box_dims = true;
                end_frame_pos = ts->pos();
            }
            else
            {
                end_of_frame();
            }
        }

        // we have now finished reading in the first frame
        bytes_per_frame = end_frame_pos - start_frame_pos;

        if (nread_bytes != bytes_per_frame)
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Disagreement over the number of read bytes... %1 vs %2. "
                                          "This indicates a program bug or IO error.")
                                          .arg(nread_bytes)
                                          .arg(bytes_per_frame),
                                      CODELOC);
        }

        // how many frames do we think we have?
        nframes = (this->_lkr_size() - start_frame_pos) / bytes_per_frame;

        if (this->_lkr_size() != (nframes * bytes_per_frame) + start_frame_pos)
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Strange file size (%1) for the Amber Traj file. It is not a clean "
                                          "multiple of the number of frames (%2) times the bytes per frame "
                                          "(%3) plus the title size (%4) = (%5)")
                                          .arg(this->_lkr_size())
                                          .arg(nframes)
                                          .arg(bytes_per_frame)
                                          .arg(start_frame_pos)
                                          .arg((nframes * bytes_per_frame) + start_frame_pos),
                                      CODELOC);
        }

        return true;
    }

    SireMol::Frame readFrame(int i, bool use_parallel = true) const
    {
        AmberTrajFile *nonconst_this = const_cast<AmberTrajFile *>(this);

        QMutexLocker lkr(&(nonconst_this->mutex));

        nonconst_this->_lkr_readFrameIntoBuffer(i);

        QStringList lines = frame_buffer;

        int nvalues_per_frame = 3 * natoms;

        if (has_box_dims)
            nvalues_per_frame += 3;

        lkr.unlock();

        if (lines.count() != lines_per_frame)
            throw SireIO::parse_error(QObject::tr(
                                          "Unexpected number of line for frame %1 from file %2. "
                                          "Expected %3 lines, but got %4.")
                                          .arg(i)
                                          .arg(this->filename())
                                          .arg(lines_per_frame)
                                          .arg(lines.count()),
                                      CODELOC);

        // this is a very simple format - just lots of doubles written
        // in a 10F8.3 format (not necessarily 10 numbers per line)
        // We will read them until we get to the specified start value,
        // and will then try to interpret them based on what we need

        QVector<double> values;

        values.reserve(nvalues_per_frame);

        for (const auto &line : lines)
        {
            const auto length = line.count();

            bool line_ok = false;

            for (int j = 0; j < 10; ++j)
            {
                int pos = j * 8;

                if (pos + 8 > length)
                    // nothing left on this line to read
                    break;

                bool ok = false;
                double val = line.midRef(pos, 8).toDouble(&ok);

                if (not ok)
                    // assume the rest of the line is corrupted
                    break;

                values.append(val);

                // we have read at least one value from this line
                line_ok = true;
            }

            if (not line_ok)
                break;
        }

        if (values.count() != nvalues_per_frame)
            throw SireIO::parse_error(QObject::tr(
                                          "Could not read the expected amount of data (%1) for frame %2 "
                                          "from the TRAJ file %3. Number of values read was only %4.")
                                          .arg(nvalues_per_frame)
                                          .arg(i)
                                          .arg(this->filename())
                                          .arg(values.count()),
                                      CODELOC);

        // ok, we now have the raw values - interpret them as coordinates
        auto coords = QVector<Vector>(natoms);
        auto coords_data = coords.data();

        const double *values_data = values.constData();

        for (int i = 0; i < natoms; ++i)
        {
            coords_data[i] = Vector(values_data[(3 * i) + 0],
                                    values_data[(3 * i) + 1],
                                    values_data[(3 * i) + 2]);
        }

        if (has_box_dims)
        {
            Vector box_dims(values_data[(3 * natoms) + 0],
                            values_data[(3 * natoms) + 1],
                            values_data[(3 * natoms) + 2]);

            if (box_dims.isZero())
            {
                // this is the infinite cartesian space
                return Frame(coords, Cartesian(), SireUnits::Dimension::Time(0));
            }
            else
            {
                return Frame(coords, PeriodicBox(box_dims), SireUnits::Dimension::Time(0));
            }
        }
        else
        {
            return Frame(coords, Cartesian(), SireUnits::Dimension::Time(0));
        }
    }

    void writeFrame(const SireMol::Frame &frame, bool use_parallel = true)
    {
        // convert the frame into a set of lines that are in TRAJ format
        auto lines = toLines(frame.coordinates(), use_parallel);

        QVector<double> box_data;

        if (frame.space().isA<PeriodicBox>())
        {
            const auto box = frame.space().asA<PeriodicBox>().dimensions();

            box_data = QVector<double>(3);
            box_data[0] = box.x();
            box_data[1] = box.y();
            box_data[2] = box.z();
        }
        else if (frame.space().isA<Cartesian>())
        {
            box_data = QVector<double>(3, 0.0);
        }
        else
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "You can only write PeriodicBox or Cartesian space information "
                                                    "to a AmberTraj file. You cannot write a %1.")
                                                    .arg(frame.space().toString()),
                                                CODELOC);
        }

        lines += writeFloatData(box_data, AmberFormat(AmberPrm::FLOAT, 10, 8, 3), 0, false, 'f');

        QMutexLocker lkr(&mutex);
        frame_buffer = lines;
        this->_lkr_writeBufferToFile();
    }

    int nAtoms() const
    {
        return natoms;
    }

    int nFrames() const
    {
        return nframes;
    }

    QString title() const
    {
        return ttle;
    }

    void writeLine(const QString &line)
    {
        QMutexLocker lkr(&mutex);

        if (f != 0 and ts != 0)
        {
            *ts << line << "\n";
        }
    }

private:
    void _lkr_reset()
    {
        frame_buffer.clear();

        natoms = 0;
        nframes = 0;

        bytes_per_frame = 0;

        ttle = QString();

        seek_frame.clear();
    }

    void _lkr_readFrameIntoBuffer(int i)
    {
        if (i < 0 or i >= nframes)
        {
            throw SireError::invalid_index(QObject::tr(
                                               "Cannot read frame %1 as the number of frames is %2.")
                                               .arg(i)
                                               .arg(nframes),
                                           CODELOC);
        }

        if (f == 0)
        {
            throw SireError::io_error(QObject::tr(
                                          "Unaable to read frame %1 for file %2 as the file is not open?")
                                          .arg(i)
                                          .arg(this->filename()),
                                      CODELOC);
        }

        if (ts == 0)
            ts = new QTextStream(f);

        int start_bytes = 0;

        if (bytes_per_frame > 0)
        {
            // there are the same number of bytes per frame,
            // so easy to calculate the seek position
            start_bytes = start_frame_pos + (i * bytes_per_frame);
        }
        else
        {
            start_bytes = seek_frame.at(i);
        }

        if (not ts->seek(start_bytes))
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Unable to seek to position %1 in TRAJ file %2. "
                                          "File size is %3.")
                                          .arg(this->filename())
                                          .arg(start_bytes)
                                          .arg(f->size()),
                                      CODELOC);
        }

        // now read in the required number of lines to the buffer
        frame_buffer.clear();

        const int max_line_length = 128;

        for (int j = 0; j < lines_per_frame; ++j)
        {
            if (ts->atEnd())
                throw SireIO::parse_error(QObject::tr(
                                              "Unable to read line %1 of frame %2 from TRAJ file %3.")
                                              .arg(j)
                                              .arg(i)
                                              .arg(this->filename()),
                                          CODELOC);

            const auto line = ts->readLine(max_line_length);
            frame_buffer.append(line);
        }
    }

    void _lkr_writeBufferToFile()
    {
        for (const auto &line : frame_buffer)
        {
            *ts << line << "\n";
        }
    }

    /** The title of this file */
    QString ttle;

    /** The current frame that has been read into the buffer */
    QStringList frame_buffer;

    /** The seek position of each frame - this is only
     *  used if the frames have different sizes
     */
    QList<qint64> seek_frame;

    /** The position in bytes in the file of the first frame */
    qint64 start_frame_pos;

    /** The size, in bytes, of each frame. This is 0
     *  if each frame has a different number of bytes
     */
    qint64 bytes_per_frame;

    /** The number of lines for each frame. The format requires
     *  that this is a constant
     */
    qint64 lines_per_frame;

    /** The number of values in the file in total */
    qint32 nvalues;

    /** The number of atoms */
    qint64 natoms;

    /** The number of frames */
    qint64 nframes;

    /** Whether or not this has box data */
    bool has_box_dims;
};

///////
///////
///////

static const RegisterMetaType<AmberTraj> r_traj;
const RegisterParser<AmberTraj> register_traj;

QDataStream &operator<<(QDataStream &ds, const AmberTraj &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);

    sds << traj.current_frame
        << traj.nframes << traj.frame_idx
        << static_cast<const MoleculeParser &>(traj);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AmberTraj &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> traj.current_frame >> traj.nframes >> traj.frame_idx >> static_cast<MoleculeParser &>(traj);
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

/** Constructor */
AmberTraj::AmberTraj()
    : ConcreteProperty<AmberTraj, MoleculeParser>(),
      nframes(0), frame_idx(0)
{
}

/** Return the format name that is used to identify this file format within Sire */
QString AmberTraj::formatName() const
{
    return "TRAJ";
}

/** Return the suffixes that RST7 files will typically have */
QStringList AmberTraj::formatSuffix() const
{
    static const QStringList suffixes = {"traj", "trj", "crd"};
    return suffixes;
}

/** Return a description of the file format */
QString AmberTraj::formatDescription() const
{
    return QObject::tr("Amber trajectory (ascii) coordinate or velocity files "
                       "supported from Amber 7 upwards.");
}

/** Scan the file to work out how many values there are,
 *  and to extract the title
 */
void AmberTraj::parse()
{
    f.reset(new AmberTrajFile(this->filename()));

    try
    {
        if (not f->open(QIODevice::ReadOnly))
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Failed to open Amber TRAJ %1")
                                          .arg(this->filename()),
                                      CODELOC);
        }

        nframes = f->nFrames();
        current_frame = f->readFrame(0, this->usesParallel());
        frame_idx = 0;
        ttle = f->title();

        this->setScore(f->nFrames() * current_frame.nAtoms());
    }
    catch (...)
    {
        this->setScore(0);
        f.reset();
        throw;
    }
}

/** Construct by parsing the passed file */
AmberTraj::AmberTraj(const QString &f, const PropertyMap &map)
    : ConcreteProperty<AmberTraj, MoleculeParser>(map),
      nframes(0), frame_idx(-1)
{
    this->setFilename(f);
    this->parse();
}

/** Construct by parsing the data in the passed text lines */
AmberTraj::AmberTraj(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<AmberTraj, MoleculeParser>(map),
      nframes(0), frame_idx(-1)
{
    throw SireIO::parse_error(QObject::tr("You cannot create a Amber Traj file from a set of text lines!"),
                              CODELOC);
}

/** Construct by extracting the necessary data from the passed System */
AmberTraj::AmberTraj(const System &system, const PropertyMap &map)
    : ConcreteProperty<AmberTraj, MoleculeParser>(system, map),
      nframes(1), frame_idx(0)
{
    ttle = system.name().value();

    if (ttle.length() > 128)
        ttle.truncate(128);
    else if (ttle.isEmpty())
        ttle = "TRAJ file create by sire";

    current_frame = MoleculeParser::createFrame(system, map);
}

/** Copy constructor */
AmberTraj::AmberTraj(const AmberTraj &other)
    : ConcreteProperty<AmberTraj, MoleculeParser>(other),
      ttle(other.ttle), current_frame(other.current_frame),
      nframes(other.nframes), frame_idx(other.frame_idx),
      f(other.f)
{
}

/** Destructor */
AmberTraj::~AmberTraj()
{
}

AmberTraj &AmberTraj::operator=(const AmberTraj &other)
{
    if (this != &other)
    {
        ttle = other.ttle;
        current_frame = other.current_frame;
        nframes = other.nframes;
        frame_idx = other.frame_idx;
        f = other.f;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool AmberTraj::operator==(const AmberTraj &other) const
{
    return MoleculeParser::operator==(other);
}

bool AmberTraj::operator!=(const AmberTraj &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *AmberTraj::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberTraj>());
}

const char *AmberTraj::what() const
{
    return AmberTraj::typeName();
}

bool AmberTraj::isFrame() const
{
    return true;
}

int AmberTraj::nFrames() const
{
    return nframes;
}

void AmberTraj::reorderLoadedFrame()
{
    this->current_frame = this->reorderFrame(this->current_frame);
}

Frame AmberTraj::getFrame(int frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    if (frame < 0)
        frame = 0;

    if (frame == frame_idx)
        return current_frame;

    if (f.get() == 0)
    {
        throw SireError::file_error(QObject::tr(
                                        "Somehow we don't have access to the underlying TRAJ text file?"),
                                    CODELOC);
    }

    return this->reorderFrame(f->readFrame(frame, this->usesParallel()));
}

AmberTraj AmberTraj::operator[](int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    AmberTraj ret(*this);

    ret.current_frame = this->getFrame(i);
    ret.frame_idx = i;

    return ret;
}

QString AmberTraj::toString() const
{
    if (this->nAtoms() == 0)
    {
        return QObject::tr("AmberTraj::null");
    }
    else
    {
        return QObject::tr("AmberTraj( title() = %1, nAtoms() = %2, nFrames() = %3 )")
            .arg(this->title())
            .arg(this->nAtoms())
            .arg(this->nFrames());
    }
}

/** Parse from the passed file */
AmberTraj AmberTraj::parse(const QString &filename)
{
    return AmberTraj(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void AmberTraj::addToSystem(System &system, const PropertyMap &map) const
{
    MoleculeParser::copyFromFrame(current_frame, system, map);

    // update the System fileformat property to record that it includes
    // data from this file format
    QString fileformat = this->formatName();

    PropertyName fileformat_property = map["fileformat"];

    try
    {
        QString last_format = system.property(fileformat_property).asA<StringProperty>().value();
        fileformat = QString("%1,%2").arg(last_format, fileformat);
    }
    catch (...)
    {
    }

    if (fileformat_property.hasSource())
    {
        system.setProperty(fileformat_property.source(), StringProperty(fileformat));
    }
    else
    {
        system.setProperty("fileformat", StringProperty(fileformat));
    }
}

/** Return the title of the file */
QString AmberTraj::title() const
{
    return ttle;
}

/** Return the number of atoms whose coordinates are contained in this restart file */
int AmberTraj::nAtoms() const
{
    return current_frame.nAtoms();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr AmberTraj::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberTraj(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr AmberTraj::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberTraj(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr AmberTraj::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberTraj(system, map));
}

/** This is not a text file that should be cached
 *  (it is potentially massive)
 */
bool AmberTraj::isTextFile() const
{
    return false;
}

/** Write this to the file 'filename' */
QStringList AmberTraj::writeToFile(const QString &filename) const
{
    if (this->nFrames() == 0 or this->nAtoms() == 0)
        return QStringList();

    auto gil = SireBase::release_gil();

    createDirectoryForFile(filename);

    AmberTrajFile outfile(filename);

    if (not outfile.open(QIODevice::WriteOnly | QIODevice::Truncate))
        throw SireError::file_error(QObject::tr(
                                        "Could not open %1 to write the Amber TRAJ file.")
                                        .arg(filename),
                                    CODELOC);

    // write the title (only once at the top of the file)
    outfile.writeLine(ttle);

    if (this->writingTrajectory())
    {
        const auto frames = this->framesToWrite();

        ProgressBar bar("Save TRAJ", frames.count());
        bar.setSpeedUnit("frames / s");

        bar = bar.enter();

        for (int i = 0; i < frames.count(); ++i)
        {
            const auto frame = this->createFrame(frames[i]);
            outfile.writeFrame(frame, usesParallel());
            bar.setProgress(i + 1);
        }

        bar.success();
    }
    else
    {
        outfile.writeFrame(current_frame, usesParallel());
    }

    outfile.close();

    return QStringList(filename);
}
