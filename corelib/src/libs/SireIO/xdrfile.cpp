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

#include "xdrfile.h"

#include "SireMol/trajectory.h"

#include <QDebug>
#include <QMutexLocker>
#include <QFileInfo>

#include "SireVol/cartesian.h"
#include "SireVol/triclinicbox.h"
#include "SireVol/periodicbox.h"

#include "SireIO/errors.h"
#include "SireError/errors.h"

#include "SireUnits/units.h"

#include "SireBase/parallel.h"
#include "SireBase/releasegil.h"
#include "SireBase/progressbar.h"
#include "SireBase/properties.h"
#include "SireBase/numberproperty.h"

#include "third_party/xdrfile.h"
#include "third_party/xdrfile_trr.h"
#include "third_party/xdrfile_xtc.h"

using namespace SireIO;
using namespace SireVol;
using namespace SireMol;
using namespace SireBase;
using namespace SireUnits;

////////
//////// Implementation of XDRFile
////////

QMutex XDRFile::mutex;

XDRFile::XDRFile()
    : boost::noncopyable(),
      f(0), sz(0)
{
}

XDRFile::XDRFile(const QString &filename)
    : boost::noncopyable(), fname(filename), f(0), sz(0)
{
}

XDRFile::~XDRFile()
{
    if (f)
    {
        // assume the mutex has been unlocked as we
        // are being destroyed
        this->_lkr_close();
        f = 0;
    }
}

QString XDRFile::filename() const
{
    return fname;
}

bool XDRFile::open(QIODevice::OpenMode mode)
{
    QMutexLocker lkr(&mutex);
    return this->_lkr_open(mode);
}

bool XDRFile::_lkr_open(QIODevice::OpenMode mode)
{
    this->_lkr_close();

    QFileInfo fullname(this->filename());

    sz = 0;

    if (mode == QIODevice::ReadOnly)
    {
        if (not fullname.exists())
        {
            throw SireError::file_error(QObject::tr(
                                            "Cannot read the file '%1' as it doesn't appear to exist?")
                                            .arg(this->filename()),
                                        CODELOC);
        }

        if (not fullname.isReadable())
        {
            throw SireError::file_error(QObject::tr(
                                            "Cannot read the file '%1' as it exists, but is not readable")
                                            .arg(this->filename()),
                                        CODELOC);
        }

        f = xdrfile_open(fullname.absoluteFilePath().toUtf8().constData(), "r");

        if (f == 0)
            throw SireError::file_error(QObject::tr(
                                            "An unknown error occurred when trying to open the XDR file '%1'")
                                            .arg(this->filename()),
                                        CODELOC);

        // seek to the end of the file and get the position
        auto ok = xdr_seek(f, 0, SEEK_END);

        if (ok != exdrOK)
        {
            xdrfile_close(f);
            f = 0;

            throw SireError::file_error(QObject::tr(
                                            "Error opening XDR file '%1' - could not get the file size!")
                                            .arg(this->filename()),
                                        CODELOC);
        }

        sz = xdr_tell(f);

        // now go back to the start of the file
        ok = xdr_seek(f, 0, SEEK_SET);

        if (ok != exdrOK)
        {
            xdrfile_close(f);
            f = 0;
            sz = 0;

            throw SireError::file_error(QObject::tr(
                                            "Error returning to the start of the file when reading '%1'.")
                                            .arg(this->filename()),
                                        CODELOC);
        }
    }
    else if (mode == QIODevice::WriteOnly or mode == QIODevice::Append)
    {
        if (fullname.exists())
        {
            if (fullname.isDir())
            {
                throw SireError::file_error(QObject::tr(
                                                "Cannot write the file '%1' as a directory with this "
                                                "name already exists.")
                                                .arg(this->filename()),
                                            CODELOC);
            }
            else if (not fullname.isWritable())
            {
                throw SireError::file_error(QObject::tr(
                                                "Cannot write the file '%1' as a file with this name "
                                                "already exists and it is not writeable.")
                                                .arg(this->filename()),
                                            CODELOC);
            }
        }

        if (mode == QIODevice::WriteOnly)
        {
            f = xdrfile_open(fullname.absoluteFilePath().toUtf8().constData(), "w");
        }
        else
        {
            f = xdrfile_open(fullname.absoluteFilePath().toUtf8().constData(), "a");
        }

        if (f == 0)
            throw SireError::file_error(QObject::tr(
                                            "An unknown error occurred when trying to open the XDR file '%1'")
                                            .arg(this->filename()),
                                        CODELOC);
    }
    else
    {
        throw SireError::file_error(QObject::tr(
                                        "You can only open a XDR file in ReadOnly, WriteOnly or Append mode."),
                                    CODELOC);
    }

    return true;
}

void XDRFile::close()
{
    QMutexLocker lkr(&mutex);
    this->_lkr_close();
}

void XDRFile::_lkr_close()
{
    if (f)
    {
        xdrfile_close(f);
        f = 0;
        sz = 0;
    }
}

qint64 XDRFile::_lkr_size() const
{
    return sz;
}

qint64 XDRFile::size() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&mutex));
    return sz;
}

///////
/////// Frame buffer used for XDRFiles
///////

namespace SireIO
{
    namespace detail
    {
        class XDRFrameBuffer
        {
        public:
            XDRFrameBuffer(const Frame &frame)
                : current_frame(-1), natoms(0),
                  coords(0), vels(0), frcs(0),
                  lambda(0), time(0), step(0),
                  precision(1000),
                  has_box(false)
            {
                if (frame.nAtoms() == 0)
                    return;

                natoms = frame.nAtoms();

                if (frame.hasCoordinates())
                    coords = new rvec[natoms];

                if (frame.hasVelocities())
                    vels = new rvec[natoms];

                if (frame.hasForces())
                    frcs = new rvec[natoms];

                if (frame.space().isPeriodic())
                    has_box = true;
            }

            XDRFrameBuffer(int num_atoms, qint32 frame_type)
                : current_frame(-1), natoms(0),
                  coords(0), vels(0), frcs(0),
                  lambda(0), time(0), step(0),
                  precision(1000),
                  has_box(false)
            {
                if (num_atoms <= 0)
                    return;

                natoms = num_atoms;

                if (frame_type & XDRFile::COORDINATES)
                {
                    coords = new rvec[natoms];
                }

                if (frame_type & XDRFile::VELOCITIES)
                {
                    vels = new rvec[natoms];
                }

                if (frame_type & XDRFile::FORCES)
                {
                    frcs = new rvec[natoms];
                }

                if (frame_type & XDRFile::BOX)
                {
                    has_box = true;
                }
            }

            ~XDRFrameBuffer()
            {
                delete[] coords;
                delete[] vels;
                delete[] frcs;
            }

            void resize(qint32 frame_type)
            {
                if (frame_type & XDRFile::COORDINATES)
                {
                    if (coords == 0)
                        coords = new rvec[natoms];
                }
                else if (coords != 0)
                {
                    delete[] coords;
                    coords = 0;
                }

                if (frame_type & XDRFile::VELOCITIES)
                {
                    if (vels == 0)
                        vels = new rvec[natoms];
                }
                else if (vels != 0)
                {
                    delete[] vels;
                    vels = 0;
                }

                if (frame_type & XDRFile::FORCES)
                {
                    if (frcs == 0)
                        frcs = new rvec[natoms];
                }
                else if (frcs != 0)
                {
                    delete[] frcs;
                    frcs = 0;
                }

                if (frame_type & XDRFile::BOX)
                {
                    has_box = true;
                }
                else
                {
                    has_box = false;
                }
            }

            int current_frame;
            int natoms;

            rvec *coords;
            rvec *vels;
            rvec *frcs;
            matrix box;

            float lambda;
            float time;
            int step;

            float precision;

            bool has_box;
        };
    }
}

////////
//////// Implementation of TRRFile
////////

TRRFile::TRRFile()
    : XDRFile(),
      natoms(0), nframes(0), bytes_per_frame(0), frame_type(0)
{
}

TRRFile::TRRFile(const QString &filename)
    : XDRFile(filename),
      natoms(0), nframes(0), bytes_per_frame(0), frame_type(0)
{
}

TRRFile::~TRRFile()
{
}

int TRRFile::nAtoms() const
{
    return natoms;
}

int TRRFile::nFrames() const
{
    return nframes;
}

void TRRFile::_lkr_reset()
{
    frame_buffer.reset();

    natoms = 0;
    nframes = 0;

    bytes_per_frame = 0;
    frame_type = 0;

    seek_frame.clear();
}

void TRRFile::_lkr_writeBufferToFile()
{
    if (f == 0)
        throw SireError::io_error(QObject::tr(
                                      "Cannot save to a file that is not open!"),
                                  CODELOC);

    if (frame_buffer.get() == 0)
        throw SireError::io_error(QObject::tr(
                                      "There is no TRR frame to save to the file..."),
                                  CODELOC);

    int ok = 0;

    if (frame_buffer->has_box)
    {
        ok = write_trr(f, frame_buffer->natoms, frame_buffer->step,
                       frame_buffer->time, frame_buffer->lambda,
                       frame_buffer->box,
                       frame_buffer->coords,
                       frame_buffer->vels,
                       frame_buffer->frcs);
    }
    else
    {
        ok = write_trr(f, frame_buffer->natoms, frame_buffer->step,
                       frame_buffer->time, frame_buffer->lambda,
                       0, // no box
                       frame_buffer->coords,
                       frame_buffer->vels,
                       frame_buffer->frcs);
    }

    if (ok != exdrOK)
        throw SireError::io_error(QObject::tr(
                                      "There was an error trying to write a TRR file to the file "
                                      "'%1'. The error was '%2'")
                                      .arg(this->filename())
                                      .arg(exdr_message[ok]),
                                  CODELOC);
}

void TRRFile::writeFrame(const Frame &frame, bool use_parallel)
{
    // create a frame buffer for this frame
    std::shared_ptr<detail::XDRFrameBuffer> buffer(new detail::XDRFrameBuffer(frame));

    // copy the data into this buffer
    if (buffer->has_box)
    {
        auto m = angstrom.to(nanometer) * frame.space().boxMatrix();

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                buffer->box[i][j] = m(i, j);
            }
        }
    }

    buffer->step = frame.property("step", NumberProperty(qint64(0))).asAnInteger();
    buffer->lambda = frame.property("lambda", NumberProperty(0.0)).asADouble();
    buffer->time = frame.time().to(picosecond);

    const int natoms = frame.nAtoms();

    if (natoms != buffer->natoms)
        throw SireError::program_bug(QObject::tr(
                                         "Memory not allocated!"),
                                     CODELOC);

    auto copy_coordinates = [&]()
    {
        if (not frame.hasCoordinates())
            return;

        const auto coords_data = frame.coordinates().constData();
        auto c_data = buffer->coords;

        if (c_data == 0)
            return;

        const double internal_to_units = (angstrom).to(nanometer);

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = coords_data[i];

                    c_data[i][0] = value.x() * internal_to_units;
                    c_data[i][1] = value.y() * internal_to_units;
                    c_data[i][2] = value.z() * internal_to_units;
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = coords_data[i];

                c_data[i][0] = value.x() * internal_to_units;
                c_data[i][1] = value.y() * internal_to_units;
                c_data[i][2] = value.z() * internal_to_units;
            }
        }
    };

    auto copy_velocities = [&]()
    {
        if (not frame.hasVelocities())
            return;

        const auto vels_data = frame.velocities().constData();
        auto v_data = buffer->vels;

        if (v_data == 0)
            return;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = vels_data[i];

                    v_data[i][0] = value.x().to(nanometer / picosecond);
                    v_data[i][1] = value.y().to(nanometer / picosecond);
                    v_data[i][2] = value.z().to(nanometer / picosecond);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = vels_data[i];

                v_data[i][0] = value.x().to(nanometer / picosecond);
                v_data[i][1] = value.y().to(nanometer / picosecond);
                v_data[i][2] = value.z().to(nanometer / picosecond);
            }
        }
    };

    auto copy_forces = [&]()
    {
        if (not frame.hasForces())
            return;

        const auto frcs_data = frame.forces().constData();
        auto f_data = buffer->frcs;

        if (f_data == 0)
            return;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = frcs_data[i];

                    f_data[i][0] = value.x().to(kilojoule / nanometer);
                    f_data[i][1] = value.y().to(kilojoule / nanometer);
                    f_data[i][2] = value.z().to(kilojoule / nanometer);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = frcs_data[i];

                f_data[i][0] = value.x().to(kilojoule / nanometer);
                f_data[i][1] = value.y().to(kilojoule / nanometer);
                f_data[i][2] = value.z().to(kilojoule / nanometer);
            }
        }
    };

    if (use_parallel)
    {
        tbb::parallel_invoke([&]()
                             { copy_coordinates(); },
                             [&]()
                             { copy_velocities(); },
                             [&]()
                             { copy_forces(); });
    }
    else
    {
        copy_coordinates();
        copy_velocities();
        copy_forces();
    }

    QMutexLocker lkr(&mutex);

    // replace the current buffer with this new buffer
    frame_buffer.reset();
    frame_buffer = buffer;

    // now write this to the file
    this->_lkr_writeBufferToFile();
}

void TRRFile::_lkr_readFrameIntoBuffer(int i)
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

    qint64 start_position = i * bytes_per_frame;
    qint32 ftype = frame_type;

    if (not seek_frame.isEmpty())
    {
        start_position = std::get<0>(seek_frame.at(i));
        ftype = std::get<1>(seek_frame.at(i));
    }

    if (frame_buffer.get() != 0)
    {
        if (frame_buffer->current_frame == i)
        {
            // we already have this frame in the buffer :-)
            return;
        }

        frame_buffer->resize(ftype);
    }
    else
    {
        frame_buffer.reset(new detail::XDRFrameBuffer(natoms, ftype));
    }

    // seek to the start position for this frame
    auto ok = xdr_seek(f, start_position, SEEK_SET);

    if (ok != exdrOK)
    {
        throw SireError::io_error(QObject::tr(
                                      "Unable to seek to the start of frame %1 in TRR file %2. "
                                      "Is the file corrupt or has it changed since opening?")
                                      .arg(i)
                                      .arg(this->filename()),
                                  CODELOC);
    }

    // read the frame into the buffer
    // (has_prop is the way read_trr communicates back whether or
    //  not it has read in coordinates, velocities, forces)
    int has_prop = 0;

    if (frame_buffer->has_box)
    {
        ok = read_trr(f,
                      frame_buffer->natoms,
                      &(frame_buffer->step),
                      &(frame_buffer->time),
                      &(frame_buffer->lambda),
                      frame_buffer->box,
                      frame_buffer->coords,
                      frame_buffer->vels,
                      frame_buffer->frcs,
                      &has_prop);
    }
    else
    {
        ok = read_trr(f,
                      frame_buffer->natoms,
                      &(frame_buffer->step),
                      &(frame_buffer->time),
                      &(frame_buffer->lambda),
                      0,
                      frame_buffer->coords,
                      frame_buffer->vels,
                      frame_buffer->frcs,
                      &has_prop);
    }

    // we could check `has_prop` to see if we have read in what we
    // expected to... (leave that for another time)

    if (ok != exdrOK)
    {
        if (seek_frame.isEmpty())
        {
            // we could have been tripped up because this trajectory
            // isn't an uniform as we expected - we should manually
            // reindex the frames...
            this->_lkr_reindexFrames();
            this->_lkr_readFrameIntoBuffer(i);
        }
        else
        {
            throw SireError::io_error(QObject::tr(
                                          "Failed to read in data for frame %1 in TRR file %2. "
                                          "Is the file corrupt or has it changed since opening?")
                                          .arg(i)
                                          .arg(this->filename()),
                                      CODELOC);
        }
    }

    frame_buffer->current_frame = i;
}

Frame TRRFile::readFrame(int i, bool use_parallel) const
{
    TRRFile *nonconst_this = const_cast<TRRFile *>(this);

    QMutexLocker lkr(&(nonconst_this->mutex));

    // load the frame into the buffer
    nonconst_this->_lkr_readFrameIntoBuffer(i);

    if (frame_buffer.get() == 0)
        // no frame has been loaded?
        return Frame();

    // copy it out from the buffer to local storage
    SpacePtr space;
    auto time = double(frame_buffer->time) * picosecond;
    float lambda = frame_buffer->lambda;
    int step = frame_buffer->step;
    bool has_box = frame_buffer->has_box;

    Matrix box(0.0);

    if (has_box)
    {
        const matrix &m = frame_buffer->box;
        box = Matrix(m[0][0], m[0][1], m[0][2],
                     m[1][0], m[1][1], m[1][2],
                     m[2][0], m[2][1], m[2][2]);
    }

    // and also the coords / vels / forces data
    QVector<float> coords;
    QVector<float> vels;
    QVector<float> frcs;

    if (frame_buffer->coords != 0)
    {
        float *start = &(frame_buffer->coords[0][0]);
        coords = QVector<float>(start, start + (3 * natoms));
    }

    if (frame_buffer->vels != 0)
    {
        float *start = &(frame_buffer->vels[0][0]);
        vels = QVector<float>(start, start + (3 * natoms));
    }

    if (frame_buffer->frcs != 0)
    {
        float *start = &(frame_buffer->frcs[0][0]);
        frcs = QVector<float>(start, start + (3 * natoms));
    }

    // we've finished with the buffer - can release the mutex
    lkr.unlock();

    // now need to convert the coords, vels and frcs into the right units
    // Can do this in parallel if allowed :-)
    QVector<Vector> c;
    QVector<Velocity3D> v;
    QVector<Force3D> f;

    auto copy_coordinates = [&]()
    {
        if (coords.isEmpty())
            return;

        c = QVector<Vector>(natoms);
        auto c_data = c.data();
        auto coords_data = coords.constData();

        const double units_to_internal = (1 * nanometer).value();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    c_data[i] = Vector(coords_data[(3 * i) + 0] * units_to_internal,
                                       coords_data[(3 * i) + 1] * units_to_internal,
                                       coords_data[(3 * i) + 2] * units_to_internal);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                c_data[i] = Vector(coords_data[(3 * i) + 0] * units_to_internal,
                                   coords_data[(3 * i) + 1] * units_to_internal,
                                   coords_data[(3 * i) + 2] * units_to_internal);
            }
        }
    };

    auto copy_velocities = [&]()
    {
        if (vels.isEmpty())
            return;

        v = QVector<Velocity3D>(natoms);
        auto v_data = v.data();
        auto vels_data = vels.constData();

        const auto units_to_internal = (1 * nanometer / picosecond);

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    v_data[i] = Velocity3D(double(vels_data[(3 * i) + 0]) * units_to_internal,
                                           double(vels_data[(3 * i) + 1]) * units_to_internal,
                                           double(vels_data[(3 * i) + 2]) * units_to_internal);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                v_data[i] = Velocity3D(double(vels_data[(3 * i) + 0]) * units_to_internal,
                                       double(vels_data[(3 * i) + 1]) * units_to_internal,
                                       double(vels_data[(3 * i) + 2]) * units_to_internal);
            }
        }
    };

    auto copy_forces = [&]()
    {
        if (frcs.isEmpty())
            return;

        f = QVector<Force3D>(natoms);
        auto f_data = f.data();
        auto frcs_data = frcs.constData();

        const auto units_to_internal = (kilojoule / nanometer);

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    f_data[i] = Force3D(double(frcs_data[(3 * i) + 0]) * units_to_internal,
                                        double(frcs_data[(3 * i) + 1]) * units_to_internal,
                                        double(frcs_data[(3 * i) + 2]) * units_to_internal);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                f_data[i] = Force3D(double(frcs_data[(3 * i) + 0]) * units_to_internal,
                                    double(frcs_data[(3 * i) + 1]) * units_to_internal,
                                    double(frcs_data[(3 * i) + 2]) * units_to_internal);
            }
        }
    };

    if (use_parallel)
    {
        tbb::parallel_invoke([&]()
                             { copy_coordinates(); },
                             [&]()
                             { copy_velocities(); },
                             [&]()
                             { copy_forces(); });
    }
    else
    {
        copy_coordinates();
        copy_velocities();
        copy_forces();
    }

    // any additional properties
    Properties props;

    props.setProperty("lambda", NumberProperty(lambda));
    props.setProperty("step", NumberProperty(qint64(step)));

    // now construct the space
    if (has_box)
    {
        box *= (1 * nanometer).value();

        if (box.isDiagonal())
        {
            return Frame(c, v, f, PeriodicBox(box.diagonal()), time, props);
        }
        else
        {
            return Frame(c, v, f,
                         TriclinicBox(box.column0(), box.column1(), box.column2()),
                         time, props);
        }
    }
    else
    {
        return Frame(c, v, f, Cartesian(), time, props);
    }

    return Frame();
}

void TRRFile::_lkr_reindexFrames()
{
    if (f == 0)
    {
        _lkr_reset();
        return;
    }

    nframes = 0;
    bytes_per_frame = 0;
    seek_frame.clear();

    // go through and manually find the start of each frame
    auto ok = xdr_seek(f, 0, SEEK_SET);

    t_trnheader header;

    auto bar = SireBase::ProgressBar("Indexing TRR trajectory frames:");
    bar.setSpeedUnit("frames / second");

    bar.enter();

    while (ok == exdrOK)
    {
        qint64 frame_start = xdr_tell(f);

        int step;
        float t;
        float lambda;
        int has_prop;

        // read in the header first, so we know what type of frame this is
        ok = do_trnheader(f, true, &header);

        if (ok != exdrOK)
        {
            // we have finished reading this trajectory
            break;
        }

        if (header.natoms != natoms)
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "Cannot read in the TRR file '%1' as the number of atoms "
                                                    "changes between frames. The first frame had natoms = %2, "
                                                    "but frame %3 has natoms = %4. Trajectories with different "
                                                    "numbers of atoms per frame are not supported yet in sire.")
                                                    .arg(this->filename())
                                                    .arg(natoms)
                                                    .arg(header.natoms)
                                                    .arg(header.step),
                                                CODELOC);
        }

        qint32 ftype = 0;

        if (header.box_size > 0)
            ftype |= XDRFile::BOX;

        if (header.x_size > 0)
            ftype |= XDRFile::COORDINATES;

        if (header.v_size > 0)
            ftype |= XDRFile::VELOCITIES;

        if (header.f_size > 0)
            ftype |= XDRFile::FORCES;

        // now seek back to the start of this frame
        ok = xdr_seek(f, frame_start, SEEK_SET);

        if (ok != exdrOK)
        {
            // something went wrong - likely end of the trajectory?
            break;
        }

        // now read in all of the data
        ok = read_trr(f, natoms, &step, &t, &lambda,
                      0, 0, 0, 0, &has_prop);

        if (ok == exdrOK)
        {
            // this was a valid frame :-)
            seek_frame.append(std::make_tuple(frame_start, ftype));
            bar.tick();
        }
    }

    bar.success();

    nframes = seek_frame.count();
}

bool TRRFile::open(QIODevice::OpenMode mode)
{
    QMutexLocker lkr(&mutex);

    this->_lkr_reset();

    if (not this->_lkr_open(mode))
        return false;

    if (f == 0)
        throw SireError::program_bug(QObject::tr(
            "The file handle should not be null if the file opened correctly..."));

    if (mode != QIODevice::ReadOnly)
    {
        // we only want to write, and so this should be an empty file
        return true;
    }

    // read the first header to see if this is really a TRR file
    t_trnheader header;

    auto ok = do_trnheader(f, true, &header);

    if (ok != exdrOK)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "The file '%1' is not a valid TRR file. The error message is '%2'")
                                      .arg(this->filename())
                                      .arg(exdr_message[ok]),
                                  CODELOC);
    }

    natoms = header.natoms;

    if (natoms <= 0)
    {
        // there are no frames either...
        nframes = 0;
        return true;
    }

    // protect against a memory DDOS or file corruption
    if (natoms > 2048 * 2048)
    {
        qint64 natoms_tmp = natoms;
        this->_lkr_close();
        this->_lkr_reset();
        throw SireError::unsupported(QObject::tr(
                                         "natoms = %1. Reading trajectory files with more than %1 atoms is not supported.")
                                         .arg(natoms_tmp)
                                         .arg(2048 * 2048),
                                     CODELOC);
    }

    // does this frame contain coordinate, velocity or force data?
    frame_type = 0;

    if (header.box_size > 0)
        frame_type |= XDRFile::BOX;

    if (header.x_size > 0)
        frame_type |= XDRFile::COORDINATES;

    if (header.v_size > 0)
        frame_type |= XDRFile::VELOCITIES;

    if (header.f_size > 0)
        frame_type |= XDRFile::FORCES;

    // ok - we now seek back to the start of the file and read the
    // frames
    ok = xdr_seek(f, 0, SEEK_SET);

    if (ok != exdrOK)
    {
        this->_lkr_close();
        this->_lkr_reset();
        throw SireError::file_error(QObject::tr(
                                        "Unable to seek to the start of '%1'")
                                        .arg(this->filename()),
                                    CODELOC);
    }

    // read in the first frame to get the frame size
    int start_pos = xdr_tell(f);

    int step;
    float t;
    float lambda;
    int has_prop;

    ok = read_trr(f, natoms, &step, &t, &lambda,
                  0, 0, 0, 0, &has_prop);

    if (ok != exdrOK)
    {
        this->_lkr_close();
        this->_lkr_reset();
        throw SireIO::parse_error(QObject::tr(
                                      "Could not parse the first TRR frame from '%1'.")
                                      .arg(this->filename()),
                                  CODELOC);
    }

    int end_pos = xdr_tell(f);

    bytes_per_frame = end_pos - start_pos;

    bool must_manually_index = false;

    // see if we can work out the number of frames,
    // assuming they all have the same size
    nframes = sz / bytes_per_frame;

    if (nframes * bytes_per_frame != sz)
    {
        // no - the size doesn't match up - we will have to
        // manually index
        bytes_per_frame = 0;
        must_manually_index = true;
    }
    else
    {
        // Looks promising - we will check by
        // trying to read the last frame...
        ok = xdr_seek(f, (nframes - 1) * bytes_per_frame, SEEK_SET);

        if (ok != exdrOK)
        {
            // nope - couldn't seek to this last frame?
            must_manually_index = true;
            bytes_per_frame = 0;
        }
        else
        {
            // try to read this last frame
            ok = read_trr(f, natoms, &step, &t, &lambda,
                          0, 0, 0, 0, &has_prop);

            if (ok != exdrOK)
            {
                // nope - something went wrong reading this frame
                // We will need to manually index
                must_manually_index = true;
                bytes_per_frame = 0;
            }
            // else we have successfully read in this last frame
            // This means we can trust that the frames are all the same size
            // and have predictable (and calculatable) seek positions
        }
    }

    if (must_manually_index)
    {
        this->_lkr_reindexFrames();
    }

    return true;
}

////////
//////// Implementation of XTCFile
////////

XTCFile::XTCFile()
    : XDRFile(), natoms(0)
{
}

XTCFile::XTCFile(const QString &filename)
    : XDRFile(filename), natoms(0)
{
}

XTCFile::~XTCFile()
{
}

int XTCFile::nAtoms() const
{
    return natoms;
}

int XTCFile::nFrames() const
{
    return seek_frame.count();
}

void XTCFile::_lkr_reset()
{
    frame_buffer.reset();

    natoms = 0;
    seek_frame.clear();
}

void XTCFile::_lkr_writeBufferToFile()
{
    if (f == 0)
        throw SireError::io_error(QObject::tr(
                                      "Cannot save to a file that is not open!"),
                                  CODELOC);

    if (frame_buffer.get() == 0)
        throw SireError::io_error(QObject::tr(
                                      "There is no XTC frame to save to the file..."),
                                  CODELOC);

    int ok = 0;

    ok = write_xtc(f, frame_buffer->natoms, frame_buffer->step,
                   frame_buffer->time,
                   frame_buffer->box,
                   frame_buffer->coords,
                   frame_buffer->precision);

    if (ok != exdrOK)
        throw SireError::io_error(QObject::tr(
                                      "There was an error trying to write a XTC file to the file "
                                      "'%1'. The error was '%2'")
                                      .arg(this->filename())
                                      .arg(exdr_message[ok]),
                                  CODELOC);
}

void XTCFile::writeFrame(const Frame &frame, bool use_parallel)
{
    // create a frame buffer for this frame
    std::shared_ptr<detail::XDRFrameBuffer> buffer(new detail::XDRFrameBuffer(frame));

    // copy the data into this buffer
    if (buffer->has_box)
    {
        auto m = angstrom.to(nanometer) * frame.space().boxMatrix();

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                buffer->box[i][j] = m(i, j);
            }
        }
    }

    buffer->step = frame.property("step", NumberProperty(qint64(0))).asAnInteger();
    buffer->time = frame.time().to(picosecond);

    const int natoms = frame.nAtoms();

    if (natoms != buffer->natoms)
        throw SireError::program_bug(QObject::tr(
                                         "Memory not allocated!"),
                                     CODELOC);

    auto copy_coordinates = [&]()
    {
        if (not frame.hasCoordinates())
            return;

        const auto coords_data = frame.coordinates().constData();
        auto c_data = buffer->coords;

        if (c_data == 0)
            return;

        const double internal_to_units = (angstrom).to(nanometer);

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = coords_data[i];

                    c_data[i][0] = value.x() * internal_to_units;
                    c_data[i][1] = value.y() * internal_to_units;
                    c_data[i][2] = value.z() * internal_to_units;
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = coords_data[i];

                c_data[i][0] = value.x() * internal_to_units;
                c_data[i][1] = value.y() * internal_to_units;
                c_data[i][2] = value.z() * internal_to_units;
            }
        }
    };

    copy_coordinates();

    QMutexLocker lkr(&mutex);

    // replace the current buffer with this new buffer
    frame_buffer.reset();
    frame_buffer = buffer;

    // now write this to the file
    this->_lkr_writeBufferToFile();
}

void XTCFile::_lkr_readFrameIntoBuffer(int i)
{
    if (i < 0 or i >= this->nFrames())
    {
        throw SireError::invalid_index(QObject::tr(
                                           "Cannot read frame %1 as the number of frames is %2.")
                                           .arg(i)
                                           .arg(this->nFrames()),
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

    auto start_position = seek_frame.at(i);

    if (frame_buffer.get() != 0)
    {
        if (frame_buffer->current_frame == i)
        {
            // we already have this frame in the buffer :-)
            return;
        }
    }
    else
    {
        frame_buffer.reset(new detail::XDRFrameBuffer(natoms, XDRFile::COORDINATES | XDRFile::BOX));
    }

    // seek to the start position for this frame
    auto ok = xdr_seek(f, start_position, SEEK_SET);

    if (ok != exdrOK)
    {
        throw SireError::io_error(QObject::tr(
                                      "Unable to seek to the start of frame %1 in TRR file %2. "
                                      "Is the file corrupt or has it changed since opening?")
                                      .arg(i)
                                      .arg(this->filename()),
                                  CODELOC);
    }

    // read the frame into the buffer
    ok = read_xtc(f,
                  frame_buffer->natoms,
                  &(frame_buffer->step),
                  &(frame_buffer->time),
                  frame_buffer->box,
                  frame_buffer->coords,
                  &(frame_buffer->precision));

    if (ok != exdrOK)
    {
        throw SireError::io_error(QObject::tr(
                                      "Failed to read in data for frame %1 in XTC file %2. "
                                      "Is the file corrupt or has it changed since opening?")
                                      .arg(i)
                                      .arg(this->filename()),
                                  CODELOC);
    }

    frame_buffer->current_frame = i;
}

Frame XTCFile::readFrame(int i, bool use_parallel) const
{
    XTCFile *nonconst_this = const_cast<XTCFile *>(this);

    QMutexLocker lkr(&(nonconst_this->mutex));

    // load the frame into the buffer
    nonconst_this->_lkr_readFrameIntoBuffer(i);

    if (frame_buffer.get() == 0)
    {
        // no frame has been loaded?
        return Frame();
    }

    // copy it out from the buffer to local storage
    SpacePtr space;
    auto time = double(frame_buffer->time) * picosecond;
    int step = frame_buffer->step;
    bool has_box = frame_buffer->has_box;

    Matrix box(0.0);

    if (has_box)
    {
        const matrix &m = frame_buffer->box;
        box = Matrix(m[0][0], m[0][1], m[0][2],
                     m[1][0], m[1][1], m[1][2],
                     m[2][0], m[2][1], m[2][2]);
    }

    // and also the coords / vels / forces data
    QVector<float> coords;

    if (frame_buffer->coords != 0)
    {
        float *start = &(frame_buffer->coords[0][0]);
        coords = QVector<float>(start, start + (3 * natoms));
    }

    // we've finished with the buffer - can release the mutex
    lkr.unlock();

    // now need to convert the coords, vels and frcs into the right units
    // Can do this in parallel if allowed :-)
    QVector<Vector> c;

    auto copy_coordinates = [&]()
    {
        if (coords.isEmpty())
            return;

        c = QVector<Vector>(natoms);
        auto c_data = c.data();
        auto coords_data = coords.constData();

        const double units_to_internal = (1 * nanometer).value();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    c_data[i] = Vector(coords_data[(3 * i) + 0] * units_to_internal,
                                       coords_data[(3 * i) + 1] * units_to_internal,
                                       coords_data[(3 * i) + 2] * units_to_internal);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                c_data[i] = Vector(coords_data[(3 * i) + 0] * units_to_internal,
                                   coords_data[(3 * i) + 1] * units_to_internal,
                                   coords_data[(3 * i) + 2] * units_to_internal);
            }
        }
    };

    copy_coordinates();

    // any additional properties
    Properties props;

    props.setProperty("step", NumberProperty(qint64(step)));

    // now construct the space
    if (has_box)
    {
        box *= (1 * nanometer).value();

        if (box.isZero())
        {
            return Frame(c, QVector<Velocity3D>(), QVector<Force3D>(),
                         Cartesian(), time, props);
        }
        else if (box.isDiagonal())
        {
            return Frame(c, QVector<Velocity3D>(), QVector<Force3D>(),
                         PeriodicBox(box.diagonal()), time, props);
        }
        else
        {
            return Frame(c, QVector<Velocity3D>(), QVector<Force3D>(),
                         TriclinicBox(box.column0(), box.column1(), box.column2()),
                         time, props);
        }
    }
    else
    {
        return Frame(c, QVector<Velocity3D>(), QVector<Force3D>(),
                     Cartesian(), time, props);
    }

    return Frame();
}

bool XTCFile::open(QIODevice::OpenMode mode)
{
    auto gil = SireBase::release_gil();

    QMutexLocker lkr(&mutex);

    this->_lkr_reset();

    if (not this->_lkr_open(mode))
        return false;

    if (f == 0)
        throw SireError::program_bug(QObject::tr(
            "The file handle should not be null if the file opened correctly..."));

    if (mode != QIODevice::ReadOnly)
    {
        // we only want to write, and so this should be an empty file
        return true;
    }

    // read the first header to see if this is really a XTC file
    int local_natoms;
    int step;
    float time;

    auto ok = xtc_header(f, &local_natoms, &step, &time, true);

    if (ok != exdrOK)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "The file '%1' is not a valid XTC file. The error message is '%2'")
                                      .arg(this->filename())
                                      .arg(exdr_message[ok]),
                                  CODELOC);
    }

    natoms = local_natoms;

    if (natoms <= 0)
    {
        return true;
    }

    // protect against a memory DDOS or file corruption
    if (natoms > 2048 * 2048)
    {
        qint64 natoms_tmp = natoms;
        this->_lkr_close();
        this->_lkr_reset();
        throw SireError::unsupported(QObject::tr(
                                         "natoms = %1. Reading trajectory files with more than %1 atoms is not supported.")
                                         .arg(natoms_tmp)
                                         .arg(2048 * 2048),
                                     CODELOC);
    }

    // we must manually index the file as each frame takes up
    // a different number of bytes
    seek_frame.clear();

    // how big is the file...
    ok = xdr_seek(f, 0, SEEK_END);
    qint64 file_size = xdr_tell(f);

    // go through and manually find the start of each frame
    ok = xdr_seek(f, 0, SEEK_SET);

    frame_buffer.reset(new detail::XDRFrameBuffer(natoms, XDRFile::COORDINATES | XDRFile::BOX));

    ProgressBar bar(file_size, "Indexing XTC");
    bar.setSpeedUnit("bytes / s");

    bar = bar.enter();

    while (ok == exdrOK)
    {
        qint64 frame_start = xdr_tell(f);

        // read in the frame - this is all we can do, because each frame
        // has a different size, and there is no shortcut available to
        // just read the frame header
        ok = read_xtc(f,
                      frame_buffer->natoms,
                      &(frame_buffer->step),
                      &(frame_buffer->time),
                      frame_buffer->box,
                      frame_buffer->coords,
                      &(frame_buffer->precision));

        if (ok != exdrOK)
        {
            // we have finished reading this trajectory
            bar.success();
            break;
        }
        else
        {
            // this was a valid frame :-)
            seek_frame.append(frame_start);
            bar.setProgress(frame_start);
        }
    }

    return true;
}
