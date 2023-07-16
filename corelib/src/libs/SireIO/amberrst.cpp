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

#include "SireIO/amberrst.h"
#include "SireIO/netcdffile.h"
#include "SireIO/xdrfile.h"

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
#include "SireBase/progressbar.h"
#include "SireBase/releasegil.h"
#include "SireBase/getinstalldir.h"

#include "SireIO/errors.h"
#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDebug>

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

//////////
////////// Implementation of AmberRstFrameBuffer
//////////

class AmberRstFrameBuffer
{
public:
    AmberRstFrameBuffer()
        : current_frame(-1), natoms(0)
    {
    }

    AmberRstFrameBuffer(const Frame &frame)
        : current_frame(-1), natoms(0),
          cell_angles(0), cell_lengths(0),
          time(0)
    {
        if (frame.nAtoms() == 0)
            return;

        natoms = frame.nAtoms();

        if (frame.hasCoordinates())
            coords = QVector<double>(3 * natoms);

        if (frame.hasVelocities())
            vels = QVector<double>(3 * natoms);

        if (frame.hasForces())
            frcs = QVector<double>(3 * natoms);

        cell_angles = QVector<double>(3, 0.0);
        cell_lengths = QVector<double>(3, 0.0);
    }

    ~AmberRstFrameBuffer()
    {
    }

    int current_frame;
    int natoms;

    QVector<double> coords;
    QVector<double> vels;
    QVector<double> frcs;

    QVector<double> cell_angles;
    QVector<double> cell_lengths;

    double time;
};

//////////
////////// Implementation of AmberRstFile
//////////

class AmberRstFile : public NetCDFFile
{
public:
    AmberRstFile();
    AmberRstFile(const QString &filename);
    ~AmberRstFile();

    bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly);

    SireMol::Frame readFrame(int i, bool use_parallel = true) const;
    void writeFrame(int frame_idx,
                    const SireMol::Frame &frame,
                    bool use_parallel = true);

    void writeHeader(const QString &title, int nframes, int natoms,
                     bool have_coords, bool have_vels, bool have_frcs,
                     bool have_space);

    int nAtoms() const;
    int nFrames() const;

private:
    void _lkr_reset();
    void _lkr_reindexFrames();
    void _lkr_readFrameIntoBuffer(int i);
    void _lkr_writeBufferToFile();

    /** The current frame that has been read into the buffer */
    std::shared_ptr<AmberRstFrameBuffer> frame_buffer;

    /** The application that created this file */
    QString application;

    /** The specific program that created this file */
    QString program;

    /** The version of the program */
    QString program_version;

    /** Title for the file */
    QString title;

    /** The dimensions info in the file */
    QHash<QString, NetCDFDataInfo> dims;

    /** The number of atoms in the frame - we assume all
     *  frames have the same number of atoms
     */
    qint64 natoms;

    /** The number of frames in the file */
    qint64 nframes;

    /** Whether or not this is a restart file */
    bool is_restart;
};

AmberRstFile::AmberRstFile()
    : NetCDFFile(),
      natoms(0), nframes(0), is_restart(false)
{
}

AmberRstFile::AmberRstFile(const QString &filename)
    : NetCDFFile(filename),
      natoms(0), nframes(0), is_restart(false)
{
}

AmberRstFile::~AmberRstFile()
{
}

int AmberRstFile::nAtoms() const
{
    return natoms;
}

int AmberRstFile::nFrames() const
{
    return nframes;
}

void AmberRstFile::_lkr_reset()
{
    frame_buffer.reset();

    natoms = 0;
    nframes = 0;
}

QVector<double> float_to_double(const QVector<float> &values)
{
    QVector<double> ret(values.count());

    const auto values_data = values.constData();
    auto ret_data = ret.data();

    for (int i = 0; i < values.count(); ++i)
    {
        ret_data[i] = values_data[i];
    }

    return ret;
}

QVector<float> double_to_float(const QVector<double> &values)
{
    QVector<float> ret(values.count());

    const auto values_data = values.constData();
    auto ret_data = ret.data();

    for (int i = 0; i < values.count(); ++i)
    {
        ret_data[i] = values_data[i];
    }

    return ret;
}

void AmberRstFile::writeHeader(const QString &title, int nframes, int natoms,
                               bool have_coordinates, bool have_velocities,
                               bool have_forces, bool have_space)
{
    if (nframes <= 0 or natoms <= 0)
        return;

    // This will be a restart file if the number of frames is 1
    is_restart = (nframes == 1);

    // create the hash of all global data
    QHash<QString, QString> globals;

    // first create all of the global variables
    if (nframes <= 1)
    {
        globals.insert("Conventions", "AMBERRESTART");
        is_restart = true;
    }
    else
    {
        globals.insert("Conventions", "AMBER");
        is_restart = false;
    }

    globals.insert("ConventionVersion", "1.0");

    globals.insert("application", "sire");

    globals.insert("program", "AmberRst");

    globals.insert("programVersion", SireBase::getReleaseVersion());

    if (title.count() > 80)
    {
        globals.insert("title", title.mid(0, 80));
    }
    else if (not title.isEmpty())
    {
        globals.insert("title", title);
    }

    // now create the descriptions of all of the dimensions
    QHash<QString, NetCDFDataInfo> dimensions;

    // start off with the label variables
    {
        dimensions.insert("spatial",
                          NetCDFDataInfo(NetCDFDataInfo::get_nc_type<char>(),
                                         "spatial", {"spatial"}, {3},
                                         QHash<QString, QVariant>()));

        if (have_space)
        {
            dimensions.insert("cell_spatial",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<char>(),
                                             "cell_spatial",
                                             {"cell_spatial"},
                                             {3},
                                             QHash<QString, QVariant>()));
            dimensions.insert("cell_angular",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<char>(),
                                             "cell_angular",
                                             {"cell_angular", "label"},
                                             {3, 5},
                                             QHash<QString, QVariant>()));
        }
    }

    // now the time
    {
        QHash<QString, QVariant> attributes;

        attributes.insert("units", QString("picosecond"));

        if (is_restart)
        {
            dimensions.insert("time",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<double>(),
                                             "time",
                                             {}, {}, attributes));
        }
        else
        {
            dimensions.insert("time",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<float>(),
                                             "time",
                                             {"frame"},
                                             {nframes}, attributes));
        }
    }

    QHash<QString, QVariant> angstrom_units;
    angstrom_units.insert("units", "angstrom");

    // coordinates
    if (have_coordinates)
    {
        if (is_restart)
        {
            dimensions.insert("coordinates",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<double>(),
                                             "coordinates",
                                             {"atom", "spatial"},
                                             {natoms, 3}, angstrom_units));
        }
        else
        {
            dimensions.insert("coordinates",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<float>(),
                                             "coordinates",
                                             {"frame", "atom", "spatial"},
                                             {nframes, natoms, 3}, angstrom_units));
        }
    }

    // cell lengths and angles
    if (have_space)
    {
        QHash<QString, QVariant> degree_units;
        degree_units.insert("units", "degree");

        if (is_restart)
        {
            dimensions.insert("cell_lengths",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<double>(),
                                             "cell_lengths",
                                             {"cell_spatial"},
                                             {3}, angstrom_units));

            dimensions.insert("cell_angles",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<double>(),
                                             "cell_angles",
                                             {"cell_angular"},
                                             {3}, degree_units));
        }
        else
        {
            dimensions.insert("cell_lengths",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<float>(),
                                             "cell_lengths",
                                             {"frame", "cell_spatial"},
                                             {nframes, 3}, angstrom_units));

            dimensions.insert("cell_angles",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<float>(),
                                             "cell_angles",
                                             {"frame", "cell_angular"},
                                             {nframes, 3}, degree_units));
        }
    }

    // velocities
    if (have_velocities)
    {
        QHash<QString, QVariant> vel_attributes;
        vel_attributes.insert("units", "angstrom/picosecond");

        if (is_restart)
        {
            vel_attributes.insert("scale_factor", double(20.455));

            dimensions.insert("velocities",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<double>(),
                                             "velocities",
                                             {"atom", "spatial"},
                                             {natoms, 3}, vel_attributes));
        }
        else
        {
            vel_attributes.insert("scale_factor", float(20.455));

            dimensions.insert("velocities",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<float>(),
                                             "velocities",
                                             {"frame", "atom", "spatial"},
                                             {nframes, natoms, 3}, vel_attributes));
        }
    }

    // forces
    if (have_forces)
    {
        QHash<QString, QVariant> force_units;
        force_units.insert("units", "amu*angstrom/picosecond^2");

        if (is_restart)
        {
            dimensions.insert("forces",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<double>(),
                                             "forces",
                                             {"atom", "spatial"},
                                             {natoms, 3}, force_units));
        }
        else
        {
            dimensions.insert("forces",
                              NetCDFDataInfo(NetCDFDataInfo::get_nc_type<float>(),
                                             "forces",
                                             {"frame", "atom", "spatial"},
                                             {nframes, natoms, 3}, force_units));
        }
    }

    QMutexLocker lkr(NetCDFFile::globalMutex());
    NetCDFFile::_lkr_writeHeader(globals, dimensions);

    dims = NetCDFFile::_lkr_getVariablesInfo();

    const auto spatial = dims["spatial"];

    if (not spatial.isNull())
    {
        QVector<char> vals(3);
        vals[0] = 'x';
        vals[1] = 'y';
        vals[2] = 'z';

        NetCDFFile::_lkr_writeData(NetCDFData(spatial, vals));
    }

    if (have_space)
    {
        const auto cell_spatial = dims["cell_spatial"];

        if (not cell_spatial.isNull())
        {
            QVector<char> vals(3);
            vals[0] = 'a';
            vals[1] = 'b';
            vals[2] = 'c';

            NetCDFFile::_lkr_writeData(NetCDFData(cell_spatial, vals));
        }

        const auto cell_angular = dims["cell_angular"];

        if (not cell_angular.isNull())
        {
            QVector<char> vals(15);
            vals[0] = 'a';
            vals[1] = 'l';
            vals[2] = 'p';
            vals[3] = 'h';
            vals[4] = 'a';
            vals[5] = 'b';
            vals[6] = 'e';
            vals[7] = 't';
            vals[8] = 'a';
            vals[9] = ' ';
            vals[10] = 'g';
            vals[11] = 'a';
            vals[12] = 'm';
            vals[13] = 'm';
            vals[14] = 'a';

            NetCDFFile::_lkr_writeData(NetCDFData(cell_angular, vals));
        }
    }
}

void AmberRstFile::_lkr_writeBufferToFile()
{
    if (is_restart)
    {
        // start with the time
        const auto time = dims["time"];

        if (not time.isNull())
            NetCDFFile::_lkr_writeData(NetCDFData(time, QVector<double>(1, frame_buffer->time)));

        // now the coordinates
        const auto coordinates = dims["coordinates"];

        if (not(coordinates.isNull() or frame_buffer->coords.isEmpty()))
            NetCDFFile::_lkr_writeData(NetCDFData(coordinates,
                                                  frame_buffer->coords));

        // now the velocities
        const auto velocities = dims["velocities"];

        if (not(velocities.isNull() or frame_buffer->vels.isEmpty()))
        {
            NetCDFFile::_lkr_writeData(NetCDFData(velocities,
                                                  frame_buffer->vels));
        }

        // now the forces
        const auto forces = dims["forces"];

        if (not(forces.isNull() or frame_buffer->frcs.isEmpty()))
            NetCDFFile::_lkr_writeData(NetCDFData(forces,
                                                  frame_buffer->frcs));

        // now the box info
        const auto cell_lengths = dims["cell_lengths"];
        const auto cell_angles = dims["cell_angles"];

        if (not(cell_lengths.isNull() or cell_angles.isNull()))
        {
            NetCDFFile::_lkr_writeData(NetCDFData(cell_lengths,
                                                  frame_buffer->cell_lengths));
            NetCDFFile::_lkr_writeData(NetCDFData(cell_angles,
                                                  frame_buffer->cell_angles));
        }
    }
    else
    {
        const auto current_frame = frame_buffer->current_frame;

        // start with the time
        const auto time = dims["time"];

        if (not time.isNull())
            NetCDFFile::_lkr_writeData(NetCDFData(time,
                                                  time.hyperslab()[current_frame],
                                                  QVector<float>(1, frame_buffer->time)));

        // now the coordinates
        const auto coordinates = dims["coordinates"];

        if (not(coordinates.isNull() or frame_buffer->coords.isEmpty()))
            NetCDFFile::_lkr_writeData(NetCDFData(coordinates,
                                                  coordinates.hyperslab()[current_frame],
                                                  double_to_float(frame_buffer->coords)));

        // now the velocities
        const auto velocities = dims["velocities"];

        if (not(velocities.isNull() or frame_buffer->vels.isEmpty()))
            NetCDFFile::_lkr_writeData(NetCDFData(velocities,
                                                  velocities.hyperslab()[current_frame],
                                                  double_to_float(frame_buffer->vels)));

        // now the forces
        const auto forces = dims["forces"];

        if (not(forces.isNull() or frame_buffer->frcs.isEmpty()))
            NetCDFFile::_lkr_writeData(NetCDFData(forces,
                                                  forces.hyperslab()[current_frame],
                                                  double_to_float(frame_buffer->frcs)));

        // now the box info
        const auto cell_lengths = dims["cell_lengths"];
        const auto cell_angles = dims["cell_angles"];

        if (not(cell_lengths.isNull() or cell_angles.isNull()))
        {
            NetCDFFile::_lkr_writeData(NetCDFData(cell_lengths,
                                                  cell_lengths.hyperslab()[current_frame],
                                                  double_to_float(frame_buffer->cell_lengths)));

            NetCDFFile::_lkr_writeData(NetCDFData(cell_angles,
                                                  cell_angles.hyperslab()[current_frame],
                                                  double_to_float(frame_buffer->cell_angles)));
        }
    }
}

void AmberRstFile::writeFrame(int frame_idx, const Frame &frame, bool use_parallel)
{
    // create a frame buffer for this frame
    std::shared_ptr<AmberRstFrameBuffer> buffer(new AmberRstFrameBuffer(frame));

    buffer->current_frame = frame_idx;
    buffer->time = frame.time().to(picosecond);

    if (frame.space().isPeriodic())
    {
        if (frame.space().isA<PeriodicBox>())
        {
            buffer->cell_angles[0] = 90;
            buffer->cell_angles[1] = 90;
            buffer->cell_angles[2] = 90;

            const auto dims = frame.space().asA<PeriodicBox>().dimensions();

            buffer->cell_lengths[0] = dims.x();
            buffer->cell_lengths[1] = dims.y();
            buffer->cell_lengths[2] = dims.z();
        }
        else
        {
            TriclinicBox box;

            if (frame.space().isA<TriclinicBox>())
                box = frame.space().asA<TriclinicBox>();
            else
            {
                const auto m = frame.space().boxMatrix();
                box = TriclinicBox(m.column0(), m.column1(), m.column2());
            }

            buffer->cell_lengths[0] = box.vector0().magnitude();
            buffer->cell_lengths[1] = box.vector1().magnitude();
            buffer->cell_lengths[2] = box.vector2().magnitude();

            buffer->cell_angles[0] = box.alpha();
            buffer->cell_angles[1] = box.beta();
            buffer->cell_angles[2] = box.gamma();
        }
    }

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
        auto c_data = buffer->coords.data();

        if (c_data == 0)
            return;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = coords_data[i];

                    c_data[(3*i) + 0] = value.x();
                    c_data[(3*i) + 1] = value.y();
                    c_data[(3*i) + 2] = value.z();
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = coords_data[i];

                c_data[(3 * i) + 0] = value.x();
                c_data[(3 * i) + 1] = value.y();
                c_data[(3 * i) + 2] = value.z();
            }
        }
    };

    auto copy_velocities = [&]()
    {
        if (not frame.hasVelocities())
            return;

        const auto vels_data = frame.velocities().constData();
        auto v_data = buffer->vels.data();

        if (v_data == 0)
            return;

        const double units = 1.0 / (angstrom / (20.455 * picosecond)).value();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = vels_data[i];

                    v_data[(3*i) + 0] = value.x().value() * units;
                    v_data[(3*i) + 1] = value.y().value() * units;
                    v_data[(3*i) + 2] = value.z().value() * units;
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = vels_data[i];

                v_data[(3 * i) + 0] = value.x().value() * units;
                v_data[(3 * i) + 1] = value.y().value() * units;
                v_data[(3 * i) + 2] = value.z().value() * units;
            }
        }
    };

    auto copy_forces = [&]()
    {
        if (not frame.hasForces())
            return;

        const auto frcs_data = frame.forces().constData();
        auto f_data = buffer->frcs.data();

        if (f_data == 0)
            return;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = frcs_data[i];

                    f_data[(3*i) + 0] = value.x().to(kcal / angstrom);
                    f_data[(3*i) + 1] = value.y().to(kcal / angstrom);
                    f_data[(3*i) + 2] = value.z().to(kcal / angstrom);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = frcs_data[i];

                f_data[(3 * i) + 0] = value.x().to(kcal / angstrom);
                f_data[(3 * i) + 1] = value.y().to(kcal / angstrom);
                f_data[(3 * i) + 2] = value.z().to(kcal / angstrom);
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

    QMutexLocker lkr(NetCDFFile::globalMutex());

    // replace the current buffer with this new buffer
    frame_buffer.reset();
    frame_buffer = buffer;

    // now write this to the file
    this->_lkr_writeBufferToFile();
}

NetCDFHyperSlab get_slab(const NetCDFHyperSlab &slab, int frame)
{
    if (slab.nDimensions() == 3)
    {
        // this is a trajectory file - we want the ith frame
        return slab[frame];
    }
    else if (slab.nDimensions() == 2)
    {
        // we can only really be asking for the 0th frame
        if (frame != 0)
            throw SireError::invalid_index(QObject::tr(
                                               "Cannot read frame %1 as the number of frames is 0.")
                                               .arg(frame),
                                           CODELOC);
        return slab;
    }
    else
    {
        throw SireIO::parse_error(QObject::tr(
                                      "Unexpected number of dimensions (%1). Only expecting 2 or 3.")
                                      .arg(slab.nDimensions()),
                                  CODELOC);
    }

    return slab;
}

NetCDFHyperSlab get_time_slab(const NetCDFHyperSlab &slab, int frame)
{
    if (slab.nDimensions() == 1)
    {
        // this is a trajectory file - we want the ith frame
        return slab[frame];
    }
    else if (slab.nDimensions() == 0)
    {
        // we can only really be asking for the 0th frame
        if (frame != 0)
            throw SireError::invalid_index(QObject::tr(
                                               "Cannot read frame %1 as the number of frames is 0.")
                                               .arg(frame),
                                           CODELOC);
        return slab;
    }
    else
    {
        throw SireIO::parse_error(QObject::tr(
                                      "Unexpected number of dimensions (%1). Only expecting 0 or 1.")
                                      .arg(slab.nDimensions()),
                                  CODELOC);
    }

    return slab;
}

NetCDFHyperSlab get_cell_slab(const NetCDFHyperSlab &slab, int frame)
{
    if (slab.nDimensions() == 2)
    {
        // this is a trajectory file - we want the ith frame
        return slab[frame];
    }
    else if (slab.nDimensions() == 1)
    {
        // we can only really be asking for the 0th frame
        if (frame != 0)
            throw SireError::invalid_index(QObject::tr(
                                               "Cannot read frame %1 as the number of frames is 0.")
                                               .arg(frame),
                                           CODELOC);
        return slab;
    }
    else
    {
        throw SireIO::parse_error(QObject::tr(
                                      "Unexpected number of dimensions (%1). Only expecting 1 or 2.")
                                      .arg(slab.nDimensions()),
                                  CODELOC);
    }

    return slab;
}

void AmberRstFile::_lkr_readFrameIntoBuffer(int frame)
{
    if (frame < 0 or frame >= nframes)
    {
        throw SireError::invalid_index(QObject::tr(
                                           "Cannot read frame %1 as the number of frames is %2.")
                                           .arg(frame)
                                           .arg(nframes),
                                       CODELOC);
    }

    if (frame_buffer.get() == 0)
    {
        frame_buffer.reset(new AmberRstFrameBuffer());
    }
    else if (frame_buffer->current_frame == frame)
    {
        // we've already read this frame
        return;
    }

    const auto dimensions = this->_lkr_getVariablesInfo();

    const auto time = dimensions.value("time");

    if (time.isNull())
    {
        frame_buffer->time = 0;
    }
    else
    {
        const auto slab = get_time_slab(time.hyperslab(), frame);
        auto time_data = this->_lkr_read(time, slab);

        if (time.typeSize() == 4)
        {
            const auto t = time_data.toFloatArray();

            if (t.count() != 1)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "Unexpected number of time values (%1). Should be 1.")
                                              .arg(t.count()),
                                          CODELOC);
            }

            frame_buffer->time = t[0];
        }
        else if (time.typeSize() == 8)
        {
            const auto t = time_data.toDoubleArray();

            if (t.count() != 1)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "Unexpected number of time values (%1). Should be 1.")
                                              .arg(t.count()),
                                          CODELOC);
            }

            frame_buffer->time = t[0];
        }
        else
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Invalid data size for time: %1")
                                          .arg(time.typeSize()),
                                      CODELOC);
        }
    }

    const auto cell_lengths = dimensions.value("cell_lengths");

    if (cell_lengths.isNull())
    {
        frame_buffer->cell_lengths = QVector<double>(3, 0.0);
    }
    else
    {
        const auto slab = get_cell_slab(cell_lengths.hyperslab(), frame);
        auto cell_lengths_data = this->_lkr_read(cell_lengths, slab);

        if (cell_lengths.typeSize() == 4)
        {
            const auto c = cell_lengths_data.toFloatArray();

            if (c.count() != 3)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "Unexpected number of cell_lengths values (%1). Should be 3.")
                                              .arg(c.count()),
                                          CODELOC);
            }

            frame_buffer->cell_lengths = float_to_double(c);
        }
        else if (cell_lengths.typeSize() == 8)
        {
            const auto c = cell_lengths_data.toDoubleArray();

            if (c.count() != 3)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "Unexpected number of cell_lengths values (%1). Should be 3.")
                                              .arg(c.count()),
                                          CODELOC);
            }

            frame_buffer->cell_lengths = c;
        }
        else
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Invalid data size for cell_lengths: %1")
                                          .arg(cell_lengths.typeSize()),
                                      CODELOC);
        }
    }

    const auto cell_angles = dimensions.value("cell_angles");

    if (cell_angles.isNull())
    {
        frame_buffer->cell_angles = QVector<double>(3, 0.0);
    }
    else
    {
        const auto slab = get_cell_slab(cell_angles.hyperslab(), frame);
        auto cell_angles_data = this->_lkr_read(cell_angles, slab);

        if (cell_angles.typeSize() == 4)
        {
            const auto c = cell_angles_data.toFloatArray();

            if (c.count() != 3)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "Unexpected number of cell_angles values (%1). Should be 3.")
                                              .arg(c.count()),
                                          CODELOC);
            }

            frame_buffer->cell_angles = float_to_double(c);
        }
        else if (cell_angles.typeSize() == 8)
        {
            const auto c = cell_angles_data.toDoubleArray();

            if (c.count() != 3)
            {
                throw SireIO::parse_error(QObject::tr(
                                              "Unexpected number of cell_angles values (%1). Should be 3.")
                                              .arg(c.count()),
                                          CODELOC);
            }

            frame_buffer->cell_angles = c;
        }
        else
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Invalid data size for cell_angles: %1")
                                          .arg(cell_angles.typeSize()),
                                      CODELOC);
        }
    }

    const auto coordinates = dimensions.value("coordinates");

    if (coordinates.isNull())
    {
        frame_buffer->coords.clear();
    }
    else
    {
        const auto slab = get_slab(coordinates.hyperslab(), frame);
        auto coordinates_data = this->_lkr_read(coordinates, slab);

        if (coordinates.typeSize() == 4)
        {
            // this is float data
            frame_buffer->coords = float_to_double(coordinates_data.toFloatArray());
        }
        else if (coordinates.typeSize() == 8)
        {
            // this is double data
            frame_buffer->coords = coordinates_data.toDoubleArray();
        }
        else
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Invalid data size for coordinates: %1")
                                          .arg(coordinates.typeSize()),
                                      CODELOC);
        }
    }

    const auto velocities = dimensions.value("velocities");

    if (velocities.isNull())
    {
        frame_buffer->vels.clear();
    }
    else
    {
        const auto slab = get_slab(velocities.hyperslab(), frame);
        auto velocities_data = this->_lkr_read(velocities, slab);

        if (velocities.typeSize() == 4)
        {
            // this is float data
            frame_buffer->vels = float_to_double(velocities_data.toFloatArray());
        }
        else if (velocities.typeSize() == 8)
        {
            // this is double data
            frame_buffer->vels = velocities_data.toDoubleArray();
        }
        else
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Invalid data size for velocities: %1")
                                          .arg(velocities.typeSize()),
                                      CODELOC);
        }
    }

    const auto forces = dimensions.value("forces");

    if (forces.isNull())
    {
        frame_buffer->frcs.clear();
    }
    else
    {
        const auto slab = get_slab(forces.hyperslab(), frame);
        auto forces_data = this->_lkr_read(forces, slab);

        if (forces.typeSize() == 4)
        {
            // this is float data
            frame_buffer->frcs = float_to_double(forces_data.toFloatArray());
        }
        else if (forces.typeSize() == 8)
        {
            // this is double data
            frame_buffer->frcs = forces_data.toDoubleArray();
        }
        else
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Invalid data size for forces: %1")
                                          .arg(forces.typeSize()),
                                      CODELOC);
        }
    }

    frame_buffer->natoms = 0;

    if (not frame_buffer->coords.isEmpty())
    {
        frame_buffer->natoms = frame_buffer->coords.count() / 3;
    }

    if (not frame_buffer->vels.isEmpty())
    {
        int nats = frame_buffer->vels.count() / 3;

        if (frame_buffer->natoms == 0)
        {
            frame_buffer->natoms = nats;
        }
        else if (frame_buffer->natoms != nats)
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Disagreement in number of atoms between coordinates (%1) "
                                          "and velocities (%2).")
                                          .arg(frame_buffer->natoms)
                                          .arg(nats),
                                      CODELOC);
        }
    }

    if (not frame_buffer->frcs.isEmpty())
    {
        int nats = frame_buffer->frcs.count() / 3;

        if (frame_buffer->natoms == 0)
        {
            frame_buffer->natoms = nats;
        }
        else if (frame_buffer->natoms != nats)
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Disagreement in number of atoms between coordinates or velocities (%1) "
                                          "and forces (%2).")
                                          .arg(frame_buffer->natoms)
                                          .arg(nats),
                                      CODELOC);
        }
    }

    frame_buffer->current_frame = frame;
}

Frame AmberRstFile::readFrame(int i, bool use_parallel) const
{
    AmberRstFile *nonconst_this = const_cast<AmberRstFile *>(this);

    QMutexLocker lkr(NetCDFFile::globalMutex());

    // load the frame into the buffer
    nonconst_this->_lkr_readFrameIntoBuffer(i);

    if (frame_buffer.get() == 0)
        // no frame has been loaded?
        return Frame();

    // copy it out from the buffer to local storage
    const auto natoms = frame_buffer->natoms;

    const Vector cell_lengths(frame_buffer->cell_lengths[0],
                              frame_buffer->cell_lengths[1],
                              frame_buffer->cell_lengths[2]);

    const Vector cell_angles(frame_buffer->cell_angles[0],
                             frame_buffer->cell_angles[1],
                             frame_buffer->cell_angles[2]);

    const auto time = double(frame_buffer->time) * picosecond;

    // and also the coords / vels / forces data
    QVector<double> coords = frame_buffer->coords;
    QVector<double> vels = frame_buffer->vels;
    QVector<double> frcs = frame_buffer->frcs;

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

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    c_data[i] = Vector(coords_data[(3 * i) + 0],
                                       coords_data[(3 * i) + 1],
                                       coords_data[(3 * i) + 2]);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                c_data[i] = Vector(coords_data[(3 * i) + 0],
                                   coords_data[(3 * i) + 1],
                                   coords_data[(3 * i) + 2]);
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

        // velocity is Angstroms per 1/20.455 ps
        const auto vel_unit = (1.0 / 20.455) * angstrom / picosecond;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    v_data[i] = Velocity3D(vels_data[(3 * i) + 0] * vel_unit,
                                           vels_data[(3 * i) + 1] * vel_unit,
                                           vels_data[(3 * i) + 2] * vel_unit);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                v_data[i] = Velocity3D(vels_data[(3 * i) + 0] * vel_unit,
                                       vels_data[(3 * i) + 1] * vel_unit,
                                       vels_data[(3 * i) + 2] * vel_unit);
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

        const auto units_to_internal = (kcal / angstrom);

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    f_data[i] = Force3D(frcs_data[(3 * i) + 0] * units_to_internal,
                                        frcs_data[(3 * i) + 1] * units_to_internal,
                                        frcs_data[(3 * i) + 2] * units_to_internal);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                f_data[i] = Force3D(frcs_data[(3 * i) + 0] * units_to_internal,
                                    frcs_data[(3 * i) + 1] * units_to_internal,
                                    frcs_data[(3 * i) + 2] * units_to_internal);
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

    SpacePtr space;

    if (cell_lengths.isZero() or cell_angles.isZero())
    {
        space = Cartesian();
    }
    else if (cell_angles == Vector(90, 90, 90))
    {
        space = PeriodicBox(cell_lengths);
    }
    else
    {
        space = TriclinicBox(cell_lengths.x(), cell_lengths.y(), cell_lengths.z(),
                             cell_angles.x() * degrees,
                             cell_angles.y() * degrees,
                             cell_angles.z() * degrees);
    }

    return Frame(c, v, f, space.read(), time);
}

bool AmberRstFile::open(QIODevice::OpenMode mode)
{
    QMutexLocker lkr(NetCDFFile::globalMutex());

    if (not this->_lkr_open(mode, true, false))
        return false;

    if (mode != QIODevice::ReadOnly)
    {
        // we only want to write, and so we don't need to check the file
        return true;
    }

    // now get the attributes from the file
    QString conventions;

    try
    {
        conventions = this->_lkr_getStringAttribute("Conventions");
    }
    catch (...)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "Cannot read the AmberRst file '%1' as it is missing the required 'Conventions' "
                                      "attribute.")
                                      .arg(this->filename()),
                                  CODELOC);
    }

    double convention_version;

    try
    {
        bool version_ok;

        convention_version = this->_lkr_getStringAttribute("ConventionVersion").toDouble(&version_ok);

        if (not version_ok)
        {
            // have to assume it is version 1...
            convention_version = 1;
        }
    }
    catch (...)
    {
        // have to assume it is version 1...
        convention_version = 1;
    }

    if (not conventions.contains("AMBER"))
    {
        // this is not an Amber RST file
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "Cannot read the AmberRst file '%1' as it doesn't match the required format. "
                                      "It is missing the 'AMBER' flag in the 'Conventions' attribute.")
                                      .arg(this->filename()),
                                  CODELOC);
    }

    is_restart = conventions.contains("AMBERRESTART");

    if (convention_version != 1)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "Cannot read the AmberRst file '%1' as it is written in version %1 of the format, "
                                      "and we only support version 1.0")
                                      .arg(convention_version),
                                  CODELOC);
    }

    try
    {
        application = this->_lkr_getStringAttribute("application");
    }
    catch (...)
    {
    }

    try
    {
        program = this->_lkr_getStringAttribute("program");
    }
    catch (...)
    {
    }

    try
    {
        program_version = this->_lkr_getStringAttribute("programVersion");
    }
    catch (...)
    {
    }

    try
    {
        title = this->_lkr_getStringAttribute("title");
    }
    catch (...)
    {
    }

    const auto dims = this->_lkr_getDimensions();

    nframes = dims.value("frame", 1);
    natoms = dims.value("atom", 0);

    int spatial = dims.value("spatial", 3);

    if (spatial != 3)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "We can only read AmberRst files with 3 spatial dimensions. The number "
                                      "of dimensions in '%1' is %2.")
                                      .arg(this->filename())
                                      .arg(spatial),
                                  CODELOC);
    }

    return true;
}

//////////
////////// Implementation of AmberRst
//////////

static const RegisterMetaType<AmberRst> r_rst;
const RegisterParser<AmberRst> register_rst;

QDataStream &operator<<(QDataStream &ds, const AmberRst &rst)
{
    writeHeader(ds, r_rst, 3);

    SharedDataStream sds(ds);

    sds << rst.current_frame << rst.parse_warnings
        << rst.nframes << rst.frame_idx
        << static_cast<const MoleculeParser &>(rst);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AmberRst &rst)
{
    VersionID v = readHeader(ds, r_rst);

    if (v == 3)
    {
        SharedDataStream sds(ds);

        sds >> rst.current_frame >> rst.parse_warnings >> rst.nframes >> rst.frame_idx >> static_cast<MoleculeParser &>(rst);
    }
    else
        throw version_error(v, "3", r_rst, CODELOC);

    return ds;
}

/** Constructor */
AmberRst::AmberRst()
    : ConcreteProperty<AmberRst, MoleculeParser>(),
      nframes(0), frame_idx(0)
{
}

/** Return the format name that is used to identify this file format within Sire */
QString AmberRst::formatName() const
{
    return "RST";
}

/** Return the suffixes that AmberRst files will typically have */
QStringList AmberRst::formatSuffix() const
{
    static const QStringList suffixes = {"rst", "crd", "trj", "traj"};
    return suffixes;
}

/** Return a description of the file format */
QString AmberRst::formatDescription() const
{
    return QObject::tr("Amber coordinate/velocity binary (netcdf) restart/trajectory files "
                       "supported since Amber 9, now default since Amber 16.");
}

/** This is not a text file */
bool AmberRst::isTextFile() const
{
    return false;
}

/** Open the file and read in all the metadata */
void AmberRst::parse()
{
    // the netcdf conventions for Amber restart/trajectory files are
    // given here - http://ambermd.org/netcdf/nctraj.xhtml

    // NOTE THAT YOU CANNOT READ A NETCDF FILE IN PARALLEL
    // THE NETCDF LIBRARY IS NOT THREAD SAFE

    f.reset(new AmberRstFile(this->filename()));

    try
    {
        if (not f->open(QIODevice::ReadOnly))
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Failed to open AmberRst file %1")
                                          .arg(this->filename()),
                                      CODELOC);
        }

        nframes = f->nFrames();

        if (nframes == 0)
        {
            this->setScore(0);
            f.reset();
            return;
        }

        current_frame = f->readFrame(0, this->usesParallel());
        frame_idx = 0;

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
AmberRst::AmberRst(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<AmberRst, MoleculeParser>(map),
      nframes(0), frame_idx(0)
{
    // this gets the absolute file path
    this->setFilename(filename);
    this->parse();
}

/** Construct by parsing the data in the passed text lines */
AmberRst::AmberRst(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<AmberRst, MoleculeParser>(lines, map)
{
    throw SireIO::parse_error(QObject::tr("You cannot create a binary AmberRst file from a set of text lines!"),
                              CODELOC);
}

/** Construct by extracting the necessary data from the passed System */
AmberRst::AmberRst(const System &system, const PropertyMap &map)
    : ConcreteProperty<AmberRst, MoleculeParser>(system, map),
      nframes(1), frame_idx(0)
{
    current_frame = MoleculeParser::createFrame(system, map);
}

/** Copy constructor */
AmberRst::AmberRst(const AmberRst &other)
    : ConcreteProperty<AmberRst, MoleculeParser>(other),
      current_frame(other.current_frame), parse_warnings(other.parse_warnings),
      nframes(other.nframes), frame_idx(other.frame_idx), f(other.f)
{
}

/** Destructor */
AmberRst::~AmberRst()
{
}

AmberRst &AmberRst::operator=(const AmberRst &other)
{
    if (this != &other)
    {
        current_frame = other.current_frame;
        parse_warnings = other.parse_warnings;
        nframes = other.nframes;
        frame_idx = other.frame_idx;
        f = other.f;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool AmberRst::operator==(const AmberRst &other) const
{
    return MoleculeParser::operator==(other);
}

bool AmberRst::operator!=(const AmberRst &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *AmberRst::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberRst>());
}

const char *AmberRst::what() const
{
    return AmberRst::typeName();
}

bool AmberRst::isFrame() const
{
    return true;
}

int AmberRst::nFrames() const
{
    return nframes;
}

int AmberRst::count() const
{
    return this->nFrames();
}

int AmberRst::size() const
{
    return this->nFrames();
}

Frame AmberRst::getFrame(int frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    if (frame < 0)
        frame = 0;

    if (frame == frame_idx)
        return current_frame;

    if (f.get() == 0)
    {
        throw SireError::file_error(QObject::tr(
                                        "Somehow we don't have access to the underlying AmberRst file?"),
                                    CODELOC);
    }

    return f->readFrame(frame, this->usesParallel());
}

AmberRst AmberRst::operator[](int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    AmberRst ret(*this);

    ret.current_frame = this->getFrame(i);
    ret.frame_idx = i;

    return ret;
}

QString AmberRst::toString() const
{
    if (this->nAtoms() == 0)
    {
        return QObject::tr("AmberRst::null");
    }
    else
    {
        return QObject::tr("AmberRst( nAtoms() = %1, nFrames() = %2 )")
            .arg(this->nAtoms())
            .arg(this->nFrames());
    }
}

/** Parse from the passed file */
AmberRst AmberRst::parse(const QString &filename)
{
    return AmberRst(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void AmberRst::addToSystem(System &system, const PropertyMap &map) const
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

/** Return the number of atoms whose coordinates are contained in this restart file */
int AmberRst::nAtoms() const
{
    return current_frame.nAtoms();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr AmberRst::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberRst(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr AmberRst::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberRst(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr AmberRst::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberRst(system, map));
}

/** Write this binary file 'filename' */
QStringList AmberRst::writeToFile(const QString &filename) const
{
    if (this->nFrames() == 0 or this->nAtoms() == 0)
        return QStringList();

    auto gil = SireBase::release_gil();

    createDirectoryForFile(filename);

    AmberRstFile outfile(filename);

    if (not outfile.open(QIODevice::WriteOnly))
        throw SireError::file_error(QObject::tr(
                                        "Could not open %1 to write the AmberRst file.")
                                        .arg(filename),
                                    CODELOC);

    if (this->writingTrajectory())
    {
        const auto frames = this->framesToWrite();

        if (frames.isEmpty())
            return QStringList();

        Frame first_frame = this->createFrame(frames[0]);

        outfile.writeHeader(this->saveTitle(), frames.count(),
                            first_frame.nAtoms(),
                            first_frame.hasCoordinates(),
                            first_frame.hasVelocities(),
                            first_frame.hasForces(),
                            first_frame.space().isPeriodic());

        ProgressBar bar("Save RST", frames.count());
        bar.setSpeedUnit("frames / s");

        bar = bar.enter();

        outfile.writeFrame(0, first_frame, usesParallel());
        first_frame = Frame();
        bar.setProgress(1);

        for (int i = 1; i < frames.count(); ++i)
        {
            const auto frame = this->createFrame(frames[i]);
            outfile.writeFrame(i, frame, usesParallel());
            bar.setProgress(i + 1);
        }

        bar.success();
    }
    else
    {
        outfile.writeHeader(this->saveTitle(), 1, current_frame.nAtoms(),
                            current_frame.hasCoordinates(),
                            current_frame.hasVelocities(),
                            current_frame.hasForces(),
                            current_frame.space().isPeriodic());

        outfile.writeFrame(0, current_frame, usesParallel());
    }

    outfile.close();

    return QStringList(filename);
}
