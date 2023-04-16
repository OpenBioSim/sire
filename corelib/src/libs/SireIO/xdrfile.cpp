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

#include <QDebug>
#include <QFileInfo>

#include "SireIO/errors.h"
#include "SireError/errors.h"

#include "third_party/xdrfile.h"
#include "third_party/xdrfile_trr.h"
#include "third_party/trr_seek.h"

using namespace SireIO;

////////
//////// Implementation of XDRFile
////////

XDRFile::XDRFile() : boost::noncopyable(), f(0), sz(0)
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
        this->close();
        f = 0;
    }
}

QString XDRFile::filename() const
{
    return fname;
}

bool XDRFile::open(QIODevice::OpenMode mode)
{
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
    if (f)
    {
        xdrfile_close(f);
        f = 0;
        sz = 0;
    }
}

qint64 XDRFile::size() const
{
    return sz;
}

////////
//////// Implementation of TRRFile
////////

TRRFile::TRRFile()
    : XDRFile(),
      coords_buffer(0), vels_buffer(0), frcs_buffer(0),
      natoms(0), nframes(0)
{
}

TRRFile::TRRFile(const QString &filename)
    : XDRFile(filename),
      coords_buffer(0), vels_buffer(0), frcs_buffer(0),
      natoms(0), nframes(0)
{
}

TRRFile::~TRRFile()
{
    delete[] coords_buffer;
    delete[] vels_buffer;
    delete[] frcs_buffer;
}

void TRRFile::reset()
{
    delete[] coords_buffer;
    delete[] vels_buffer;
    delete[] frcs_buffer;

    coords_buffer = 0;
    vels_buffer = 0;
    frcs_buffer = 0;

    natoms = 0;
    nframes = 0;
}

bool TRRFile::open(QIODevice::OpenMode mode)
{
    if (not XDRFile::open(mode))
        return false;

    qDebug() << "OPENED" << this->size();

    // read the first header to see if this is really a TRR file
    t_trnheader header;

    auto ok = do_trnheader(f, true, &header);

    if (ok != exdrOK)
    {
        this->close();
        throw SireIO::parse_error(QObject::tr(
                                      "The file '%1' is not a valid TRR file. The error message is '%2'")
                                      .arg(this->filename())
                                      .arg(exdr_message[ok]),
                                  CODELOC);
    }

    natoms = header.natoms;

    qDebug() << "natoms =" << natoms;

    if (natoms <= 0)
    {
        // there are no frames either...
        nframes = 0;
        return true;
    }

    // protect against a memory DDOS or file corruption
    if (natoms > 2048 * 2048)
    {
        this->close();
        qint64 natoms_tmp = natoms;
        natoms = 0;
        nframes = 0;
        throw SireError::unsupported(QObject::tr(
                                         "natoms = %1. Reading trajectory files with more than %1 atoms is not supported.")
                                         .arg(natoms_tmp)
                                         .arg(2048 * 2048),
                                     CODELOC);
    }

    // ok - we now seek back to the start of the file and read the
    // frames
    ok = xdr_seek(f, 0, SEEK_SET);

    if (ok != exdrOK)
    {
        this->close();
        natoms = 0;
        nframes = 0;
        throw SireError::file_error(QObject::tr(
                                        "Unable to seek to the start of '%1'")
                                        .arg(this->filename()),
                                    CODELOC);
    }

    // allocate buffers for the coordinates, velocities and forces
    coords_buffer = new rvec[natoms];
    vels_buffer = new rvec[natoms];
    frcs_buffer = new rvec[natoms];

    current_frame = -1;

    while (true)
    {
        matrix box;
        int step;
        float t;
        float lambda;
        int has_prop;
        int nats = natoms;
        ok = read_trr(f, nats, &step, &t, &lambda, box,
                      coords_buffer, vels_buffer, frcs_buffer,
                      &has_prop);

        if (ok != exdrOK)
            break;

        current_frame += 1;

        qDebug() << "READ" << current_frame << step << t;

        qDebug() << "POSITION" << xdr_tell(f);
    }

    nframes = current_frame + 1;

    qDebug() << "number of frames equals" << nframes;

    return true;
}
