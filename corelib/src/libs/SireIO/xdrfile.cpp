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

#include "SireError/errors.h"

#include "third_party/xdrfile.h"

using namespace SireIO;

XDRFile::XDRFile() : boost::noncopyable(), f(0)
{
}

XDRFile::XDRFile(const QString &filename)
    : boost::noncopyable(), fname(filename), f(0)
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
    }
}
