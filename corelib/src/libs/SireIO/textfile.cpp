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

#include "textfile.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include <QFileInfo>

using namespace SireIO;

TextFile::TextFile()
    : boost::noncopyable(), f(0), ts(0), sz(0)
{
}

TextFile::TextFile(const QString &filename)
    : boost::noncopyable(), fname(filename),
      f(0), ts(0), sz(0)
{
}

TextFile::~TextFile()
{
    if (ts)
    {
        delete ts;
        ts = 0;
    }

    if (f)
    {
        this->_lkr_close();
        f = 0;
    }
}

bool TextFile::open(QIODevice::OpenMode mode)
{
    QMutexLocker lkr(&mutex);
    return this->_lkr_open(mode);
}

void TextFile::close()
{
    QMutexLocker lkr(&mutex);
    this->_lkr_close();
}

QString TextFile::filename() const
{
    return fname;
}

qint64 TextFile::size() const
{
    QMutexLocker lkr(const_cast<QMutex *>(&mutex));
    return this->_lkr_size();
}

bool TextFile::_lkr_open(QIODevice::OpenMode mode)
{
    this->_lkr_close();

    QFileInfo fullname(this->filename());

    sz = 0;

    if (mode & QIODevice::ReadOnly)
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

        f = new QFile(fullname.absoluteFilePath());

        if (not f->open(mode))
        {
            throw SireError::file_error(*f, CODELOC);
        }

        sz = f->size();
    }
    else if (mode & QIODevice::WriteOnly or mode & QIODevice::Append)
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

        f = new QFile(fullname.absoluteFilePath());

        if (not f->open(mode))
        {
            throw SireError::file_error(*f, CODELOC);
        }
    }

    ts = new QTextStream(f);
    ts->setCodec("UTF-8");
    ts->setAutoDetectUnicode(true);

    return true;
}

void TextFile::_lkr_close()
{
    if (f)
    {
        f->close();
        delete f;
        f = 0;
        sz = 0;

        if (ts)
        {
            delete ts;
            ts = 0;
        }
    }
}

qint64 TextFile::_lkr_size() const
{
    return sz;
}
