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

#ifndef SIREIO_TEXTFILE_H
#define SIREIO_TEXTFILE_H

#include "sireglobal.h"

#include <boost/noncopyable.hpp>

#include <QFile>
#include <QTextStream>
#include <QMutex>

#include <memory>

SIRE_BEGIN_HEADER

namespace SireIO
{
    /** This is a handle for a text file that is designed to
     *  make it fast to access in a line-by-line way. The lines
     *  are cached. There are functions that let you quickly
     *  skip ahead to a specified line
     */
    class TextFile : public boost::noncopyable
    {
    public:
        TextFile();
        TextFile(const QString &filename);

        ~TextFile();

        bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly);

        void close();

        QString filename() const;

        qint64 size() const;

    protected:
        bool _lkr_open(QIODevice::OpenMode mode);
        void _lkr_close();
        qint64 _lkr_size() const;

        /** The name of the file */
        QString fname;

        /** Mutex to serialise access to this file */
        QMutex mutex;

        /** Handle to the text file */
        QFile *f;

        /** The QTextStream used to read the file */
        QTextStream *ts;

        /** The size of the file, in bytes */
        qint64 sz;
    };

}

SIRE_END_HEADER

#endif
