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

#ifndef SIREIO_XDRFILE_H
#define SIREIO_XDRFILE_H

#include "sireglobal.h"

#include <QIODevice>

#include <boost/noncopyable.hpp>

#include <functional>

SIRE_BEGIN_HEADER

/** This is the XDRFILE struct defined in third_party/xdrfile.h */
#ifdef __cplusplus
extern "C"
{
#endif
    typedef struct XDRFILE XDRFILE;
    typedef float rvec[3];
#ifdef __cplusplus
}
#endif

namespace SireIO
{
    /** This class provides an internal interface to XDR files */
    class SIREIO_EXPORT XDRFile : public boost::noncopyable
    {
    public:
        XDRFile();

        XDRFile(const QString &filename);

        ~XDRFile();

        bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly);

        void close();

        QString filename() const;

        qint64 size() const;

    protected:
        /** The name of the file */
        QString fname;

        /** Handle to the XDR file */
        XDRFILE *f;

        /** The size of the file, in bytes */
        qint64 sz;
    };

    /** This class builds on the XDRFile to provide a higher-level
     *  interface for TRR files
     */
    class SIREIO_EXPORT TRRFile : public XDRFile
    {
    public:
        TRRFile();
        TRRFile(const QString &filename);
        ~TRRFile();

        bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly);

        qint64 nAtoms() const;
        qint64 nFrames() const;

    private:
        void reset();

        rvec *coords_buffer;
        rvec *vels_buffer;
        rvec *frcs_buffer;

        /** The number of atoms in the frame - we assume all
         *  frames have the same number of atoms
         */
        qint64 natoms;

        /** The number of frames in the file */
        qint64 nframes;

        /** The index of the current read frame.
         *  This is -1 if no frames have been read
         */
        qint64 current_frame;
    };

} // namespace SireIO

SIRE_END_HEADER

#endif
