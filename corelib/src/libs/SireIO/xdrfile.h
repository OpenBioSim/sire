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
#include <QMutex>

#include <boost/noncopyable.hpp>

#include <functional>
#include <memory>

SIRE_BEGIN_HEADER

namespace SireMol
{
    class Frame;
}

#ifdef __cplusplus
extern "C"
{
#endif
    typedef struct XDRFILE XDRFILE;
#ifdef __cplusplus
}
#endif

namespace SireIO
{
    /** This class provides an internal interface to XDR files */
    class SIREIO_EXPORT XDRFile : public boost::noncopyable
    {
    public:
        enum FrameDataType
        {
            COORDINATES = 0x0001,
            VELOCITIES = 0x0010,
            FORCES = 0x0100,
            BOX = 0x1000
        };

        XDRFile();

        XDRFile(const QString &filename);

        virtual ~XDRFile();

        virtual bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly);

        void close();

        QString filename() const;

        qint64 size() const;

    protected:
        bool _lkr_open(QIODevice::OpenMode mode = QIODevice::ReadOnly);
        void _lkr_close();
        qint64 _lkr_size() const;

        /** The name of the file */
        QString fname;

        /** Mutex to serialise access to this file */
        QMutex mutex;

        /** Handle to the XDR file */
        XDRFILE *f;

        /** The size of the file, in bytes */
        qint64 sz;
    };

    namespace detail
    {
        class TRRFrameBuffer;
    }

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

        SireMol::Frame readFrame(int i, bool use_parallel = true) const;
        void writeFrame(const SireMol::Frame &frame,
                        bool use_parallel = true);

        int nAtoms() const;
        int nFrames() const;

    private:
        void _lkr_reset();
        void _lkr_reindexFrames();
        void _lkr_readFrameIntoBuffer(int i);
        void _lkr_writeBufferToFile();

        /** The current frame that has been read into the buffer */
        std::shared_ptr<detail::TRRFrameBuffer> frame_buffer;

        /** The number of atoms in the frame - we assume all
         *  frames have the same number of atoms
         */
        qint64 natoms;

        /** The number of frames in the file */
        qint64 nframes;

        /** The seek position of each frame - this is only
         *  used if the frames have different sizes
         */
        QList<std::tuple<qint64, qint32>> seek_frame;

        /** The size, in bytes, of each frame. This is 0
         *  if each frame has a different number of bytes
         */
        qint64 bytes_per_frame;

        /** The data type if all frames are the same */
        qint32 frame_type;
    };

} // namespace SireIO

SIRE_END_HEADER

#endif
