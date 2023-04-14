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

    private:
        /** The name of the file */
        QString fname;

        /** Handle to the XDR file */
        XDRFILE *f;
    };

} // namespace SireIO

SIRE_END_HEADER

#endif
