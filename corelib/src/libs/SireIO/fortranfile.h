/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREIO_FORTRANFILE_H
#define SIREIO_FORTRANFILE_H

#include "sireglobal.h"

#include <QByteArray>
#include <QFile>
#include <memory>

SIRE_BEGIN_HEADER

class FortranFileHandle;

namespace SireIO
{

    /** This represents a single Fortran record, from which you
     *  can extract data
     */
    class FortranRecord
    {
    public:
        FortranRecord();
        FortranRecord(bool is_little_endian);
        FortranRecord(const QByteArray &data, bool is_little_endian);
        FortranRecord(const FortranRecord &other);

        ~FortranRecord();

        FortranRecord &operator=(const FortranRecord &other);

        int size() const;
        const char *constData() const;

        QString readChar(int n);
        QVector<double> readFloat64(int n);
        QVector<float> readFloat32(int n);
        QVector<qint32> readInt32(int n);
        QVector<qint64> readInt64(int n);

        char readCharAt(int pos) const;
        double readFloat64At(int pos) const;
        float readFloat32At(int pos) const;
        qint32 readInt32At(int pos) const;
        qint64 readInt64At(int pos) const;

        void writeChar(const QString &text, int n);

        void writeFloat64(const QVector<double> &values, int n);
        void writeFloat32(const QVector<float> &values, int n);
        void writeInt32(const QVector<qint32> &values, int n);
        void writeInt64(const QVector<qint64> &values, int n);

        void writeFloat64(double value);
        void writeFloat32(float value);
        void writeInt32(qint32 value);
        void writeInt64(qint64 value);

    private:
        void _assertPosValid(int pos, int size) const;
        void _assertValid(int n, int size) const;

        QByteArray data;
        int cursor;
        bool is_little_endian;
    };

    /** This class is used to read and write fortran binary
     *  unformatted files (written as record files, not
     *  streaming files).
     *
     *  This automatically detects the endianness of the file
     */
    class FortranFile
    {
    public:
        FortranFile();
        FortranFile(const QString &filename,
                    QIODevice::OpenMode mode = QIODevice::ReadOnly);
        FortranFile(const FortranFile &other);
        ~FortranFile();

        FortranFile &operator=(const FortranFile &other);

        int nRecords() const;

        FortranRecord operator[](int i) const;

        void write(const FortranRecord &record);

        bool isLittleEndian() const;

    private:
        bool try_read();

        QString abs_filename;

        QVector<qint64> record_pointers;
        QVector<qint64> record_sizes;

        std::shared_ptr<FortranFileHandle> f;

        int int_size;

        bool is_little_endian;
    };

} // namespace SireIO

SIRE_END_HEADER

#endif
