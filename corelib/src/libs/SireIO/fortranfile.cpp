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

#include "fortranfile.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include <QDataStream>
#include <QDebug>
#include <QFile>
#include <QFileInfo>
#include <QtEndian>
#include <QMutex>

using namespace SireIO;

//////
////// Implementatin of FortranFile
//////

class FortranFileHandle
{
public:
    FortranFileHandle(const QString &name, QIODevice::OpenMode mode = QIODevice::ReadOnly)
    {
        filename = name;

        f.reset(new QFile(filename));

        if (not f->open(QIODevice::ReadOnly))
        {
            throw SireError::io_error(
                QObject::tr("Could not open file %1. Please check it exists and is readable.").arg(filename), CODELOC);
        }

        ds.reset(new QDataStream(f.get()));
        current_pos = 0;
    }

    ~FortranFileHandle()
    {
    }

    qint64 _lkr_skip_to(qint64 position)
    {
        // re-read from the beginning
        f.reset(new QFile(filename));

        if (not f->open(QIODevice::ReadOnly))
        {
            throw SireError::io_error(
                QObject::tr("Could not open file %1. Please check it exists and is readable.").arg(filename), CODELOC);
        }

        ds.reset(new QDataStream(f.get()));
        current_pos = ds->skipRawData(position);
        return current_pos;

        /*
        if (position > current_pos)
        {
            current_pos += ds->skipRawData(position - current_pos);
        }
        else if (position < current_pos)
        {
            ds.reset(new QDataStream(f.get()));
            current_pos = ds->skipRawData(position);
        }*/

        return current_pos;
    }

    qint64 _lkr_read(char *data, qint64 size)
    {
        if (size > 0)
        {
            int read = ds->readRawData(data, size);
            current_pos += read;
            return read;
        }
        else
            return 0;
    }

    void _lkr_write(const char *data, qint64 size)
    {
        if (size > 0)
        {
            int written = ds->writeRawData(data, size);

            if (written != size)
                throw SireError::io_error(QObject::tr(
                                              "Failed to write %1 bytes to the file. IO Error")
                                              .arg(size),
                                          CODELOC);

            current_pos += written;
        }
    }

    std::shared_ptr<QFile> f;
    std::shared_ptr<QDataStream> ds;

    QString filename;
    QMutex mutex;
    qint64 current_pos;
};

FortranFile::FortranFile() : int_size(4), is_little_endian(true)
{
}

void FortranFile::write(const FortranRecord &record)
{
    if (record.size() <= 0)
        return;

    if (f.get() == 0)
        throw SireError::io_error(
            QObject::tr("Problem writing file '%1'. File pointer has been closed?").arg(abs_filename), CODELOC);

    QMutexLocker lkr(&(f->mutex));

    // write the record size
    if (int_size == 4)
    {
        qint32 sz;

        if (is_little_endian)
        {
            sz = qToLittleEndian<qint32>(record.size());
        }
        else
        {
            sz = qToBigEndian<qint32>(record.size());
        }

        f->_lkr_write(reinterpret_cast<char *>(&sz), int_size);
    }
    else
    {
        qint64 sz;

        if (is_little_endian)
        {
            sz = qToLittleEndian<qint64>(record.size());
        }
        else
        {
            sz = qToBigEndian<qint64>(record.size());
        }

        f->_lkr_write(reinterpret_cast<char *>(&sz), int_size);
    }

    // now write the record
    f->_lkr_write(record.constData(), record.size());
}

bool FortranFile::try_read()
{
    record_pointers.clear();
    record_sizes.clear();

    QFile file(abs_filename);

    if (not file.open(QIODevice::ReadOnly))
    {
        throw SireError::io_error(
            QObject::tr("Could not open file %1. Please check it exists and is readable.").arg(abs_filename), CODELOC);
    }

    QDataStream ds(&file);

    QByteArray start_buffer(int_size, 0);
    QByteArray end_buffer(int_size, 0);

    qint64 read_count = 0;

    // each fortran record starts and ends with an integer that
    // gives the size in bytes of the record. We will now scan
    // through the file, and if we can consistently get all of the
    // sizes, then we much have the right int_size and endianness
    while (not ds.atEnd())
    {
        int read_size = ds.readRawData(start_buffer.data(), int_size);

        if (read_size != int_size)
        {
            return false;
        }

        qint64 start_size;

        if (int_size == 4)
        {
            if (is_little_endian)
                start_size = qFromLittleEndian<qint32>(start_buffer.data());
            else
                start_size = qFromBigEndian<qint32>(start_buffer.data());
        }
        else
        {
            if (is_little_endian)
                start_size = qFromLittleEndian<qint64>(start_buffer.data());
            else
                start_size = qFromBigEndian<qint64>(start_buffer.data());
        }

        read_count += int_size;

        read_size = ds.skipRawData(start_size);

        if (read_size != start_size)
        {
            // could not read this much!
            return false;
        }

        read_size = ds.readRawData(end_buffer.data(), int_size);

        if (read_size != int_size)
        {
            return false;
        }

        qint64 end_size;

        if (int_size == 4)
        {
            if (is_little_endian)
                end_size = qFromLittleEndian<qint32>(end_buffer.data());
            else
                end_size = qFromBigEndian<qint32>(end_buffer.data());
        }
        else
        {
            if (is_little_endian)
                end_size = qFromLittleEndian<qint64>(end_buffer.data());
            else
                end_size = qFromBigEndian<qint64>(end_buffer.data());
        }

        if (start_size != end_size)
        {
            // disagreement - cannot be a valid record
            return false;
        }

        record_pointers.append(read_count);
        record_sizes.append(start_size);

        read_count += start_size + int_size;
    }

    file.close();

    // we can read it - so return a handle to the file
    f.reset(new FortranFileHandle(abs_filename));

    return true;
}

FortranFile::FortranFile(const QString &filename,
                         QIODevice::OpenMode mode)
    : int_size(4), is_little_endian(true)
{
    abs_filename = QFileInfo(filename).absoluteFilePath();

    if (mode == QIODevice::ReadOnly)
    {
        // try to read using 4 byte header and native endian
        int_size = 4;
        is_little_endian = true;

        if (try_read())
            return;

        // try to read using 8 byte header and native endian
        int_size = 8;
        is_little_endian = true;

        if (try_read())
            return;

        // try to read using 4 byte header and swapped endian
        int_size = 4;
        is_little_endian = false;

        if (try_read())
            return;

        // try to read using 8 byte header and swapped endian
        int_size = 8;
        is_little_endian = false;

        if (try_read())
            return;

        throw SireError::io_error(QObject::tr("Could not read a consistent set of records from %1. "
                                              "It could not be read as a record-based unformatted "
                                              "Fortran binary file.")
                                      .arg(abs_filename),
                                  CODELOC);
    }
    else
    {
        // open the file in writing mode
        f.reset(new FortranFileHandle(abs_filename, mode));
    }
}

FortranFile::FortranFile(const FortranFile &other)
    : abs_filename(other.abs_filename), record_pointers(other.record_pointers), record_sizes(other.record_sizes),
      f(other.f), int_size(other.int_size), is_little_endian(other.is_little_endian)
{
}

FortranFile::~FortranFile()
{
}

FortranFile &FortranFile::operator=(const FortranFile &other)
{
    if (this != &other)
    {
        abs_filename = other.abs_filename;
        record_pointers = other.record_pointers;
        record_sizes = other.record_sizes;
        f = other.f;
        int_size = other.int_size;
        is_little_endian = other.is_little_endian;
    }

    return *this;
}

int FortranFile::nRecords() const
{
    return record_pointers.count();
}

void FortranFile::write(const QByteArray &values)
{
    this->write(FortranRecord(values, is_little_endian));
}

void FortranFile::write(double value)
{
    this->write(QVector<double>(value, 1));
}

void FortranFile::write(float value)
{
    this->write(QVector<float>(value, 1));
}

void FortranFile::write(qint32 value)
{
    this->write(QVector<qint32>(value, 1));
}

void FortranFile::write(qint64 value)
{
    this->write(QVector<qint64>(value, 1));
}

void FortranFile::write(const QVector<double> &values)
{
    int n = values.count();

    if (n <= 0)
        return;

    if (sizeof(double) != 8)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-64bit doubles..."),
                                         CODELOC);

    QByteArray data(8 * n, char('\0'));

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 8);
    }
    else
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian doubles..."),
                                         CODELOC);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian doubles..."),
                                         CODELOC);
    }
    else
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 8);
    }
#endif

    this->write(FortranRecord(data, is_little_endian));
}

void FortranFile::write(const QVector<float> &values)
{
    int n = values.count();

    if (n <= 0)
        return;

    if (sizeof(float) != 4)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-32bit floats..."),
                                         CODELOC);

    QByteArray data(4 * n, char('\0'));

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 4);
    }
    else
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian floats..."),
                                         CODELOC);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian floats..."),
                                         CODELOC);
    }
    else
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 4);
    }
#endif

    this->write(FortranRecord(data, is_little_endian));
}

void FortranFile::write(const QVector<qint32> &values)
{
    int n = values.count();

    if (n <= 0)
        return;

    QByteArray data(4 * n, char('\0'));

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 4);
    }
    else
    {
        // need to reverse the data
        qint32 *int_data = reinterpret_cast<qint32 *>(data.data());

        for (int i = 0; i < n; ++i)
        {
            int_data[i] = qToBigEndian<qint32>(values[i]);
        }

        memcpy(int_data, values.constData(), n * 4);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        qint32 *int_data = reinterpret_cast<qint32 *>(data.data());

        for (int i = 0; i < n; ++i)
        {
            int_data[i] = qToLittleEndian<qint32>(values[i]);
        }

        memcpy(int_data, values.constData(), n * 4);
    }
    else
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 4);
    }
#endif

    this->write(FortranRecord(data, is_little_endian));
}

void FortranFile::write(const QVector<qint64> &values)
{
    int n = values.count();

    if (n <= 0)
        return;

    QByteArray data(8 * n, char('\0'));

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 8);
    }
    else
    {
        // need to reverse the data
        qint64 *int_data = reinterpret_cast<qint64 *>(data.data());

        for (int i = 0; i < n; ++i)
        {
            int_data[i] = qToBigEndian<qint64>(values[i]);
        }

        memcpy(int_data, values.constData(), n * 8);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        qint64 *int_data = reinterpret_cast<qint64 *>(data.data());

        for (int i = 0; i < n; ++i)
        {
            int_data[i] = qToLittleEndian<qint64>(values[i]);
        }

        memcpy(int_data, values.constData(), n * 8);
    }
    else
    {
        // just copy the data
        memcpy(data.data(), values.constData(), n * 8);
    }
#endif

    this->write(FortranRecord(data, is_little_endian));
}

FortranRecord FortranFile::operator[](int i) const
{
    i = SireID::Index(i).map(this->nRecords());

    if (f.get() == 0)
    {
        throw SireError::io_error(
            QObject::tr("Problem opening file '%1'. File pointer has been closed?").arg(abs_filename), CODELOC);
    }

    FortranFileHandle *handle = const_cast<FortranFileHandle *>(f.get());

    QMutexLocker lkr(&(handle->mutex));

    auto pointer = this->record_pointers[i];

    auto skipped = handle->_lkr_skip_to(pointer);

    if (pointer != skipped)
    {
        throw SireError::io_error(QObject::tr("Problem reading record %1. Needed to skip %2 bytes, but could "
                                              "only skip %3. Is the file corrupted?")
                                      .arg(abs_filename)
                                      .arg(pointer)
                                      .arg(skipped),
                                  CODELOC);
    }

    qint64 size = this->record_sizes[i];

    QByteArray array(size, 0);

    int read = handle->_lkr_read(array.data(), size);

    if (size != read)
    {
        throw SireError::io_error(QObject::tr("Problem reading record %1. Needed to read %3 bytes, but could "
                                              "only read %3. Is the file corrupted?")
                                      .arg(abs_filename)
                                      .arg(size)
                                      .arg(read),
                                  CODELOC);
    }

    lkr.unlock();

    return FortranRecord(array, is_little_endian);
}

//////
////// Implementatin of FortranRecord
//////

FortranRecord::FortranRecord() : cursor(0)
{
    is_little_endian = true;
}

FortranRecord::FortranRecord(const QByteArray &d, bool le) : data(d), cursor(0), is_little_endian(le)
{
}

FortranRecord::FortranRecord(const FortranRecord &other)
    : data(other.data), cursor(other.cursor), is_little_endian(other.is_little_endian)
{
}

FortranRecord::~FortranRecord()
{
}

FortranRecord &FortranRecord::operator=(const FortranRecord &other)
{
    if (this != &other)
    {
        data = other.data;
        cursor = other.cursor;
        is_little_endian = other.is_little_endian;
    }

    return *this;
}

int FortranRecord::size() const
{
    return data.count();
}

const char *FortranRecord::constData() const
{
    return data.constData();
}

void FortranRecord::_assertValid(int n, int size) const
{
    if (cursor + (n * size) > data.count())
    {
        throw SireError::io_error(QObject::tr("Cannot read %1 x %2 bytes of data as only %3 bytes remain!")
                                      .arg(n)
                                      .arg(size)
                                      .arg(data.count() - cursor),
                                  CODELOC);
    }
}

void FortranRecord::_assertPosValid(int pos, int size) const
{
    if (pos < 0)
    {
        throw SireError::io_error(
            QObject::tr("Cannot read %1 bytes at position %2 as position is negative!").arg(size).arg(pos), CODELOC);
    }
    else if (pos + size > data.count())
    {
        throw SireError::io_error(QObject::tr("Cannot read %1 bytes at position %2 as only %3 bytes remain!")
                                      .arg(size)
                                      .arg(pos)
                                      .arg(data.count() - pos),
                                  CODELOC);
    }
}

char FortranRecord::readCharAt(int pos) const
{
    _assertPosValid(pos, 1);
    return data.constData()[pos];
}

double FortranRecord::readFloat64At(int pos) const
{
    _assertPosValid(pos, 8);

    double ret;

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 8);
    }
    else
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian doubles..."),
                                         CODELOC);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian doubles..."),
                                         CODELOC);
    }
    else
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 8);
    }
#endif

    return ret;
}

float FortranRecord::readFloat32At(int pos) const
{
    _assertPosValid(pos, 4);

    float ret;

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 4);
    }
    else
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian doubles..."),
                                         CODELOC);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian doubles..."),
                                         CODELOC);
    }
    else
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 4);
    }
#endif

    return ret;
}

qint32 FortranRecord::readInt32At(int pos) const
{
    _assertPosValid(pos, 4);

    qint32 ret;

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 4);
    }
    else
    {
        ret = qFromBigEndian<qint32>(data.constData() + pos);
    }
#else
    if (is_little_endian)
    {
        ret = qFromLittleEndian<qint32>(data.constData() + pos);
    }
    else
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 4);
    }
#endif

    return ret;
}

qint64 FortranRecord::readInt64At(int pos) const
{
    _assertPosValid(pos, 8);

    qint64 ret;

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 8);
    }
    else
    {
        ret = qFromBigEndian<qint64>(data.constData() + pos);
    }
#else
    if (is_little_endian)
    {
        ret = qFromLittleEndian<qint64>(data.constData() + pos);
    }
    else
    {
        // just copy the data
        memcpy(&ret, data.constData() + pos, 8);
    }
#endif

    return ret;
}

QString FortranRecord::readChar(int n)
{
    if (n <= 0)
        return QString();

    _assertValid(n, 1);

    auto ret = QString::fromUtf8(data.constData() + cursor, n);

    cursor += n;

    return ret;
}

QVector<double> FortranRecord::readFloat64(int n)
{
    if (n <= 0)
        return QVector<double>();

    if (sizeof(double) != 8)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-64bit doubles..."),
                                         CODELOC);

    _assertValid(n, 8);

    QVector<double> ret(n);

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 8);
        cursor += n * 8;
    }
    else
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian doubles..."),
                                         CODELOC);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian doubles..."),
                                         CODELOC);
    }
    else
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 8);
        cursor += n * 8;
    }
#endif

    return ret;
}

QVector<float> FortranRecord::readFloat32(int n)
{
    if (n <= 0)
        return QVector<float>();

    if (sizeof(float) != 4)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-32bit floats..."), CODELOC);

    _assertValid(n, 4);

    QVector<float> ret(n);

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 4);
        cursor += n * 4;
    }
    else
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian floats..."),
                                         CODELOC);
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian floats..."),
                                         CODELOC);
    }
    else
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 4);
        cursor += n * 4;
    }
#endif

    return ret;
}

QVector<qint32> FortranRecord::readInt32(int n)
{
    if (n <= 0)
        return QVector<qint32>();

    _assertValid(n, 4);

    QVector<qint32> ret(n);

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 4);
        cursor += n * 4;
    }
    else
    {
        // need to reverse the data
        for (int i = 0; i < n; ++i)
        {
            ret[i] = qFromBigEndian<qint32>(data.constData() + cursor);
            cursor += 4;
        }
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        for (int i = 0; i < n; ++i)
        {
            ret[i] = qFromLittleEndian<qint32>(data.constData() + cursor);
            cursor += 4;
        }
    }
    else
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 4);
        cursor += n * 4;
    }
#endif

    return ret;
}

QVector<qint64> FortranRecord::readInt64(int n)
{
    if (n <= 0)
        return QVector<qint64>();

    _assertValid(n, 8);

    QVector<qint64> ret(n);

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 8);
        cursor += n * 8;
    }
    else
    {
        // need to reverse the data
        for (int i = 0; i < n; ++i)
        {
            ret[i] = qFromBigEndian<qint64>(data.constData() + cursor);
            cursor += 8;
        }
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        for (int i = 0; i < n; ++i)
        {
            ret[i] = qFromLittleEndian<qint64>(data.constData() + cursor);
            cursor += 8;
        }
    }
    else
    {
        // just copy the data
        memcpy(ret.data(), data.constData() + cursor, n * 8);
        cursor += n * 8;
    }
#endif

    return ret;
}
