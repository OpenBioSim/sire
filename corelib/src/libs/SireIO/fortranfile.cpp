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

#include "SireBase/releasegil.h"
#include "SireBase/progressbar.h"

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
    FortranFileHandle(const QString &filename, QIODevice::OpenMode mode = QIODevice::ReadOnly) : f(filename)
    {
        if (not f.open(mode))
        {
            throw SireError::file_error(f, CODELOC);
        }
    }

    ~FortranFileHandle()
    {
    }

    qint64 _lkr_skip_to(qint64 position)
    {
        if (not f.seek(position))
            throw SireError::io_error(QObject::tr("Could not seek to position %1.").arg(position), CODELOC);

        return f.pos();
    }

    qint64 _lkr_read(char *data, qint64 size)
    {
        if (size > 0)
        {
            int read = f.read(data, size);

            if (read != size)
                throw SireError::io_error(QObject::tr(
                                              "Failed to read %1 bytes from the file. Only %2 bytes were read.")
                                              .arg(size)
                                              .arg(read),
                                          CODELOC);

            return read;
        }
        else
            return 0;
    }

    void _lkr_write(const char *data, qint64 size)
    {
        if (size > 0)
        {
            int written = f.write(data, size);

            if (written <= 0)
                throw SireError::io_error(QObject::tr(
                                              "Failed to write %1 bytes to the file. IO Error")
                                              .arg(size),
                                          CODELOC);
            else if (written != size)
                throw SireError::io_error(QObject::tr(
                                              "Failed to write %1 bytes to the file. Only %2 bytes were written.")
                                              .arg(size)
                                              .arg(written),
                                          CODELOC);
        }
    }

    QMutex mutex;
    QFile f;
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

    // now write the size of the record again (it is written
    // at the start and end of the record)
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
}

bool FortranFile::try_read()
{
    auto gil = SireBase::release_gil();

    record_pointers.clear();
    record_sizes.clear();

    QFile file(abs_filename);

    if (not file.open(QIODevice::ReadOnly))
    {
        throw SireError::io_error(
            QObject::tr("Could not open file %1. Please check it exists and is readable.").arg(abs_filename), CODELOC);
    }

    // remember the total file size
    const auto file_size = file.size();

    QDataStream ds(&file);

    QByteArray start_buffer(int_size, 0);
    QByteArray end_buffer(int_size, 0);

    qint64 read_count = 0;

    SireBase::ProgressBar bar(file_size, "Indexing Fortran File");
    bar.setSpeedUnit("bytes / s");

    bar = bar.enter();

    // each fortran record starts and ends with an integer that
    // gives the size in bytes of the record. We will now scan
    // through the file, and if we can consistently get all of the
    // sizes, then we much have the right int_size and endianness
    while (not ds.atEnd())
    {
        int read_size = ds.readRawData(start_buffer.data(), int_size);

        if (read_size != int_size)
        {
            bar.failure();
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

        if (start_size < 0 or start_size > file_size)
        {
            // this is not a fortran file as the sizes don't make sense
            bar.failure();
            return false;
        }

        read_count += int_size;

        read_size = ds.skipRawData(start_size);

        if (read_size != start_size)
        {
            // could not read this much!
            bar.failure();
            return false;
        }

        read_size = ds.readRawData(end_buffer.data(), int_size);

        if (read_size != int_size)
        {
            bar.failure();
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
            bar.failure();
            return false;
        }

        record_pointers.append(read_count);
        record_sizes.append(start_size);

        read_count += start_size + int_size;

        bar.setProgress(read_count);
    }

    bar.success();

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
        // try to read using 4 byte header and little endian
        int_size = 4;
        is_little_endian = true;

        if (try_read())
            return;

        // try to read using 8 byte header and little endian
        int_size = 8;
        is_little_endian = true;

        if (try_read())
            return;

        // try to read using 4 byte header and big endian
        int_size = 4;
        is_little_endian = false;

        if (try_read())
            return;

        // try to read using 8 byte header and big endian
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

        int_size = 8;

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
        is_little_endian = true;
#else
        is_little_endian = false;
#endif
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

bool FortranFile::isLittleEndian() const
{
    return is_little_endian;
}

//////
////// Implementatin of FortranRecord
//////

FortranRecord::FortranRecord() : cursor(0)
{
    is_little_endian = true;
}

FortranRecord::FortranRecord(bool le) : cursor(0), is_little_endian(le)
{
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

void assert_doubles_ok(bool is_little_endian)
{
    if (sizeof(double) != 8)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-64bit doubles..."));

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (not is_little_endian)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian doubles..."));
#else
    if (not is_big_endian)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian doubles..."));
#endif
}

void assert_floats_ok(bool is_little_endian)
{
    if (sizeof(float) != 4)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-32bit floats..."));

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (not is_little_endian)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with big endian floats..."));
#else
    if (not is_big_endian)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with little endian floats..."));
#endif
}

QVector<double> FortranRecord::readFloat64(int n)
{
    if (n <= 0)
        return QVector<double>();

    assert_doubles_ok(is_little_endian);

    _assertValid(n, 8);

    QVector<double> ret(n);

    // just copy the data
    memcpy(ret.data(), data.constData() + cursor, n * 8);
    cursor += n * 8;

    return ret;
}

QVector<float> FortranRecord::readFloat32(int n)
{
    if (n <= 0)
        return QVector<float>();

    assert_floats_ok(is_little_endian);

    _assertValid(n, 4);

    QVector<float> ret(n);

    // just copy the data
    memcpy(ret.data(), data.constData() + cursor, n * 4);
    cursor += n * 4;

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

void FortranRecord::writeFloat64(const QVector<double> &values, int n)
{
    if (n <= 0)
        return;

    if (sizeof(double) != 8)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-64bit doubles..."),
                                         CODELOC);

    QByteArray d(8 * n, char('\0'));

    int n_to_copy = std::min(n, values.count());

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(d.data(), values.constData(), n_to_copy * 8);
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
        memcpy(d.data(), values.constData(), n_to_copy * 8);
    }
#endif

    data.append(d);
    cursor = data.count();
}

void FortranRecord::writeFloat32(const QVector<float> &values, int n)
{
    if (n <= 0)
        return;

    if (sizeof(float) != 4)
        throw SireError::incomplete_code(QObject::tr("Haven't written code to deal with non-32bit floats..."),
                                         CODELOC);

    int n_to_copy = std::min(n, values.count());

    QByteArray d(4 * n, char('\0'));

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(d.data(), values.constData(), n_to_copy * 4);
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
        memcpy(d.data(), values.constData(), n_to_copy * 4);
    }
#endif

    data.append(d);
    cursor = data.count();
}

void FortranRecord::writeInt32(const QVector<qint32> &values, int n)
{
    if (n <= 0)
        return;

    QByteArray d(4 * n, char('\0'));

    int n_to_copy = std::min(n, values.count());

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(d.data(), values.constData(), n_to_copy * 4);
    }
    else
    {
        // need to reverse the data
        qint32 *int_data = reinterpret_cast<qint32 *>(d.data());

        for (int i = 0; i < n_to_copy; ++i)
        {
            int_data[i] = qToBigEndian<qint32>(values[i]);
        }
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        qint32 *int_data = reinterpret_cast<qint32 *>(d.data());

        for (int i = 0; i < n_to_copy; ++i)
        {
            int_data[i] = qToLittleEndian<qint32>(values[i]);
        }
    }
    else
    {
        // just copy the data
        memcpy(d.data(), values.constData(), n_to_copy * 4);
    }
#endif

    data.append(d);
    cursor = data.count();
}

void FortranRecord::writeInt64(const QVector<qint64> &values, int n)
{
    if (n <= 0)
        return;

    QByteArray d(8 * n, char('\0'));
    int n_to_copy = std::min(n, values.count());

#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (is_little_endian)
    {
        // just copy the data
        memcpy(d.data(), values.constData(), n_to_copy * 8);
    }
    else
    {
        // need to reverse the data
        qint64 *int_data = reinterpret_cast<qint64 *>(d.data());

        for (int i = 0; i < n_to_copy; ++i)
        {
            int_data[i] = qToBigEndian<qint64>(values[i]);
        }
    }
#else
    if (is_little_endian)
    {
        // need to reverse the data
        qint64 *int_data = reinterpret_cast<qint64 *>(d.data());

        for (int i = 0; i < n_to_copy; ++i)
        {
            int_data[i] = qToLittleEndian<qint64>(values[i]);
        }
    }
    else
    {
        // just copy the data
        memcpy(d.data(), values.constData(), n_to_copy * 8);
    }
#endif

    data.append(d);
    cursor = data.count();
}

void FortranRecord::writeChar(const QString &text, int n)
{
    if (n <= 0)
        return;

    // we will utf8 encode the strings
    const auto utf8 = text.toUtf8();

    if (utf8.count() >= n)
    {
        // only write the first n characters
        data.append(utf8.constData(), n);
    }
    else
    {
        // make sure that we pad with zeros
        data.append(utf8.constData(), utf8.count());
        data.append(n - utf8.count(), char('\0'));
    }

    cursor = data.count();
}

void FortranRecord::writeFloat64(double value)
{
    assert_doubles_ok(is_little_endian);
    data.append(reinterpret_cast<char *>(&value), 8);
    cursor = data.count();
}

void FortranRecord::writeFloat32(float value)
{
    assert_floats_ok(is_little_endian);
    data.append(reinterpret_cast<char *>(&value), 4);
    cursor = data.count();
}

void FortranRecord::writeInt32(qint32 value)
{
#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (not is_little_endian)
        value = qToBigEndian<qint32>(value);
#else
    if (is_little_endian)
        value = qToLittleEndian<qint32>(value);
#endif

    data.append(reinterpret_cast<char *>(&value), 4);
    cursor = data.count();
}

void FortranRecord::writeInt64(qint64 value)
{
#if Q_BYTE_ORDER == Q_LITTLE_ENDIAN
    if (not is_little_endian)
        value = qToBigEndian<qint64>(value);
#else
    if (is_little_endian)
        value = qToLittleEndian<qint64>(value);
#endif

    data.append(reinterpret_cast<char *>(&value), 8);
    cursor = data.count();
}
