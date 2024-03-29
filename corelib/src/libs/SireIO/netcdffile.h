/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef SIREIO_NETCDFFILE_H
#define SIREIO_NETCDFFILE_H

#include "sireglobal.h"

#include <QMutex>
#include <QString>
#include <QStringList>
#include <QVariant>
#include <QIODevice>

#include <boost/noncopyable.hpp>

#include <functional>

SIRE_BEGIN_HEADER

namespace SireIO
{

    /** This class represents a hyperslab in a NetCDF file */
    class SIREIO_EXPORT NetCDFHyperSlab
    {
    public:
        NetCDFHyperSlab();
        NetCDFHyperSlab(const QVector<size_t> &starts,
                        const QVector<size_t> &counts);
        NetCDFHyperSlab(const NetCDFHyperSlab &other);
        ~NetCDFHyperSlab();

        NetCDFHyperSlab &operator=(const NetCDFHyperSlab &other);

        const size_t *starts() const;
        const size_t *counts() const;

        QString toString() const;

        int nDimensions() const;

        bool isNull() const;

        NetCDFHyperSlab operator[](int i) const;

        NetCDFHyperSlab operator()(int i) const;
        NetCDFHyperSlab operator()(int i, int j) const;
        NetCDFHyperSlab operator()(int i, int j, int k) const;

    private:
        QVector<size_t> sts;
        QVector<size_t> cts;
    };

    /** This class provides information about a data variable in
        a NetCDF file

        @author Christopher Woods
    */
    class SIREIO_EXPORT NetCDFDataInfo
    {

        friend class NetCDFFile;

    public:
        template <class T>
        static QString get_nc_type()
        {
            if (std::is_same<T, float>::value)
            {
                return "NC_FLOAT";
            }
            else if (std::is_same<T, double>::value)
            {
                return "NC_DOUBLE";
            }
            else if (std::is_same<T, char>::value)
            {
                return "NC_CHAR";
            }
            else if (std::is_same<T, qint32>::value)
            {
                return "NC_INT";
            }
            else if (std::is_same<T, qint64>::value)
            {
                return "NC_INT64";
            }
            else
            {
                return "unknown";
            }
        }

        NetCDFDataInfo();

        NetCDFDataInfo(const NetCDFDataInfo &other);
        NetCDFDataInfo(const NetCDFDataInfo &other,
                       const NetCDFHyperSlab &slab);

        NetCDFDataInfo(const QString &nc_type,
                       const QString &name,
                       const QStringList &dim_names,
                       const QList<int> &dim_sizes,
                       const QHash<QString, QVariant> &attributes);

        ~NetCDFDataInfo();

        int ID() const;
        QString name() const;
        QString type() const;

        int typeSize() const;
        int dataSize() const;

        QStringList dimensions() const;
        QList<int> dimensionSizes() const;

        int nValues() const;

        int nAttributes() const;

        QStringList attributeNames() const;

        QVariant attribute(const QString &name) const;
        QString attributeType(const QString &name) const;

        QHash<QString, QVariant> attributes() const;
        QHash<QString, QString> attributeTypes() const;

        NetCDFHyperSlab hyperslab() const;

        QString toString() const;

        bool isNull() const;

        void assertNValuesEquals(int nvalues) const;

    protected:
        NetCDFDataInfo(int idnum, QString name, int xtyp, QStringList dim_names, QList<int> dim_sizes,
                       QStringList att_names, QList<int> att_types, QList<QVariant> att_values);

        NetCDFDataInfo(int idnum, QString name, const QString &xtyp, QStringList dim_names, QList<int> dim_sizes,
                       QStringList att_names, QList<int> att_types, QList<QVariant> att_values);

        static QStringList get_attribute_names(const QHash<QString, QVariant> &attributes);
        static QList<int> get_attribute_types(const QHash<QString, QVariant> &attributes);
        static QList<QVariant> get_attribute_values(const QHash<QString, QVariant> &attributes);

        /** The name of the variable */
        QString nme;

        /** The names of each of the dimensions of the variable */
        QStringList dim_names;

        /** The size of each of the dimensions */
        QList<int> dim_sizes;

        /** The names of all attributes associated with the variable */
        QStringList att_names;
        /** The types of all of the attributes */
        QList<int> att_types;
        /** The values of all of the attributes */
        QList<QVariant> att_values;

        /** The hyperslab info */
        NetCDFHyperSlab slab;

        /** The ID number of the variable in the data file */
        int idnum;

        /** The type of the data in the data file */
        int xtyp;
    };

    /** This class holds the actual data read from a NetCDF file

        @author Christopher Woods
    */
    class SIREIO_EXPORT NetCDFData : public NetCDFDataInfo
    {

        friend class NetCDFFile;

    public:
        NetCDFData();

#ifndef SIRE_SKIP_INLINE_FUNCTIONS
        /** Construct a piece of NetCDF data called 'name', with the passed 'values',
            using the specified dimensions and dimension sizes, and optionally
            with the associated attributes */
        template <class T>
        NetCDFData(const QString &name, const QVector<T> &values, const QStringList &dimensions,
                   const QList<int> &dimension_sizes,
                   const QHash<QString, QVariant> &attributes = QHash<QString, QVariant>())
            : NetCDFDataInfo(0, name, get_nc_type<T>(), dimensions, dimension_sizes, get_attribute_names(attributes),
                             get_attribute_types(attributes), get_attribute_values(attributes))
        {
            // ensure that there is sufficient data in 'values' for the dimensions
            if (values.count() != this->nValues())
            {
                this->assertNValuesEquals(values.count());
            }

            memdata = QByteArray();
            memdata.resize(values.count() * sizeof(T));

            char *data = memdata.data();
            const char *orig = reinterpret_cast<const char *>(values.data());

            for (int i = 0; i < memdata.count(); ++i)
            {
                data[i] = orig[i];
            }
        }

        template <class T>
        NetCDFData(const NetCDFDataInfo &info,
                   const QVector<T> &values)
            : NetCDFDataInfo(info)
        {
            if (values.count() != this->nValues())
            {
                this->assertNValuesEquals(values.count());
            }

            memdata = QByteArray();
            memdata.resize(values.count() * sizeof(T));

            char *data = memdata.data();
            const char *orig = reinterpret_cast<const char *>(values.data());

            for (int i = 0; i < memdata.count(); ++i)
            {
                data[i] = orig[i];
            }
        }

        template <class T>
        NetCDFData(const NetCDFDataInfo &info,
                   const NetCDFHyperSlab &hyperslab,
                   const QVector<T> &values)
            : NetCDFDataInfo(info, hyperslab)
        {
            if (values.count() != this->nValues())
            {
                this->assertNValuesEquals(values.count());
            }

            memdata = QByteArray();
            memdata.resize(values.count() * sizeof(T));

            char *data = memdata.data();
            const char *orig = reinterpret_cast<const char *>(values.data());

            for (int i = 0; i < memdata.count(); ++i)
            {
                data[i] = orig[i];
            }
        }

#endif

        NetCDFData(const NetCDFData &other);

        ~NetCDFData();

        QVector<QVariant> toArray() const;

        QVector<float> toFloatArray() const;
        QVector<double> toDoubleArray() const;

        QVector<qint32> toInt32Array() const;
        QVector<qint64> toInt64Array() const;

        const void *data() const;

    protected:
        NetCDFData(const NetCDFDataInfo &info);
        NetCDFData(const NetCDFDataInfo &info,
                   const NetCDFHyperSlab &slab);

        void setData(const QByteArray &data);

    private:
        /** Raw memory containing the data */
        QByteArray memdata;
    };

    /** This class provides an internal interface to NetCDF files

        @author Christopher Woods
    */
    class SIREIO_EXPORT NetCDFFile : public boost::noncopyable
    {
    public:
        NetCDFFile();
        NetCDFFile(const QString &filename);

        ~NetCDFFile();

        bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly,
                  bool use_64bit_offset = true, bool use_netcdf4 = true);

        static QMutex *globalMutex();

        QString filename() const;

        QString getStringAttribute(const QString &name) const;

        QHash<QString, int> getDimensions() const;

        QHash<QString, NetCDFDataInfo> getVariablesInfo() const;

        NetCDFData read(const NetCDFDataInfo &variable) const;
        NetCDFData read(const NetCDFDataInfo &variable,
                        const NetCDFHyperSlab &slab) const;

        void close();

    protected:
        bool _lkr_open(QIODevice::OpenMode mode = QIODevice::ReadOnly,
                       bool use_64bit_offset = true, bool use_netcdf4 = true);

        void _lkr_writeHeader(const QHash<QString, QString> &globals,
                              const QHash<QString, NetCDFDataInfo> &dimensions);

        void _lkr_writeData(const NetCDFData &data);

        QString _lkr_getStringAttribute(const QString &name) const;

        QHash<QString, int> _lkr_getDimensions() const;

        QHash<QString, NetCDFDataInfo> _lkr_getVariablesInfo() const;

        NetCDFData _lkr_read(const NetCDFDataInfo &variable) const;
        NetCDFData _lkr_read(const NetCDFDataInfo &variable,
                             const NetCDFHyperSlab &slab) const;

        void _lkr_close();

    private:
        NetCDFFile(const QString &filename, bool overwrite_file, bool use_64bit_offset = true, bool use_netcdf4 = false);

        void _lkr_writeData(const QHash<QString, QString> &globals, const QHash<QString, NetCDFData> &data);

        int call_netcdf_function(std::function<int()> func, int ignored_error = 0) const;

        /** The name of the file */
        QString fname;

        /** Handle to the NetCDF file */
        int hndl;

        /** Global mutex used to ensure that all NetCDF operations
         *  are completely serialised */
        static QMutex mutex;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

    /** Return the ID number of this piece of data */
    SIRE_ALWAYS_INLINE int NetCDFDataInfo::ID() const
    {
        return idnum;
    }

    /** Return the name of this piece of data */
    SIRE_ALWAYS_INLINE QString NetCDFDataInfo::name() const
    {
        return nme;
    }

    /** Return the names of the dimensions of this data */
    SIRE_ALWAYS_INLINE QStringList NetCDFDataInfo::dimensions() const
    {
        return dim_names;
    }

    /** Return the number of values for each of the dimensions of this data */
    SIRE_ALWAYS_INLINE QList<int> NetCDFDataInfo::dimensionSizes() const
    {
        return dim_sizes;
    }

    /** Return the number of attributes of this data in the file */
    SIRE_ALWAYS_INLINE int NetCDFDataInfo::nAttributes() const
    {
        return att_names.count();
    }

#endif

} // namespace SireIO

SIRE_END_HEADER

#endif
