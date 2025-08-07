/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2025  Christopher Woods
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

#include "cmapparameter.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QReadWriteLock>

using namespace SireMM;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<CMAPParameter> r_cmapparam(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const CMAPParameter &cmapparam)
{
    writeHeader(ds, r_cmapparam, 1);

    SharedDataStream sds(ds);

    sds << cmapparam.param;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, CMAPParameter &cmapparam)
{
    VersionID v = readHeader(ds, r_cmapparam);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> cmapparam.param;
    }
    else
        throw version_error(v, "1", r_cmapparam, CODELOC);

    return ds;
}

CMAPParameter::CMAPParameter()
{
}

Q_GLOBAL_STATIC(QVector<Array2D<double>>, shared_cache)

Q_GLOBAL_STATIC(QReadWriteLock, shared_cache_lock)

Array2D<double> deduplicate_parameter(const Array2D<double> &grid)
{
    QReadLocker locker(shared_cache_lock());

    int num_params = 0;

    for (const Array2D<double> &cached_grid : *shared_cache())
    {
        if (cached_grid == grid)
        {
            return cached_grid;
        }

        num_params += 1;
    }

    locker.unlock();

    QWriteLocker write_locker(shared_cache_lock());

    // make sure no-one else has added the parameter while we were waiting
    if (num_params == shared_cache()->size())
    {
        shared_cache()->append(grid);
        return grid;
    }
    else
    {
        write_locker.unlock();
        // try again
        return deduplicate_parameter(grid);
    }
}

namespace SireMM
{
    uint qHash(const CMAPParameter &param)
    {
        // because the parameters are de-duplicated, we can just hash the
        // pointer to the shared cache
        return ::qHash(quintptr(param.param.constData()));
    }
} // namespace SireMM

CMAPParameter::CMAPParameter(const Array2D<double> &grid)
{
    // put the grid into the shared cache, so that we only have
    // one copy of each different parameter - this will massively
    // reduce both memory consumption and the size of the stream
    // files
    this->param = deduplicate_parameter(grid);
}

CMAPParameter::CMAPParameter(const CMAPParameter &other)
{
    this->param = other.param;
}

CMAPParameter::~CMAPParameter()
{
}

CMAPParameter &CMAPParameter::operator=(const CMAPParameter &other)
{
    this->param = other.param;
    return *this;
}

bool CMAPParameter::operator==(const CMAPParameter &other) const
{
    // this works as the parameters have been de-duplicated
    return this->param.constData() == other.param.constData();
}

bool CMAPParameter::operator!=(const CMAPParameter &other) const
{
    return not CMAPParameter::operator==(other);
}

bool CMAPParameter::operator<(const CMAPParameter &other) const
{
    // the ordering is based on the smaller grid, and then comparing
    // the numbers one by one to find the first smaller value
    if (this->param.nRows() < other.param.nRows())
        return true;
    else if (this->param.nRows() == other.param.nRows() and this->param.nColumns() < other.param.nColumns())
        return true;
    else if (this->param.nRows() == other.param.nRows() and this->param.nColumns() == other.param.nColumns())
    {
        for (int i = 0; i < this->param.nRows(); i++)
        {
            for (int j = 0; j < this->param.nColumns(); j++)
            {
                if (this->param(i, j) < other.param(i, j))
                    return true;
                else if (this->param(i, j) > other.param(i, j))
                    return false;
            }
        }

        // must be equal if we get here
        return false;
    }
    else
        return false;
}

bool CMAPParameter::operator<=(const CMAPParameter &other) const
{
    return CMAPParameter::operator<(other) or CMAPParameter::operator==(other);
}

bool CMAPParameter::operator>(const CMAPParameter &other) const
{
    return not CMAPParameter::operator<=(other);
}

bool CMAPParameter::operator>=(const CMAPParameter &other) const
{
    return not CMAPParameter::operator<(other);
}

const char *CMAPParameter::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CMAPParameter>());
}

const char *CMAPParameter::what() const
{
    return CMAPParameter::typeName();
}

CMAPParameter *CMAPParameter::clone() const
{
    return new CMAPParameter(*this);
}

QString CMAPParameter::toString() const
{
    return QString("CMAPParameter(\n%1\n)").arg(param.toString());
}

const Array2D<double> &CMAPParameter::grid() const
{
    return param;
}

bool CMAPParameter::isEmpty() const
{
    return param.nRows() == 0 or param.nColumns() == 0;
}

bool CMAPParameter::isNull() const
{
    return this->isEmpty();
}

QVector<double> CMAPParameter::values() const
{
    return param.toRowMajorVector();
}

int CMAPParameter::nRows() const
{
    return param.nRows();
}

int CMAPParameter::nColumns() const
{
    return param.nColumns();
}
