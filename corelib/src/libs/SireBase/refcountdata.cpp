/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
  *
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published
  *  by the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  *
  *  For full details of the license please see the LICENSE file
  *  that should have come with this distribution.
  *
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include "refcountdata.h"

#include "SireError/getbacktrace.h"

#include <QDebug>

using namespace SireBase;

Q_GLOBAL_STATIC(tbb::spin_mutex, get_spin_mutex);

namespace SireBase
{
    namespace detail
    {
        tbb::spin_mutex *get_shared_null_mutex()
        {
            return get_spin_mutex();
        }
    } // namespace detail
} // namespace SireBase

/** Constructor */
RefCountData::RefCountData()
{
}

/** Copying is forbidden */
RefCountData::RefCountData(const RefCountData &o)
{
}

/** Destructor */
RefCountData::~RefCountData()
{
}

/** Copy assignment is forbidden */
RefCountData &RefCountData::operator=(RefCountData &)
{
    return *this;
}

/** Only the same if exactly the same object */
bool RefCountData::operator==(const RefCountData &other) const
{
    return this == &other;
}

/** Only the same if exactly the same object */
bool RefCountData::operator!=(const RefCountData &other) const
{
    return not operator==(other);
}

/** Constructor */
RefCountData::Counter::Counter()
{
    refcount = 0;
    // qDebug() << "NEW REFCOUNTDATA" << qintptr(this);
}

/** Destructor */
RefCountData::Counter::~Counter()
{
    // qDebug() << "DELETING REFCOUNTDATA" << qintptr(this);
}
