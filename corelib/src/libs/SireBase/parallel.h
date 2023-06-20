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

#ifndef SIREBASE_PARALLEL_H
#define SIREBASE_PARALLEL_H

#include "sireglobal.h"

#include "SireBase/booleanproperty.h"
#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

#ifndef GCCXML_PARSE

#include <QMutex>
#include <QVector>

// We have to undef the 'emit' from Qt as this is a function
// name used in TBB! This should be safe for Qt as that
// code uses Q_EMIT, and we don't use signals and slots in Sire
#ifdef emit
#undef emit
#endif

#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_reduce.h>

#include <memory>

#endif // skip TBB when parsing with GCCXML

namespace SireBase
{
    SIREBASE_EXPORT bool should_run_in_parallel(int count, const PropertyMap &map = PropertyMap());

    SIREBASE_EXPORT int get_max_num_threads();

    SIREBASE_EXPORT void set_max_num_threads(int n);

    SIREBASE_EXPORT void set_default_num_threads();

#ifndef GCCXML_PARSE
    /** This function runs the passed array T of functions in parallel, if
        the optional 'run_parallel' is true. Otherwise, it runs the functions
        serially, one after another */
    template <class T>
    void parallel_invoke(const T &functions, bool run_parallel = true)
    {
        if (run_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, functions.count()), [&](const tbb::blocked_range<int> &r)
                              {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                functions[i]();
            } });
        }
        else
        {
            for (int i = 0; i < functions.count(); ++i)
            {
                functions[i]();
            }
        }
    }
#endif

} // end of namespace SireBase

SIRE_EXPOSE_FUNCTION(SireBase::should_run_in_parallel)
SIRE_EXPOSE_FUNCTION(SireBase::get_max_num_threads)
SIRE_EXPOSE_FUNCTION(SireBase::set_max_num_threads)
SIRE_EXPOSE_FUNCTION(SireBase::set_default_num_threads)

SIRE_END_HEADER

#endif
