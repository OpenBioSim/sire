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

#include "lazyevaluator.h"

#include "SireBase/property.h"

#include <boost/noncopyable.hpp>

// We have to undef the 'emit' from Qt as this is a function
// name used in TBB! This should be safe for Qt as that
// code uses Q_EMIT, and we don't use signals and slots in Sire
#ifdef emit
#undef emit
#endif

#include <tbb/collaborative_call_once.h>
#include <tbb/spin_rw_mutex.h>

#include <memory>

#include <QDebug>

namespace SireBase
{
    namespace detail
    {
        class LazyEvaluatorData : boost::noncopyable
        {
        public:
            LazyEvaluatorData() : boost::noncopyable()
            {
            }

            ~LazyEvaluatorData()
            {
            }

            PropertyPtr evaluate(const QString &key, std::function<PropertyPtr()> func);

        private:
            typedef std::pair<tbb::collaborative_once_flag, PropertyPtr> ResultBuffer;

            QHash<QString, std::shared_ptr<ResultBuffer>> results;
            tbb::spin_rw_mutex results_mutex;
        };
    }
}

using namespace SireBase;
using namespace SireBase::detail;

PropertyPtr LazyEvaluatorData::evaluate(const QString &key,
                                        std::function<PropertyPtr()> func)
{
    std::shared_ptr<ResultBuffer> buffer;

    // get a pointer to the buffer - this is not thread safe,
    // so we need to quickly look for it, and create it
    // if it doesn't exist
    {
        tbb::spin_rw_mutex::scoped_lock lock(results_mutex);

        auto it = results.find(key);

        if (it == results.end())
        {
            // need a write lock
            if (not lock.upgrade_to_writer())
            {
                // the lock has been upgraded to a writer, but it
                // had to be released. We need to double-check that
                // another thread hasn't already added the buffer
                it = results.find(key);
            }

            if (it == results.end())
            {
                buffer.reset(new ResultBuffer());
                results.insert(key, buffer);
            }
            else
            {
                buffer = it.value();
            }
        }
        else
        {
            buffer = it.value();
        }
    }

    // now we have the buffer, let's call the function once only!
    tbb::collaborative_call_once(buffer->first, [&]()
                                 { buffer->second = func(); });

    return buffer->second;
}

LazyEvaluator::LazyEvaluator() : d(new LazyEvaluatorData())
{
}

LazyEvaluator::LazyEvaluator(const LazyEvaluator &other)
    : d(other.d)
{
}

LazyEvaluator::~LazyEvaluator()
{
}

LazyEvaluator &LazyEvaluator::operator=(const LazyEvaluator &other)
{
    d = other.d;
    return *this;
}

PropertyPtr LazyEvaluator::evaluate(const QString &key,
                                    std::function<PropertyPtr()> func) const
{
    return const_cast<LazyEvaluatorData *>(d.get())->evaluate(key, func);
}
