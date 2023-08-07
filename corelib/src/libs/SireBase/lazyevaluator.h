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

#ifndef SIREBASE_LAZYEVALUATOR_H
#define SIREBASE_LAZYEVALUATOR_H

#include "SireBase/property.h"

#include <functional>
#include <memory>

SIRE_BEGIN_HEADER

namespace SireBase
{
    namespace detail
    {
        class LazyEvaluatorData;
    }

    /** This class provides a LazyEvaluator that can be used to
     *  lazily evalute and cache the result of calling a function
     *  once only. This uses tbb::collaborative_call_once to ensure
     *  that other threads that block on this are able to go down
     *  and help the calling thread with any downstream parallel
     *  sections
     */
    class SIREBASE_EXPORT LazyEvaluator
    {
    public:
        LazyEvaluator();
        LazyEvaluator(const LazyEvaluator &other);

        ~LazyEvaluator();

        LazyEvaluator &operator=(const LazyEvaluator &other);

        PropertyPtr evaluate(const QString &key, std::function<PropertyPtr()> func) const;

    private:
        std::shared_ptr<detail::LazyEvaluatorData> d;
    };

} // namespace SireBase

SIRE_EXPOSE_CLASS(SireBase::LazyEvaluator)

SIRE_END_HEADER

#endif
