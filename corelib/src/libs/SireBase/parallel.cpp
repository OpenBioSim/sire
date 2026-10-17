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

#include "parallel.h"

#include <QMutex>

#include <tbb/global_control.h>

namespace SireBase
{
    static tbb::global_control *global_tbb_control = 0;

    static QMutex global_mutex;

    int get_max_num_threads()
    {
        return tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    }

    void set_default_num_threads()
    {
        QMutexLocker lkr(&global_mutex);

        if (global_tbb_control != 0)
        {
            delete global_tbb_control;
            global_tbb_control = 0;
        }
    }

    bool should_run_in_parallel(int count, const PropertyMap &map)
    {
        if (count < 8)
            return false;
        else if (get_max_num_threads() <= 1)
            return false;
        else if (map["parallel"].hasValue())
            return map["parallel"].value().asA<BooleanProperty>().value();
        else
            return true;
    }

    void set_max_num_threads(int n)
    {
        if (n <= 0)
            return;

        QMutexLocker lkr(&global_mutex);

        if (global_tbb_control != 0)
        {
            delete global_tbb_control;
            global_tbb_control = 0;
        }

        global_tbb_control = new tbb::global_control(tbb::global_control::max_allowed_parallelism, n);
    }
}
