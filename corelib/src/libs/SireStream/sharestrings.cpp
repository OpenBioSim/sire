/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "sharestrings.h"

#include <QMutex>
#include <QSet>

Q_GLOBAL_STATIC(QMutex, stringRegistryMutex)

static QSet<QString> *string_registry(0);

namespace SireStream
{

    /** This function adds the string 'string' to shared storage, and returns
        a copy of the shared-stored value. Use this function to ensure that
        there is only one copy of duplicate strings */
    QString shareString(const QString &string)
    {
        QMutex *mutex = stringRegistryMutex();

        if (not mutex)
            return string;

        QMutexLocker lkr(mutex);

        if (not string_registry)
        {
            string_registry = new QSet<QString>();
        }

        QSet<QString>::const_iterator it = string_registry->constFind(string);

        if (it != string_registry->constEnd())
            return *it;
        else
        {
            string_registry->insert(string);
            return string;
        }
    }

    /** This function adds all of the strings in 'strings' to shared storage */
    void shareStrings(QStringList &strings)
    {
        for (QStringList::iterator it = strings.begin(); it != strings.end(); ++it)
        {
            *it = shareString(*it);
        }
    }

} // namespace SireStream
