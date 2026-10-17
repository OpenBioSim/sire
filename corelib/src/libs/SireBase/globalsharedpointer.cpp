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

#include <QHash>

#include "globalsharedpointer.hpp"

using namespace SireBase;

Q_GLOBAL_STATIC(QMutex, globalPointerMutex);

/** The registry of all global shared pointers */
static QHash<QString, QSet<const void *>> *global_pointers(0);

GlobalSharedPointerBase::GlobalSharedPointerBase()
{
}

GlobalSharedPointerBase::~GlobalSharedPointerBase()
{
}

QMutex *GlobalSharedPointerBase::registryMutex()
{
    return globalPointerMutex();
}

QSet<const void *> &GlobalSharedPointerBase::getRegistry(const char *typname)
{
    if (not global_pointers)
        global_pointers = new QHash<QString, QSet<const void *>>();

    return (*global_pointers)[QLatin1String(typname)];
}
