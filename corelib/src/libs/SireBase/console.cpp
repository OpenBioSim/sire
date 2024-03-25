/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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

#include "SireBase/console.h"

using namespace SireBase;

ConsoleBase::ConsoleBase()
{
}

ConsoleBase::~ConsoleBase()
{
}

ConsoleBase *Console::c(0);

Console::Console()
{
}

Console::~Console()
{
}

/** Write the passed debug message to the console */
void Console::debug(const QString &message)
{
    if (c)
        c->debug(message);
}

/** Write the passed warning message to the console */
void Console::warning(const QString &message)
{
    if (c)
        c->warning(message);
}

/** Write the passed error message to the console */
void Console::error(const QString &message)
{
    if (c)
        c->error(message);
}

/** Write the passed info message to the console */
void Console::info(const QString &message)
{
    if (c)
        c->info(message);
}

/** Set the driver for the global console - this will delete
 *  any existing console - and will take ownership of the pointer!
 */
void Console::setConsole(ConsoleBase *console)
{
    if (c)
        delete c;

    c = console;
}
