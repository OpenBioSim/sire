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

#ifndef SIREBASE_CONSOLE_H
#define SIREBASE_CONSOLE_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
    /** Virtual base class of the actual console implementation */
    class SIREBASE_EXPORT ConsoleBase
    {
    public:
        ConsoleBase();
        virtual ~ConsoleBase();

        virtual void debug(const QString &message) const = 0;
        virtual void warning(const QString &message) const = 0;
        virtual void error(const QString &message) const = 0;
        virtual void info(const QString &message) const = 0;
    };

    /** This class provides static functions that can be used
     *  to (sparingly) output messages to the console. This is
     *  controlled by the tools in sire.utils.console
     */
    class SIREBASE_EXPORT Console
    {
    public:
        static void debug(const QString &message);
        static void warning(const QString &message);
        static void error(const QString &message);
        static void info(const QString &message);

        static void setConsole(ConsoleBase *console);

    private:
        Console();
        ~Console();

        static ConsoleBase *c;
    };
}

SIRE_END_HEADER

#endif // SIREBASE_CONSOLE_H
