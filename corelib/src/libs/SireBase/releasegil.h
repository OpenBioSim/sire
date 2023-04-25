/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2023  Christopher Woods
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
  *  at http://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREBASE_RELEASEGIL_H
#define SIREBASE_RELEASEGIL_H

#include "sireglobal.h"
#include <memory>

SIRE_BEGIN_HEADER

namespace SireBase
{
    namespace detail
    {
        class ReleaseGILBase;
    }

    typedef std::shared_ptr<detail::ReleaseGILBase> GILHandle;

    SIREBASE_EXPORT GILHandle release_gil();

    SIREBASE_EXPORT void ipython_clear_output(bool wait = true);

    SIREBASE_EXPORT bool sys_stdout_is_ipython();
    SIREBASE_EXPORT bool sys_stdout_is_atty();
    SIREBASE_EXPORT void sys_stdout_write(const QString &text, bool flush = false);
    SIREBASE_EXPORT void sys_stdout_move_up(int n);

    struct ANSI
    {
        enum Color
        {
            BLACK = 30,
            RED = 31,
            GREEN = 32,
            YELLOW = 33,
            BLUE = 34,
            MAGENTA = 35,
            CYAN = 36,
            WHITE = 37,
            DEFAULT = 39,
            BRIGHT_BLACK = 90,
            BRIGHT_RED = 91,
            BRIGHT_GREEN = 92,
            BRIGHT_YELLOW = 93,
            BRIGHT_BLUE = 94,
            BRIGHT_MAGENTA = 95,
            BRIGHT_CYAN = 96,
            BRIGHT_WHITE = 97
        };
    };

    SIREBASE_EXPORT QString esc_color(ANSI::Color fg, ANSI::Color bg = ANSI::DEFAULT, bool bold = false, bool underline = false);
    SIREBASE_EXPORT QString esc_color(ANSI::Color fg, bool bold, bool underline = false);
    SIREBASE_EXPORT QString esc_reset();

    namespace detail
    {
        /** This is the base class that will be sub-classed in the
         *  Python wrapper to provide the real code to release
         *  (and re-aquire) the GIL. We have to do this as
         *  the C++ libraries DO NOT link directly to Python. Only
         *  the wrappers link to Python.
         */
        class SIREBASE_EXPORT ReleaseGILBase
        {
            friend SIREBASE_EXPORT GILHandle SireBase::release_gil();

            friend SIREBASE_EXPORT void SireBase::ipython_clear_output(bool wait);
            friend SIREBASE_EXPORT void SireBase::sys_stdout_write(const QString &text, bool flush);
            friend SIREBASE_EXPORT bool SireBase::sys_stdout_is_ipython();
            friend SIREBASE_EXPORT void SireBase::sys_stdout_move_up(int n);
            friend SIREBASE_EXPORT bool SireBase::sys_stdout_is_atty();

        public:
            ReleaseGILBase();
            virtual ~ReleaseGILBase();

        protected:
            virtual bool stdout_is_atty() const = 0;
            virtual void stdout_write(const QString &text, bool flush) const = 0;
            virtual bool is_ipython() const = 0;
            virtual void ipython_clear(bool wait) const = 0;
            virtual void move_up(int n) const = 0;

            static void registerReleaseGIL(ReleaseGILBase *handle);
            virtual GILHandle releaseGIL() const = 0;

        private:
            /** Pointer to the concrete class that is registered by
             *  the wrappers
             */
            static ReleaseGILBase *handle;
        };
    }
}

SIRE_END_HEADER

#endif
