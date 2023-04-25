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

#include "releasegil.h"

#include <QTextStream>
#include <QMutex>

#include <QDebug>

namespace SireBase
{
    /** Release the GIL. This returns a handle which, when destroyed,
     *  will automatically re-aquire the GIL (RAII pattern)
     */
    GILHandle release_gil()
    {
        if (detail::ReleaseGILBase::handle != 0)
        {
            return detail::ReleaseGILBase::handle->releaseGIL();
        }
        else
        {
            return GILHandle();
        }
    }

    /** Print the passed string on the Python STDOUT. This is needed
     *  when STDOUT is a different stream or buffered, e.g. when running
     *  in a Jupyter notebook cell
     */
    void sys_stdout_write(const QString &text, bool flush)
    {
        if (detail::ReleaseGILBase::handle != 0)
        {
            detail::ReleaseGILBase::handle->stdout_write(text, flush);
        }
        else
        {
            static QTextStream ts(stdout);
            ts << text;

            if (flush)
                ts.flush();
        }
    }

    void ipython_clear_output(bool wait)
    {
        if (sys_stdout_is_ipython() and detail::ReleaseGILBase::handle != 0)
        {
            detail::ReleaseGILBase::handle->ipython_clear(wait);
        }
    }

    bool sys_stdout_is_ipython()
    {
        if (detail::ReleaseGILBase::handle != 0)
        {
            return detail::ReleaseGILBase::handle->is_ipython();
        }
        else
            return false;
    }

    void sys_stdout_move_up(int n)
    {
        if (n <= 0)
            return;
        else if (n > 3200)
            n = 3200;

        if (sys_stdout_is_ipython())
        {
            // we don't have cursor control, so all we can do is
            // completely clear the cell
            detail::ReleaseGILBase::handle->ipython_clear(true);
        }
        else
        {
            sys_stdout_write(QString("\x1b[%1A").arg(n));
        }
    }

    QString esc_color(ANSI::Color fg, ANSI::Color bg, bool bold, bool underline)
    {
        QString code = QString::number(int(fg));

        if (bg != ANSI::DEFAULT)
            code += QString(";%1").arg(int(bg) + 10);

        if (bold)
            code += ";1";

        if (underline)
            code += ";4";

        return QString("\x1b[%1m").arg(code);
    }

    QString esc_color(ANSI::Color fg, bool bold, bool underline)
    {
        return esc_color(fg, ANSI::DEFAULT, bold, underline);
    }

    QString esc_reset()
    {
        return QString("\x1b[0m");
    }

    namespace detail
    {
        ReleaseGILBase::ReleaseGILBase()
        {
        }

        ReleaseGILBase::~ReleaseGILBase()
        {
        }

        ReleaseGILBase *ReleaseGILBase::handle(0);

        void ReleaseGILBase::registerReleaseGIL(ReleaseGILBase *h)
        {
            if (handle != 0)
            {
                delete handle;
            }

            handle = h;
        }
    }
}
