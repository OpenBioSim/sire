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

#ifndef SIREBASE_ATEXIT_H
#define SIREBASE_ATEXIT_H

#include "sireglobal.h"

#include <functional>

SIRE_BEGIN_HEADER

namespace SireBase
{
    SIREBASE_EXPORT void clean_up();

    SIREBASE_EXPORT void register_clean_up_function(std::function<void()> func);

    /** This class makes it easier to register a function to be called
     *  when the program exits. This is useful for cleaning up resources
     *  that are allocated during the program's execution.
     *
     *  Simply define your function (should be void func() { ... }) and
     *  then create a static instance of this class with the function as the
     *  argument. The constructor will be called at library load
     *  (static initialisation) and the function will be registered to be
     *  called at exit.
     *
     *  e.g.
     *
     *  void my_exit_function()
     *  {
     *     // clean up code here
     *  }
     *
     *  static RegisterExitFunction my_exit_function_instance(my_exit_function);
     *
     */
    class SIREBASE_EXPORT RegisterExitFunction
    {
    public:
        RegisterExitFunction(std::function<void()> func)
        {
            SireBase::register_clean_up_function(func);
        }

        ~RegisterExitFunction() {}
    };

} // namespace SireBase

SIRE_EXPOSE_FUNCTION(SireBase::clean_up);

SIRE_END_HEADER

#endif