/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#include <Python.h>
#include <boost/python.hpp>

#include "sireglobal.h"

#include "SireError/getbacktrace.h"
#include "SireError/printerror.h"
#include "SireError/exception.h"

#include "wrap_SireError.h"

namespace SireError
{

void export_exceptions();

QHash<QString, QString> get_last_error_details();

void export_SireError()
{
    export_exceptions();

    boost::python::def( "getBackTrace", &SireError::getBackTrace );

    {
        typedef void (*printError_type)(const QString&);
        printError_type printError_value( &SireError::printError );

        boost::python::def( "printError", printError_value );
    }

    boost::python::def( "get_last_error_details", &get_last_error_details );
}

}
