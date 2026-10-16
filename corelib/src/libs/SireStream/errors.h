/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIRESTREAM_ERRORS_H
#define SIRESTREAM_ERRORS_H

#include "SireError/exception.h"

SIRE_BEGIN_HEADER

namespace SireStream
{

    /** This exception is thrown when corrupted data is detected

        @author Christopher Woods
    */
    class SIRESTREAM_EXPORT corrupted_data : public SireError::exception
    {
    public:
        corrupted_data() : SireError::exception()
        {
        }

        corrupted_data(QString err, QString place = QString()) : SireError::exception(err, place)
        {
        }

        corrupted_data(const corrupted_data &other) : SireError::exception(other)
        {
        }

        ~corrupted_data() throw()
        {
        }

        static const char *typeName();

        const char *what() const throw()
        {
            return corrupted_data::typeName();
        }

        void throwSelf() const
        {
            throw corrupted_data(*this);
        }
    };

} // namespace SireStream

Q_DECLARE_METATYPE(SireStream::corrupted_data)

SIRE_END_HEADER

#endif
