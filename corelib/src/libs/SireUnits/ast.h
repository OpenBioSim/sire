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

#ifndef SIREUNITS_AST_H
#define SIREUNITS_AST_H

#include "SireUnits/generalunit.h"

#include "SireError/errors.h"

SIRE_BEGIN_HEADER

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>

#include <boost/phoenix.hpp>
#include <boost/phoenix/fusion.hpp>
#include <boost/phoenix/object.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/stl.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

#include <memory>

// A lot of the below code is heavily inspired by
// https://medium.com/@alinakipoglu/parsing-with-spirit-qi-fcaeaf4357b3

namespace SireUnits
{

    /** Namespace holding the objects used in the abstract syntax tree */
    namespace AST
    {
        struct Node;

        /** The root node of the AST - this holds a single Expression */
        struct Node
        {
            QString toString() const;
        };

    } // namespace AST

} // namespace SireUnits

SIRE_END_HEADER

#endif
