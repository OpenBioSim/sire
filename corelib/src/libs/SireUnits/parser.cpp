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

#define BOOST_SPIRIT_USE_PHOENIX_V3
#define BOOST_SPIRIT_UNICODE

#include "SireError/errors.h"

#include "SireUnits/generalunit.h"
#include "SireUnits/ast.h"

#include <QDebug>
#include <QHash>
#include <QMutex>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/fusion.hpp>
#include <boost/phoenix/object.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/phoenix/stl.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

using namespace SireUnits;
using namespace SireUnits::Dimension;

namespace spirit = boost::spirit;
namespace qi = spirit::qi;
namespace unicode = boost::spirit::unicode;
namespace phoenix = boost::phoenix;

#include "SireUnits/grammar.h" //file containing the actual grammar - separated to ease reading

template <typename IteratorT>
QString toString(IteratorT begin, IteratorT end)
{
    QStringList lines;
    for (; begin != end; ++begin)
    {
        lines.append(QString(*begin));
    }

    return lines.join("");
}

/** Function that parses the passed string (represented via its iterators) into an AST::Node */
template <typename IteratorT>
SireUnits::AST::Node parse(const IteratorT &begin, const IteratorT &end)
{
    using LinePosIteratorT = spirit::line_pos_iterator<IteratorT>;

    using SkipperGrammarT = SkipperGrammar<LinePosIteratorT>;
    using ParserGrammarT = Grammar<LinePosIteratorT, SkipperGrammarT>;

    SkipperGrammarT skipper;
    ParserGrammarT grammar;
    LinePosIteratorT posIterBegin(begin);
    LinePosIteratorT posIterEnd(end);

    SireUnits::AST::Node result;

    const bool parseResult = qi::phrase_parse(posIterBegin, posIterEnd, grammar, skipper, result);

    if (not(parseResult && posIterBegin == posIterEnd))
    {
        QString line = toString(LinePosIteratorT(begin), LinePosIteratorT(end));
        QString left = toString(posIterBegin, posIterEnd);

        if (line == left)
            throw SireError::io_error(QObject::tr("Failed to parse any of the unit '%1'.").arg(line), CODELOC);
        else
            throw SireError::io_error(QObject::tr("Failed to parse the unit '%1'. "
                                                  "Successfully parsed the beginning, but failed to parse '%2'.")
                                          .arg(line)
                                          .arg(left),
                                      CODELOC);
    }

    return result;
}

/** Function used internally to parse a string into an AST::Node */
static SireUnits::AST::Node parse_main(const std::string &str)
{
    // Read file contents.
    return parse(str.begin(), str.end());
}

/** Construct a unit (or value * unit) from the passed string */
GeneralUnit GeneralUnit::fromString(const QString &value)
{
    return GeneralUnit(value);
}

/** Construct a unit (or a value * unit) from the passed string */
GeneralUnit::GeneralUnit(const QString &value) : Unit(0)
{
}

/** Construct a value as 'value' * 'unit', where 'unit' is
 *  interpreted from the passed string
 */
GeneralUnit::GeneralUnit(double value, const QString &unit) : Unit(0)
{
    QString processed_unit = unit.simplified();
    processed_unit = processed_unit.toLower();

    if (processed_unit.length() == 0)
    {
        // this is a dimensionless value
        this->operator=(GeneralUnit(value));
        return;
    }

    // try to get the unit from the passed string
    auto ast = ::parse_main(processed_unit.toStdString());

    // this->operator=(value * parsed_unit);
}
