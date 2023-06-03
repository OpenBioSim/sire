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

#ifndef SIREUNITS_GRAMMAR_H
#define SIREUNITS_GRAMMAR_H

using qi::_1;
using qi::double_;
using qi::eps;
using qi::fail;
using qi::int_;
using qi::lexeme;
using qi::lit;
using qi::on_error;
using namespace qi::labels;
using qi::as_string;

using phoenix::construct;
using phoenix::val;

using boost::spirit::ascii::char_;

#include "SireUnits/units.h"

/** This is the grammar that enables skipping of spaces, newlines and comments */
template <typename IteratorT>
class SkipperGrammar : public qi::grammar<IteratorT>
{
public:
    SkipperGrammar() : SkipperGrammar::base_type(rule)
    {
        lineCommentRule = qi::lit("//") >> *(qi::char_ - qi::eol) >> qi::eol;
        blockCommentRule = qi::lit("/*") >> *(qi::char_ - qi::lit("*/")) >> qi::lit("*/");
        spaceRule = qi::space;
        rule = spaceRule | lineCommentRule | blockCommentRule;
    }

    qi::rule<IteratorT> lineCommentRule;
    qi::rule<IteratorT> blockCommentRule;
    qi::rule<IteratorT> spaceRule;
    qi::rule<IteratorT> rule;
};

/** This is a quoted string grammar that will parse quoted strings and also
    auto-escape characters */
template <typename IteratorT, typename SkipperT>
class ValueGrammar : public qi::grammar<IteratorT, std::string(), SkipperT>
{
public:
    ValueGrammar() : ValueGrammar::base_type(rule, "String")
    {
        escapedStringRule %=
            qi::lexeme[qi::lit("'") >> *(escapeCharSymbols | (qi::char_ - qi::char_("'"))) >> qi::lit("'")];

        rawStringRule %= qi::lexeme[+(qi::alnum | qi::char_('.') | qi::char_('/') | qi::char_('_') | qi::char_('-'))];

        rule %= rawStringRule | escapedStringRule;

        escapeCharSymbols.add("\\a", '\a')("\\b", '\b')("\\f", '\f')("\\n", '\n')("\\r", '\r')("\\t", '\t')(
            "\\v", '\v')("\\\\", '\\')("\\\'", '\'')("\\\"", '\"');

        escapedStringRule.name("Escaped String");
        rawStringRule.name("Escaped String");

        escapeCharSymbols.name("Escaped Chars");
    }

    qi::rule<IteratorT, std::string(), SkipperT> escapedStringRule;
    qi::rule<IteratorT, std::string(), SkipperT> rawStringRule;
    qi::rule<IteratorT, std::string(), SkipperT> rule;
    qi::symbols<const char, const char> escapeCharSymbols;
};

/** This the main grammar for the selection statements */
template <typename IteratorT, typename SkipperT>
class Grammar : public qi::grammar<IteratorT, AST::Node(), SkipperT>
{
public:
    Grammar() : Grammar::base_type(nodeRule, "Node")
    {
        unit_token.add(
            "calorie", AST::Unit(SireUnits::cal))(
            "joule", AST::Unit(SireUnits::joule))(
            "mole", AST::Unit(SireUnits::mole))(
            "mol", AST::Unit(SireUnits::mole))(
            "dozen", AST::Unit(SireUnits::dozen))(
            "radian", AST::Unit(SireUnits::radian))(
            "degree", AST::Unit(SireUnits::degree))(
            "angstrom", AST::Unit(SireUnits::angstrom))(
            "meter", AST::Unit(SireUnits::meter));

        // root rule to read a node as a single expression and no further input (eoi)
        nodeRule = expressionRule >> qi::eoi;

        // convenient shortcuts for the brackets
        static const auto leftB = qi::lit("(");
        static const auto rightB = qi::lit(")");

        powerRule = eps[_val = AST::Power()] >>
                        double_[_val *= _1] |
                    (qi::lit("**") >> double_[_val *= _1]) |
                    (qi::lit("^") >> double_[_val *= _1]);

        fullUnitRule = eps[_val = AST::FullUnit()] >>
                           (unit_token[_val += _1]) |
                       (leftB >> unit_token[_val += _1] >> rightB) |
                       (unit_token[_val += _1] >> -powerRule[_val *= _1]) |
                       (leftB >> unit_token[_val += _1] >> rightB >> -powerRule[_val *= _1]);

        expressionRule = fullUnitRule;

        nodeRule.name("Node");
        expressionRule.name("Expression");
        fullUnitRule.name("FullUnit");
        powerRule.name("PowerRule");

        // action on failure to parse the string using the grammar
        on_error<fail>(nodeRule, std::cout << val("Error! Expecting ") << _4 // what failed?
                                           << val(" here: \"")
                                           << construct<std::string>(_3, _2) // iterators to error-pos, end
                                           << val("\"") << std::endl);
    }

    qi::symbols<char, AST::Unit> unit_token;

    qi::rule<IteratorT, AST::Node(), SkipperT> nodeRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> expressionRule;
    qi::rule<IteratorT, AST::FullUnit(), SkipperT> fullUnitRule;
    qi::rule<IteratorT, AST::Power(), SkipperT> powerRule;
};

#endif
