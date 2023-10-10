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
  *  at https://sire.openbiosim.org
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
#include "SireUnits/temperature.h"

/** This is the grammar that enables skipping of spaces, newlines and comments */
template <typename IteratorT>
class SkipperGrammar : public qi::grammar<IteratorT>
{
public:
    SkipperGrammar() : SkipperGrammar::base_type(rule)
    {
        lineCommentRule = qi::lit("//") >> *(qi::char_);
        blockCommentRule = qi::lit("/*") >> *(qi::char_ - qi::lit("*/")) >> qi::lit("*/");
        rule = lineCommentRule | blockCommentRule;
    }

    qi::rule<IteratorT> lineCommentRule;
    qi::rule<IteratorT> blockCommentRule;
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
        short_unit_token.add(
            "cal", AST::Unit(SireUnits::cal))(
            "J", AST::Unit(SireUnits::joule))(
            "Ha", AST::Unit(SireUnits::hartree))(
            "mol", AST::Unit(SireUnits::mole))(
            "rad", AST::Unit(SireUnits::radian))(
            "°", AST::Unit(SireUnits::degree))(
            "Å", AST::Unit(SireUnits::angstrom))(
            "A", AST::Unit(SireUnits::angstrom))( // this replaces A for amp, but is more normal in MD
            "m", AST::Unit(SireUnits::meter))(
            "\"", AST::Unit(SireUnits::inch))(
            "in", AST::Unit(SireUnits::inch))(
            "'", AST::Unit(SireUnits::foot))(
            "ft", AST::Unit(SireUnits::foot))(
            "mph", AST::Unit(SireUnits::miles_per_hour))(
            "kph", AST::Unit(SireUnits::kilometers_per_hour))(
            "s", AST::Unit(SireUnits::second))(
            "g", AST::Unit(SireUnits::gram))(
            "N", AST::Unit(SireUnits::newton))(
            "Pa", AST::Unit(SireUnits::pascal))(
            "K", AST::Unit(SireUnits::kelvin))(
            "°K", AST::Unit(SireUnits::kelvin))(
            "V", AST::Unit(SireUnits::volt))(
            "F", AST::Unit(SireUnits::farad))(
            "W", AST::Unit(SireUnits::watt))(
            "e", AST::Unit(SireUnits::e_charge))(
            "|e|", AST::Unit(SireUnits::mod_electron))(
            "C", AST::Unit(SireUnits::coulomb)); // this is Coulomb, but lots of people use it for celsius

        unit_token.add(
            "calorie", AST::Unit(SireUnits::cal))(
            "joule", AST::Unit(SireUnits::joule))(
            "hartree", AST::Unit(SireUnits::hartree))(
            "mole", AST::Unit(SireUnits::mole))(
            "dozen", AST::Unit(SireUnits::dozen))(
            "radian", AST::Unit(SireUnits::radian))(
            "degree", AST::Unit(SireUnits::degree))(
            "angstrom", AST::Unit(SireUnits::angstrom))(
            "meter", AST::Unit(SireUnits::meter))(
            "bohr", AST::Unit(SireUnits::bohr_radii))(
            "inch", AST::Unit(SireUnits::inch))(
            "inches", AST::Unit(SireUnits::inch))(
            "foot", AST::Unit(SireUnits::foot))(
            "feet", AST::Unit(SireUnits::foot))(
            "yard", AST::Unit(SireUnits::yard))(
            "mile", AST::Unit(SireUnits::mile))(
            "second", AST::Unit(SireUnits::second))(
            "minute", AST::Unit(SireUnits::minute))(
            "hour", AST::Unit(SireUnits::hour))(
            "day", AST::Unit(SireUnits::day))(
            "week", AST::Unit(SireUnits::week))(
            "fortnight", AST::Unit(SireUnits::fortnight))(
            "akma", AST::Unit(SireUnits::akma_time))(
            "dalton", AST::Unit(SireUnits::g_per_mol * SireUnits::mole))(
            "gram", AST::Unit(SireUnits::gram))(
            "tonne", AST::Unit(SireUnits::tonne))(
            "newton", AST::Unit(SireUnits::newton))(
            "ounce", AST::Unit(SireUnits::ounce))(
            "pound", AST::Unit(SireUnits::pound))(
            "stone", AST::Unit(SireUnits::stone))(
            "hundredweight", AST::Unit(SireUnits::hundredweight))(
            "pascal", AST::Unit(SireUnits::pascal))(
            "bar", AST::Unit(SireUnits::bar))(
            "atm", AST::Unit(SireUnits::atm))(
            "atmosphere", AST::Unit(SireUnits::atm))(
            "psi", AST::Unit(SireUnits::psi))(
            "mmHg", AST::Unit(SireUnits::mmHg))(
            "kelvin", AST::Unit(SireUnits::kelvin))(
            "amp", AST::Unit(SireUnits::amp))(
            "ampere", AST::Unit(SireUnits::amp))(
            "volt", AST::Unit(SireUnits::volt))(
            "farad", AST::Unit(SireUnits::farad))(
            "watt", AST::Unit(SireUnits::watt))(
            "electron", AST::Unit(SireUnits::e_charge))(
            "e_charge", AST::Unit(SireUnits::e_charge))(
            "mod_electron", AST::Unit(SireUnits::mod_electron))(
            "faraday", AST::Unit(SireUnits::faraday))(
            "coulomb", AST::Unit(SireUnits::coulomb))(
            "kcal_per_mol", AST::Unit(SireUnits::kcal_per_mol))( // this is commonly used
            "kJ_per_mol", AST::Unit(SireUnits::kJ_per_mol));     // this is commonly used

        short_prefix_token.add(
            "d", AST::Prefix(1e-1))(
            "c", AST::Prefix(1e-2))(
            "m", AST::Prefix(1e-3))(
            "u", AST::Prefix(1e-6))(
            "µ", AST::Prefix(1e-6))(
            "μ", AST::Prefix(1e-6))(
            "n", AST::Prefix(1e-9))(
            "p", AST::Prefix(1e-12))(
            "f", AST::Prefix(1e-15))(
            "a", AST::Prefix(1e-18))(
            "z", AST::Prefix(1e-21))(
            "y", AST::Prefix(1e-24))(
            "r", AST::Prefix(1e-27))(
            "q", AST::Prefix(1e-30))(
            "da", AST::Prefix(10.0))(
            "h", AST::Prefix(1e2))(
            "k", AST::Prefix(1e3))(
            "M", AST::Prefix(1e6))(
            "G", AST::Prefix(1e9))(
            "T", AST::Prefix(1e12))(
            "P", AST::Prefix(1e15))(
            "E", AST::Prefix(1e18))(
            "Z", AST::Prefix(1e21))(
            "Y", AST::Prefix(1e24))(
            "R", AST::Prefix(1e27))(
            "Q", AST::Prefix(1e30));

        prefix_token.add(
            "deci", AST::Prefix(1e-1))(
            "centi", AST::Prefix(1e-2))(
            "milli", AST::Prefix(1e-3))(
            "micro", AST::Prefix(1e-6))(
            "nano", AST::Prefix(1e-9))(
            "pico", AST::Prefix(1e-12))(
            "femto", AST::Prefix(1e-15))(
            "atto", AST::Prefix(1e-18))(
            "zepto", AST::Prefix(1e-21))(
            "yocto", AST::Prefix(1e-24))(
            "ronto", AST::Prefix(1e-27))(
            "quecto", AST::Prefix(1e-30))(
            "deca", AST::Prefix(10.0))(
            "hecto", AST::Prefix(1e2))(
            "kilo", AST::Prefix(1e3))(
            "mega", AST::Prefix(1e6))(
            "giga", AST::Prefix(1e9))(
            "tera", AST::Prefix(1e12))(
            "peta", AST::Prefix(1e15))(
            "exa", AST::Prefix(1e18))(
            "zetta", AST::Prefix(1e21))(
            "yotta", AST::Prefix(1e24))(
            "ronna", AST::Prefix(1e27))(
            "quetta", AST::Prefix(1e30));

        // convenient shortcuts for the unit separator
        static const auto sep = qi::blank;

        // convenient shortcuts for the brackets
        static const auto leftB = qi::lit("(");
        static const auto rightB = qi::lit(")");

        // root rule to read a node as a single expression and no further input (eoi)
        nodeRule = (expressionRule >> qi::eoi) |
                   (expressionRule >> sep >> qi::eoi) |
                   (sep >> expressionRule >> qi::eoi) |
                   (sep >> expressionRule >> sep >> qi::eoi);

        unitRule = eps[_val = AST::Unit()] >>
                       (prefix_token[_val = _1] >> unit_token[_val *= _1] >> qi::lit("s")) |
                   (prefix_token[_val = _1] >> unit_token[_val *= _1]) |
                   (unit_token[_val = _1] >> qi::lit("s")) |
                   (unit_token[_val = _1]) |
                   (short_prefix_token[_val = _1] >> short_unit_token[_val *= _1]) |
                   (short_unit_token[_val = _1]);

        powerRule = eps[_val = AST::Power()] >>
                        (*sep >> (qi::lit("**") | qi::lit("^")) >> *sep >> int_[_val *= _1]) |
                    int_[_val *= _1];

        fullUnitRule = eps[_val = AST::FullUnit()] >>
                           (leftB >> *sep >> unitRule[_val = _1] >> *sep >> rightB >> -powerRule[_val *= _1]) |
                       (leftB >> *sep >> unitRule[_val = _1] >> *sep >> rightB) |
                       (unitRule[_val = _1] >> sep >> qi::lit("per") >> sep >> unitRule[_val /= _1]) |
                       (unitRule[_val = _1] >> -powerRule[_val *= _1]) |
                       (unitRule[_val = _1]);

        expressionPartRule = eps[_val = AST::Expression()] >>
                                 (fullUnitRule[_val = _1] >> powerRule[_val *= _1]) |
                             (leftB >> *sep >> fullUnitRule[_val = _1] >> *sep >> rightB >> powerRule[_val *= _1]) |
                             (fullUnitRule[_val = _1]) |
                             (leftB >> *sep >> fullUnitRule[_val = _1] >> *sep >> rightB);

        // rule for the left hand side of an expression - this can only
        // be an ExpressionPartRule, or an ExpressionRule that is enclosed
        // in brackets
        lhsRule = eps[_val = AST::Expression()] >>
                      (expressionPartRule[_val = _1] >> powerRule[_val *= _1]) |
                  (leftB >> *sep >> expressionRule[_val = _1] >> *sep >> rightB >> powerRule[_val *= _1]) |
                  (expressionPartRule[_val = _1]) |
                  (leftB >> *sep >> expressionRule[_val = _1] >> *sep >> rightB);

        // rule for the right hand side of an expression - this can only
        // be an ExpressionRule
        rhsRule = eps[_val = AST::Expression()] >>
                      (expressionRule[_val = _1] >> powerRule[_val *= _1]) |
                  (leftB >> *sep >> expressionRule[_val = _1] >> *sep >> rightB >> powerRule[_val *= _1]) |
                  (expressionRule[_val = _1]) |
                  (leftB >> *sep >> expressionRule[_val = _1] >> *sep >> rightB);

        binaryRule = eps[_val = AST::Expression()] >>
                         (lhsRule[_val = _1] >> qi::lit(".") >> rhsRule[_val *= _1]) |
                     (lhsRule[_val = _1] >> *sep >> qi::lit("*") >> *sep >> rhsRule[_val *= _1]) |
                     (lhsRule[_val = _1] >> *sep >> qi::lit("/") >> *sep >> rhsRule[_val /= _1]) |
                     (lhsRule[_val = _1] >> sep >> rhsRule[_val *= _1]) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> qi::lit(".") >> rhsRule[_val *= _1] >> *sep >> rightB >> powerRule[_val *= _1]) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> *sep >> qi::lit("*") >> *sep >> rhsRule[_val *= _1] >> *sep >> rightB >> powerRule[_val *= _1]) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> *sep >> qi::lit("/") >> *sep >> rhsRule[_val /= _1] >> *sep >> rightB >> powerRule[_val *= _1]) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> sep >> rhsRule[_val *= _1] >> *sep >> rightB >> powerRule[_val *= _1]) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> qi::lit(".") >> rhsRule[_val *= _1] >> *sep >> rightB) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> *sep >> qi::lit("*") >> *sep >> rhsRule[_val *= _1] >> *sep >> rightB) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> *sep >> qi::lit("/") >> *sep >> rhsRule[_val /= _1] >> *sep >> rightB) |
                     (leftB >> *sep >> lhsRule[_val = _1] >> sep >> rhsRule[_val *= _1] >> *sep >> rightB);

        expressionRule = eps[_val = AST::Expression()] >>
                             (binaryRule[_val = _1]) |
                         (expressionPartRule[_val = _1]) |
                         (leftB >> *sep >> expressionRule[_val = _1] >> *sep >> rightB);

        nodeRule.name("Node");
        expressionRule.name("Expression");
        fullUnitRule.name("FullUnit");
        unitRule.name("Unit");
        powerRule.name("PowerRule");
        lhsRule.name("LHSRule");
        rhsRule.name("RHSRule");
        binaryRule.name("BinaryRule");

        // action on failure to parse the string using the grammar
        on_error<fail>(nodeRule, std::cout << val("Error! Expecting ") << _4 // what failed?
                                           << val(" here: \"")
                                           << construct<std::string>(_3, _2) // iterators to error-pos, end
                                           << val("\"") << std::endl);
    }

    qi::symbols<char, AST::Unit> unit_token;
    qi::symbols<char, AST::Prefix> prefix_token;
    qi::symbols<char, AST::Unit> short_unit_token;
    qi::symbols<char, AST::Prefix> short_prefix_token;

    qi::rule<IteratorT, AST::Node(), SkipperT> nodeRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> lhsRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> rhsRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> binaryRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> expressionRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> expressionPartRule;
    qi::rule<IteratorT, AST::Unit> unitRule;
    qi::rule<IteratorT, AST::FullUnit(), SkipperT> fullUnitRule;
    qi::rule<IteratorT, AST::Power(), SkipperT> powerRule;
};

#endif
