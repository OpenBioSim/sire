/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#ifndef SIRESEARCH_PARSER_H
#define SIRESEARCH_PARSER_H

#include "SireMol/parser.h"
#include "SireMol/select.h"

SIRE_BEGIN_HEADER

namespace SireSearch
{

    namespace parser
    {
        void set_token(const QString &token, const QString &selection);
        void reset_tokens();

        SireMol::parser::SelectEnginePtr parse(const QString &str);

        class SIRESEARCH_EXPORT SearchParser : public SireMol::parser::Parser
        {
        public:
            SearchParser();
            ~SearchParser();

            static void install();

            void set_token(const QString &token, const QString &selection);
            bool has_token(const QString &token);
            QString get_token(const QString &token);
            void delete_token(const QString &token);
            void delete_all_tokens();

            void reset_tokens();

            SireMol::parser::SelectEnginePtr parse(const QString &str);
        };
    } // namespace parser

} // namespace SireSearch

SIRE_END_HEADER

#endif
