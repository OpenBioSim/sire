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

#ifndef SIREMOL_PARSER_H
#define SIREMOL_PARSER_H

#include "select.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

    namespace parser
    {
        class SIREMOL_EXPORT Parser
        {
        public:
            Parser()
            {
            }

            virtual ~Parser()
            {
            }

            static void install_parser(Parser *parser);

            virtual bool has_token(const QString &token) = 0;

            virtual QString get_token(const QString &token) = 0;

            virtual void set_token(const QString &token, const QString &selection) = 0;

            virtual void delete_token(const QString &token) = 0;

            virtual void delete_all_tokens() = 0;

            virtual void reset_tokens() = 0;

            virtual SireMol::parser::SelectEnginePtr parse(const QString &str) = 0;

            static Parser &globalParser();

        private:
            static Parser *global_parser;
        };

        void set_token(const QString &token, const QString &selection);
        void reset_tokens();

        SireMol::parser::SelectEnginePtr parse(const QString &str);
    } // namespace parser

} // namespace SireMol

SIRE_END_HEADER

#endif
