/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIRECAS_EXPRESSIONS_H
#define SIRECAS_EXPRESSIONS_H

#include <QList>

#include "expression.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{

  class Expression;

  /**
  Expressions provides a list of expressions.

  @author Christopher Woods
  */
  class SIRECAS_EXPORT Expressions : public QList<Expression>
  {
  public:
    Expressions();

    Expressions(const Expression &expression);

    Expressions(const QList<Expression> &expressions);

    ~Expressions();

    Expressions differentiate(const Symbol &symbol) const;
    Expressions integrate(const Symbol &symbol) const;
  };

} // namespace SireCAS

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

namespace SireCAS
{
  /** Return a list of all children of type 'T' in this expression */
  template <class T>
  SIRE_OUTOFLINE_TEMPLATE QList<T> Expression::children() const
  {
    Expressions exs = this->children();

    QList<T> children_t;

    for (Expressions::const_iterator it = exs.constBegin(); it != exs.constEnd(); ++it)
    {
      const ExpressionBase &base = it->base();

// gccxml doesn't like this section, so remove it
// when we are generating the python wrappers
#ifndef SKIP_BROKEN_GCCXML_PARTS
      if (base.isA<T>())
        children_t.append(base.asA<T>());
#endif
    }

    return children_t;
  }
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

SIRE_END_HEADER

#endif
