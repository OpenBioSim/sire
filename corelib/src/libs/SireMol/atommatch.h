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

#ifndef SIREMOL_ATOMMATCH_H
#define SIREMOL_ATOMMATCH_H

#include "SireMol/atom.h"
#include "SireMol/selector.hpp"
#include "SireMol/selectorm.hpp"
#include "SireMol/partialmolecule.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class AtomMatch;
    class AtomMatchM;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::AtomMatch &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::AtomMatch &);

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::AtomMatchM &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::AtomMatchM &);

namespace SireMol
{
    /** This class holds the results of performing a match on a molecule. */
    class SIREMOL_EXPORT AtomMatch : public SireBase::ConcreteProperty<AtomMatch, PartialMolecule>
    {

        friend QDataStream & ::operator<<(QDataStream &, const AtomMatch &);
        friend QDataStream & ::operator>>(QDataStream &, AtomMatch &);

    public:
        AtomMatch();
        AtomMatch(const MoleculeView &molview);
        AtomMatch(const Selector<Atom> &molview, const QList<qint64> &matches);
        AtomMatch(const Selector<Atom> &molview, const QList<QList<qint64>> &matches);

        AtomMatch(const AtomMatch &other);
        virtual ~AtomMatch();

        static const char *typeName();

        virtual const char *what() const
        {
            return AtomMatch::typeName();
        }

        virtual AtomMatch *clone() const
        {
            return new AtomMatch(*this);
        }

        AtomMatch &operator=(const AtomMatch &other);

        bool operator==(const AtomMatch &other) const;
        bool operator!=(const AtomMatch &other) const;

        QString toString() const;

        int nGroups() const;

        Selector<Atom> group(int i) const;

        QList<Selector<Atom>> groups() const;

    protected:
        Selector<Atom> reference;
        QList<QList<qint64>> matches;
    };

    /** This class holds the result of performing a match on multiple
     *  molecules */
    class SIREMOL_EXPORT AtomMatchM : public SireBase::ConcreteProperty<AtomMatchM, SireBase::Property>
    {

        friend QDataStream & ::operator<<(QDataStream &, const AtomMatchM &);
        friend QDataStream & ::operator>>(QDataStream &, AtomMatchM &);

    public:
        AtomMatchM();
        AtomMatchM(const AtomMatch &match);
        AtomMatchM(const QList<AtomMatch> &matches);
        AtomMatchM(const SelectResult &mols);

        AtomMatchM(const AtomMatchM &other);
        virtual ~AtomMatchM();

        static const char *typeName();

        virtual const char *what() const
        {
            return AtomMatchM::typeName();
        }

        virtual AtomMatchM *clone() const
        {
            return new AtomMatchM(*this);
        }

        AtomMatchM &operator=(const AtomMatchM &other);

        bool operator==(const AtomMatchM &other) const;
        bool operator!=(const AtomMatchM &other) const;

        QString toString() const;

        bool isEmpty() const;

        int nGroups() const;

        Selector<Atom> group(int i) const;

        QList<Selector<Atom>> groups() const;

    protected:
        QList<AtomMatch> matches;
    };

} // namespace SireMM

Q_DECLARE_METATYPE(SireMol::AtomMatch)
Q_DECLARE_METATYPE(SireMol::AtomMatchM)

SIRE_EXPOSE_CLASS(SireMol::AtomMatch)
SIRE_EXPOSE_CLASS(SireMol::AtomMatchM)

SIRE_END_HEADER

#endif
