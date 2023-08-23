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

        operator Selector<Atom>() const;

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
        typedef typename QList<AtomMatch>::const_iterator iterator;
        typedef typename QList<AtomMatch>::const_iterator const_iterator;

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

        Atom operator[](int i) const;
        SelectorM<Atom> operator[](const SireBase::Slice &slice) const;
        SelectorM<Atom> operator[](const QList<qint64> &idxs) const;
        Atom operator[](const QString &name) const;
        Atom operator[](const typename Atom::ID &id) const;

        Atom operator()(int i) const;
        Atom operator()(const QString &name) const;
        Atom operator()(const typename Atom::ID &id) const;

        QList<MolViewPtr> toList() const;
        Molecules toMolecules() const;

        int count() const;
        int size() const;

        void update(const MoleculeView &molview);
        void update(const MoleculeData &moldata);
        void update(const Molecules &molecules);
        void update(const SelectorMol &molecules);

        EvaluatorM evaluate() const;

        MoleculeGroup toMoleculeGroup() const;
        SelectResult toSelectResult() const;

        bool isSelector() const;

        SelectorMol extract() const;

        QList<qint64> find(const Atom &view) const;
        QList<qint64> find(const Selector<Atom> &views) const;
        QList<qint64> find(const SelectorM<Atom> &views) const;

        SelectorM<Atom> intersection(const SelectorM<Atom> &other) const;
        SelectorM<Atom> intersection(const Selector<Atom> &views) const;
        SelectorM<Atom> intersection(const Atom &view) const;

        SelectorM<Atom> invert() const;

        bool intersects(const SelectorM<Atom> &other) const;
        bool intersects(const Selector<Atom> &view) const;
        bool intersects(const Atom &view) const;

        bool contains(const SelectorM<Atom> &other) const;
        bool contains(const Selector<Atom> &view) const;
        bool contains(const Atom &view) const;

        Molecule molecule(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Molecule molecule(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Molecule molecule(const MolID &molid, const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SelectorMol molecules() const;
        SelectorMol molecules(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorMol molecules(const SireBase::Slice &slice,
                              const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorMol molecules(const QList<qint64> &idxs, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorMol molecules(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorMol molecules(const MolID &molid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        Atom atom(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Atom atom(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Atom atom(const AtomID &atomid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        Residue residue(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Residue residue(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Residue residue(const ResID &resid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        Chain chain(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Chain chain(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Chain chain(const ChainID &chainid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        Segment segment(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Segment segment(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        Segment segment(const SegID &segid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        CutGroup cutGroup(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        CutGroup cutGroup(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        CutGroup cutGroup(const CGID &cgid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SelectorM<Atom> atoms() const;
        SelectorM<Atom> atoms(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Atom> atoms(const SireBase::Slice &slice,
                              const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Atom> atoms(const QList<qint64> &idxs, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Atom> atoms(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Atom> atoms(const AtomID &atomid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SelectorM<Residue> residues() const;
        SelectorM<Residue> residues(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Residue> residues(const SireBase::Slice &slice,
                                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Residue> residues(const QList<qint64> &idxs,
                                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Residue> residues(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Residue> residues(const ResID &resid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SelectorM<Chain> chains() const;
        SelectorM<Chain> chains(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Chain> chains(const SireBase::Slice &slice,
                                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Chain> chains(const QList<qint64> &idxs,
                                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Chain> chains(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Chain> chains(const ChainID &chainid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SelectorM<Segment> segments() const;
        SelectorM<Segment> segments(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Segment> segments(const SireBase::Slice &slice,
                                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Segment> segments(const QList<qint64> &idxs,
                                    const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Segment> segments(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<Segment> segments(const SegID &segid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SelectorM<CutGroup> cutGroups() const;
        SelectorM<CutGroup> cutGroups(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<CutGroup> cutGroups(const SireBase::Slice &slice,
                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<CutGroup> cutGroups(const QList<qint64> &idxs,
                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<CutGroup> cutGroups(const QString &name,
                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SelectorM<CutGroup> cutGroups(const CGID &cgid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SelectResult search(const QString &search_string) const;
        SelectResult search(const QString &search_string, const SireBase::PropertyMap &map) const;

        QList<Atom::Index> IDs() const;
        QList<Atom::Index> indexes() const;
        QList<Atom::Number> numbers() const;
        QList<Atom::Name> names() const;

        int nAtoms() const;
        int nResidues() const;
        int nChains() const;
        int nSegments() const;
        int nCutGroups() const;
        int nMolecules() const;

        bool isEmpty() const;

        int nFrames() const;
        int nFrames(const SireBase::PropertyMap &map) const;

        void loadFrame(int frame);
        void loadFrame(int frame, const SireBase::LazyEvaluator &evaluator);
        void saveFrame(int frame);
        void saveFrame();
        void deleteFrame(int frame);
        void deleteAllFrames();

        void loadFrame(int frame, const SireBase::PropertyMap &map);
        void loadFrame(int frame, const SireBase::LazyEvaluator &evaluator,
                       const SireBase::PropertyMap &map);
        void saveFrame(int frame, const SireBase::PropertyMap &map);
        void saveFrame(const SireBase::PropertyMap &map);
        void deleteFrame(int frame, const SireBase::PropertyMap &map);
        void deleteAllFrames(const SireBase::PropertyMap &map);

        const_iterator begin() const;
        const_iterator end() const;

        const_iterator constBegin() const;
        const_iterator constEnd() const;

        virtual QString toString() const;

        bool hasProperty(const PropertyName &key) const;
        bool hasMetadata(const PropertyName &metakey) const;
        bool hasMetadata(const PropertyName &key, const PropertyName &metakey) const;

        QStringList propertyKeys() const;
        QStringList metadataKeys() const;
        QStringList metadataKeys(const PropertyName &key) const;

        template <class V>
        QList<V> property(const PropertyName &key) const;

        template <class V>
        QList<V> metadata(const PropertyName &metakey) const;

        template <class V>
        QList<V> metadata(const PropertyName &key, const PropertyName &metakey) const;

        int nGroups() const;

        Selector<Atom> group(int i) const;

        QList<Selector<Atom>> groups() const;

        operator SelectorM<Atom>() const;

    protected:
        SelectorM<Atom> _atoms() const;

        QList<AtomMatch> matches;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

#ifndef GCCXML_PARSE
    template <class V>
    SIRE_OUTOFLINE_TEMPLATE QList<V> AtomMatchM::property(const PropertyName &key) const
    {
        return this->_atoms().property<V>(key);
    }

    template <class V>
    SIRE_OUTOFLINE_TEMPLATE QList<V> AtomMatchM::metadata(const PropertyName &metakey) const
    {
        return this->_atoms().metadata<V>(metakey);
    }

    template <class V>
    SIRE_OUTOFLINE_TEMPLATE QList<V> AtomMatchM::metadata(const PropertyName &key,
                                                          const PropertyName &metakey) const
    {
        return this->_atoms().metadata<V>(key, metakey);
    }

#endif // GCCXML_PARSE

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // namespace SireMM

Q_DECLARE_METATYPE(SireMol::AtomMatch)
Q_DECLARE_METATYPE(SireMol::AtomMatchM)

SIRE_EXPOSE_CLASS(SireMol::AtomMatch)
SIRE_EXPOSE_CLASS(SireMol::AtomMatchM)

SIRE_END_HEADER

#endif
