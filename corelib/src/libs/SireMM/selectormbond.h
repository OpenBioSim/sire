/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMM_SELECTORMBOND_H
#define SIREMM_SELECTORMBOND_H

#include "SireMol/selectorm.hpp"
#include "SireMol/selectormol.h"

#include "selectorbond.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class Atom;
    class Residue;
    class Chain;
    class Segment;
    class CutGroup;
    class Molecule;
    class MoleculeGroup;
    class Molecules;
    class SelectResult;

    class SelectorMol;
    class EvaluatorM;

} // namespace SireMol

namespace SireMM
{
    class SelectorMBond;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::SelectorMBond &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::SelectorMBond &);

namespace SireMM
{

    /** Multi-molecule selector for bonds */
    class SIREMM_EXPORT SelectorMBond : public SireBase::ConcreteProperty<SelectorMBond, SireBase::Property>
    {

        friend QDataStream & ::operator<<(QDataStream &, const SelectorMBond &);
        friend QDataStream & ::operator>>(QDataStream &, SelectorMBond &);

    public:
        typedef QList<SelectorBond>::const_iterator const_iterator;
        typedef QList<SelectorBond>::const_iterator iterator;

        SelectorMBond();
        SelectorMBond(const Bond &view);

        SelectorMBond(const SireMol::Molecules &mols, const SireBase::PropertyMap &map = SireBase::PropertyMap());
        SelectorMBond(const SireMol::MoleculeGroup &mols, const SireBase::PropertyMap &map = SireBase::PropertyMap());
        SelectorMBond(const SireMol::MolGroupsBase &mols, const SireBase::PropertyMap &map = SireBase::PropertyMap());
        SelectorMBond(const SireMol::SelectResult &mols, const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SelectorMBond(const SelectorBond &bonds);

        SelectorMBond(const SireMol::SelectorMol &mols, const SireBase::PropertyMap &map = SireBase::PropertyMap());

        template <class U>
        SelectorMBond(const SireMol::SelectorM<U> &other, const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SelectorMBond(const SelectorMBond &bonds, const SireBase::Slice &slice);
        SelectorMBond(const SelectorMBond &bonds, const QList<qint64> &idxs);

        SelectorMBond(const SireMol::SelectResult &mols, const SireMol::BondID &bond,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SelectorMBond(const SireMol::SelectorM<SireMol::Atom> &atoms,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SelectorMBond(const SireMol::SelectorM<SireMol::Atom> &atoms0, const SireMol::SelectorM<SireMol::Atom> &atoms1,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SelectorMBond(const SelectorMBond &other);

        virtual ~SelectorMBond();

        static const char *typeName();

        virtual SelectorMBond *clone() const
        {
            return new SelectorMBond(*this);
        }

        SelectorMBond &operator=(const SelectorMBond &other);

        bool operator==(const SelectorMBond &other) const;
        bool operator!=(const SelectorMBond &other) const;

        Bond operator[](int i) const;
        SelectorMBond operator[](const SireBase::Slice &slice) const;
        SelectorMBond operator[](const QList<qint64> &idxs) const;
        SelectorMBond operator[](const SireMol::BondID &id) const;

        Bond operator()(int i) const;
        SelectorMBond operator()(const SireBase::Slice &slice) const;
        SelectorMBond operator()(const QList<qint64> &idxs) const;
        SelectorMBond operator()(const SireMol::BondID &id) const;

        QList<SireMol::MolViewPtr> toList() const;
        SireMol::Molecules toMolecules() const;

        int count() const;
        int size() const;

        void update(const SireMol::MoleculeView &molecule);
        void update(const SireMol::MoleculeData &moldata);
        void update(const SireMol::Molecules &molecules);
        void update(const SireMol::SelectorMol &molecules);

        SireMol::EvaluatorM evaluate() const;

        SireMol::SelectorMol extract() const;

        bool isSelector() const;

        SireMol::MoleculeGroup toMoleculeGroup() const;
        SireMol::SelectResult toSelectResult() const;

        SireMol::Molecule molecule(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Molecule molecule(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Molecule molecule(const SireMol::MolID &molid, const SireBase::PropertyMap &map = SireBase::PropertyMap());

        SireMol::SelectorMol molecules() const;
        SireMol::SelectorMol molecules(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorMol molecules(const SireBase::Slice &slice,
                                       const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorMol molecules(const QList<qint64> &idxs,
                                       const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorMol molecules(const QString &name,
                                       const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorMol molecules(const SireMol::MolID &molid,
                                       const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::Atom atom(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Atom atom(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Atom atom(const SireMol::AtomID &atomid, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::Residue residue(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Residue residue(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Residue residue(const SireMol::ResID &resid,
                                 const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::Chain chain(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Chain chain(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Chain chain(const SireMol::ChainID &chainid,
                             const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::Segment segment(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Segment segment(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::Segment segment(const SireMol::SegID &segid,
                                 const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::CutGroup cutGroup(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::CutGroup cutGroup(const QString &name, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::CutGroup cutGroup(const SireMol::CGID &cgid,
                                   const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::SelectorM<SireMol::Atom> atoms() const;
        SireMol::SelectorM<SireMol::Atom> atoms(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Atom> atoms(const SireBase::Slice &slice,
                                                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Atom> atoms(const QList<qint64> &idxs,
                                                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Atom> atoms(const QString &name,
                                                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Atom> atoms(const SireMol::AtomID &atomid,
                                                const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::SelectorM<SireMol::Residue> residues() const;
        SireMol::SelectorM<SireMol::Residue> residues(int i,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Residue> residues(const SireBase::Slice &slice,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Residue> residues(const QList<qint64> &idxs,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Residue> residues(const QString &name,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Residue> residues(const SireMol::ResID &resid,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::SelectorM<SireMol::Chain> chains() const;
        SireMol::SelectorM<SireMol::Chain> chains(int i, const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Chain> chains(const SireBase::Slice &slice,
                                                  const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Chain> chains(const QList<qint64> &idxs,
                                                  const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Chain> chains(const QString &name,
                                                  const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Chain> chains(const SireMol::ChainID &chainid,
                                                  const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::SelectorM<SireMol::Segment> segments() const;
        SireMol::SelectorM<SireMol::Segment> segments(int i,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Segment> segments(const SireBase::Slice &slice,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Segment> segments(const QList<qint64> &idxs,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Segment> segments(const QString &name,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::Segment> segments(const SireMol::SegID &segid,
                                                      const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SireMol::SelectorM<SireMol::CutGroup> cutGroups() const;
        SireMol::SelectorM<SireMol::CutGroup> cutGroups(int i,
                                                        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::CutGroup> cutGroups(const SireBase::Slice &slice,
                                                        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::CutGroup> cutGroups(const QList<qint64> &idxs,
                                                        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::CutGroup> cutGroups(const QString &name,
                                                        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;
        SireMol::SelectorM<SireMol::CutGroup> cutGroups(const SireMol::CGID &cgid,
                                                        const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

        SelectorMBond add(const SelectorMBond &other) const;

        SelectorMBond intersection(const SelectorMBond &other) const;

        SelectorMBond invert(const SireBase::PropertyMap &map) const;
        SelectorMBond invert() const;

        bool hasProperty(const SireBase::PropertyName &key) const;
        bool hasMetadata(const SireBase::PropertyName &key) const;
        bool hasMetadata(const SireBase::PropertyName &key, const SireBase::PropertyName &metakey) const;

        QStringList propertyKeys() const;
        QStringList metadataKeys() const;
        QStringList metadataKeys(const SireBase::PropertyName &key) const;

        QList<SireBase::Properties> properties() const;

        QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key) const;
        QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key, const Property &default_value) const;

        QList<SireUnits::Dimension::Length> lengths() const;
        QList<SireUnits::Dimension::Length> lengths(const SireBase::PropertyMap &map) const;

        QList<SireUnits::Dimension::Length> measures() const;
        QList<SireUnits::Dimension::Length> measures(const SireBase::PropertyMap &map) const;

        QList<SireCAS::Expression> potentials() const;
        QList<SireCAS::Expression> potentials(const SireBase::PropertyMap &map) const;

        QList<SireUnits::Dimension::GeneralUnit> energies() const;
        QList<SireUnits::Dimension::GeneralUnit> energies(const SireBase::PropertyMap &map) const;

        SireUnits::Dimension::GeneralUnit energy() const;
        SireUnits::Dimension::GeneralUnit energy(const SireBase::PropertyMap &map) const;

        SireMol::SelectResult search(const QString &search_string) const;
        SireMol::SelectResult search(const QString &search_string, const SireBase::PropertyMap &map) const;

        QList<SireMol::BondID> IDs() const;

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
        void saveFrame(int frame);
        void saveFrame();
        void deleteFrame(int frame);
        void deleteAllFrames();

        void loadFrame(int frame, const SireBase::PropertyMap &map);
        void saveFrame(int frame, const SireBase::PropertyMap &map);
        void saveFrame(const SireBase::PropertyMap &map);
        void deleteFrame(int frame, const SireBase::PropertyMap &map);
        void deleteAllFrames(const SireBase::PropertyMap &map);

        const_iterator begin() const;
        const_iterator end() const;

        const_iterator constBegin() const;
        const_iterator constEnd() const;

        virtual QString toString() const;

    protected:
        void _append(const Bond &bond);
        void _append(const SelectorBond &bonds);

        /** The actual bonds */
        QList<SelectorBond> bnds;
    };

} // end of namespace SireMM

#include "SireMol/atom.h"
#include "SireMol/chain.h"
#include "SireMol/cutgroup.h"
#include "SireMol/residue.h"
#include "SireMol/segment.h"
#include "SireMol/selectorm.hpp"

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

namespace SireMM
{

    template <class T>
    SIRE_OUTOFLINE_TEMPLATE SelectorMBond::SelectorMBond(const SireMol::SelectorM<T> &other,
                                                         const SireBase::PropertyMap &map)
        : SireBase::ConcreteProperty<SelectorMBond, SireBase::Property>()
    {
        if (not other.isEmpty())
        {
            this->bnds.reserve(other.nMolecules());

            for (int i = 0; i < other.nMolecules(); ++i)
            {
                this->bnds.append(SelectorBond(other.molecule(i), map));
            }
        }
    }

} // namespace SireMM

#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireMM::SelectorMBond)

SIRE_EXPOSE_CLASS(SireMM::SelectorMBond)

SIRE_END_HEADER

#endif
