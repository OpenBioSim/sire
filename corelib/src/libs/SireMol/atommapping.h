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

#ifndef SIREMOL_ATOMMAPPING_H
#define SIREMOL_ATOMMAPPING_H

#include "core.h"
#include "atom.h"
#include "selectorm.hpp"

#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class AtomMapping;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::AtomMapping &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::AtomMapping &);

namespace SireMol
{

    /** This class holds the mapping from one set of atoms to another.
     *  This enables you associate, atom by atom, atoms in one set to
     *  atoms in another set. This is useful, e.g. for building perturbations,
     *  or for specifying mappings for alignments or RMSD calculations etc.
     */
    class SIREMOL_EXPORT AtomMapping
        : public SireBase::ConcreteProperty<AtomMapping, SireBase::Property>
    {

        friend QDataStream & ::operator<<(QDataStream &, const AtomMapping &);
        friend QDataStream & ::operator>>(QDataStream &, AtomMapping &);

    public:
        AtomMapping();
        AtomMapping(const SelectorM<Atom> &atoms0,
                    const SelectorM<Atom> &atoms1,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap());

        AtomMapping(const SelectorM<Atom> &atoms0,
                    const SelectorM<Atom> &atoms1,
                    const SireBase::PropertyMap &map0,
                    const SireBase::PropertyMap &map1);

        AtomMapping(const MoleculeView &mol0,
                    const MoleculeView &mol1,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap());

        AtomMapping(const MoleculeView &mol0,
                    const MoleculeView &mol1,
                    const SireBase::PropertyMap &map0,
                    const SireBase::PropertyMap &map1);

        AtomMapping(const SelectorMol &mols0,
                    const SelectorMol &mols1,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap());

        AtomMapping(const SelectorMol &mols0,
                    const SelectorMol &mols1,
                    const SireBase::PropertyMap &map0,
                    const SireBase::PropertyMap &map1);

        AtomMapping(const SelectorM<Atom> &atoms0,
                    const SelectorM<Atom> &atoms1,
                    const SelectorM<Atom> &matched_atoms0,
                    const SelectorM<Atom> &matched_atoms1,
                    const SireBase::PropertyMap &map = SireBase::PropertyMap());

        AtomMapping(const SelectorM<Atom> &atoms0,
                    const SelectorM<Atom> &atoms1,
                    const SelectorM<Atom> &matched_atoms0,
                    const SelectorM<Atom> &matched_atoms1,
                    const SireBase::PropertyMap &map0,
                    const SireBase::PropertyMap &map1);

        AtomMapping(const AtomMapping &other);

        ~AtomMapping();

        AtomMapping &operator=(const AtomMapping &other);

        bool operator==(const AtomMapping &other) const;
        bool operator!=(const AtomMapping &other) const;

        Atom operator[](int i) const;
        Atom operator[](const Atom &atom) const;

        SelectorM<Atom> operator[](const SireBase::Slice &slice) const;
        SelectorM<Atom> operator[](const QList<qint64> &idxs) const;

        SelectorM<Atom> operator[](const Selector<Atom> &atoms) const;
        SelectorM<Atom> operator[](const SelectorM<Atom> &atoms) const;

        virtual AtomMapping *clone() const;

        const char *what() const;
        static const char *typeName();

        QString toString() const;

        int count() const;
        int size() const;

        bool isEmpty() const;

        bool contains(const Atom &atom) const;

        const SelectorM<Atom> &atoms0() const;
        const SelectorM<Atom> &atoms1() const;

        SelectorM<Atom> mappedAtoms0() const;
        SelectorM<Atom> mappedAtoms1() const;

        SelectorM<Atom> unmappedAtoms0() const;
        SelectorM<Atom> unmappedAtoms1() const;

        bool isMapped(const Atom &atom) const;

        bool isSingleMolecule() const;

        void assertSingleMolecule() const;

        AtomMapping swap() const;

        AtomMapping align() const;

        AtomMapping alignTo0() const;
        AtomMapping alignTo1() const;

        Atom map(const Atom &atom, bool find_all = true) const;

        SelectorM<Atom> map(const MoleculeView &atoms,
                            bool find_all = true) const;
        SelectorM<Atom> map(const SelectorM<Atom> &atoms,
                            bool find_all = true) const;

        Atom find(const Atom &atom, const MoleculeView &container,
                  bool find_all = true) const;
        Atom find(const Atom &atom, const SelectorM<Atom> &container,
                  bool find_all = true) const;

        Selector<Atom> find(const MoleculeView &atoms,
                            const MoleculeView &container,
                            bool find_all = true) const;

        SelectorM<Atom> find(const SelectorM<Atom> &atoms,
                             const MoleculeView &container,
                             bool find_all = true) const;

        SelectorM<Atom> find(const MoleculeView &atoms,
                             const SelectorM<Atom> &container,
                             bool find_all = true) const;

        SelectorM<Atom> find(const SelectorM<Atom> &atoms,
                             const SelectorM<Atom> &container,
                             bool find_all = true) const;

        const SireBase::PropertyMap &propertyMap0() const;
        const SireBase::PropertyMap &propertyMap1() const;

    private:
        /** The reference atoms - we map from these to the other atoms */
        SelectorM<Atom> atms0;

        /** The mapped atoms - we map from the reference atoms to these atoms */
        SelectorM<Atom> atms1;

        /** The original set of reference atoms - this includes the atoms
         *  in the same order as input, including any unmatched atoms
         */
        SelectorM<Atom> orig_atms0;

        /** The original set of mapped atoms - this includes the atoms
         *  in the same order as input, including any unmatched atoms
         */
        SelectorM<Atom> orig_atms1;

        /** The indexes of reference atoms in `orig_atms0` that are not mapped,
         *  in the order they appear in `orig_atms0`
         */
        QList<qint64> unmapped_atms0;

        /** The indexes of mapped atoms in `orig_atms1` that are not mapped,
         *  in the order they appear in `orig_atms1`
         */
        QList<qint64> unmapped_atms1;

        /** Property map for the reference atoms */
        SireBase::PropertyMap map0;

        /** Property map for the mapped atoms */
        SireBase::PropertyMap map1;
    };
}

Q_DECLARE_METATYPE(SireMol::AtomMapping)

SIRE_EXPOSE_CLASS(SireMol::AtomMapping)

SIRE_END_HEADER

#endif
