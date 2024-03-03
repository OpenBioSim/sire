/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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

#ifndef SIREMOL_ATOMIDXMAPPING_H
#define SIREMOL_ATOMIDXMAPPING_H

#include "SireBase/property.h"

#include "SireMol/atomidx.h"
#include "SireMol/cgatomidx.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class AtomIdxMapping;
    class AtomIdxMappingEntry;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::AtomIdxMapping &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::AtomIdxMapping &);

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::AtomIdxMappingEntry &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::AtomIdxMappingEntry &);

namespace SireMol
{
    /** This is an individual mapping for a single atom */
    class SIREMOL_EXPORT AtomIdxMappingEntry
    {
        friend QDataStream & ::operator<<(QDataStream &, const AtomIdxMappingEntry &);
        friend QDataStream & ::operator>>(QDataStream &, AtomIdxMappingEntry &);

    public:
        AtomIdxMappingEntry();
        AtomIdxMappingEntry(const AtomIdx &index0, const AtomIdx &index1,
                            const MoleculeInfoData &molinfo0,
                            const MoleculeInfoData &molinfo1,
                            bool is_unmapped_in_reference);

        AtomIdxMappingEntry(const AtomIdxMappingEntry &other);
        ~AtomIdxMappingEntry();

        AtomIdxMappingEntry &operator=(const AtomIdxMappingEntry &other);

        bool operator==(const AtomIdxMappingEntry &other) const;
        bool operator!=(const AtomIdxMappingEntry &other) const;

        AtomIdxMappingEntry *clone() const;

        bool isNull() const;

        QString toString() const;

        static const char *typeName();
        const char *what() const;

        bool isUnmappedIn0() const;
        bool isUnmappedIn1() const;

        bool isMappedIn0() const;
        bool isMappedIn1() const;

        AtomIdx atomIdx0() const;
        AtomIdx atomIdx1() const;

        CGAtomIdx cgAtomIdx0() const;
        CGAtomIdx cgAtomIdx1() const;

    private:
        /** The atomidx in the reference state */
        AtomIdx atomidx0;

        /** The atomidx in the perturbed state */
        AtomIdx atomidx1;

        /** The CGAtomIdx in the reference state */
        CGAtomIdx cgatomidx0;

        /** The CGAtomIdx in the perturbed state */
        CGAtomIdx cgatomidx1;

        /** Whether or not this atom is unmapped in the reference state */
        bool unmapped0;
    };

    /** This class holds the mapping from one set of atom indices to another.
     *  This enables you to associate, atom by atom, atom indices in one set to
     *  atom indices in another set. This is useful, e.g. for building perturbations,
     *  or for specifying mappings for alignments or RMSD calculations etc.
     *
     *  This is mainly designed to provide sufficient information to merge
     *  properties together. It lists not just the atoms that map, but
     *  also which atoms are the ghost atoms that map (i.e. were created
     *  as ghost equivalents)
     */
    class SIREMOL_EXPORT AtomIdxMapping
        : public SireBase::ConcreteProperty<AtomIdxMapping, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const AtomIdxMapping &);
        friend QDataStream & ::operator>>(QDataStream &, AtomIdxMapping &);

    public:
        AtomIdxMapping();
        AtomIdxMapping(const AtomIdxMappingEntry &entry);
        AtomIdxMapping(const QList<AtomIdxMappingEntry> &entries);

        AtomIdxMapping(const AtomIdxMapping &other);

        ~AtomIdxMapping();

        AtomIdxMapping &operator=(const AtomIdxMapping &other);

        bool operator==(const AtomIdxMapping &other) const;
        bool operator!=(const AtomIdxMapping &other) const;

        AtomIdxMapping *clone() const;

        QString toString() const;

        static const char *typeName();
        const char *what() const;

        AtomIdxMapping &operator+=(const AtomIdxMapping &other);
        AtomIdxMapping &operator+=(const AtomIdxMappingEntry &entry);

        AtomIdxMapping operator+(const AtomIdxMapping &other) const;
        AtomIdxMapping operator+(const AtomIdxMappingEntry &entry) const;

        QList<AtomIdxMappingEntry>::const_iterator begin() const;
        QList<AtomIdxMappingEntry>::const_iterator constBegin() const;

        QList<AtomIdxMappingEntry>::const_iterator end() const;
        QList<AtomIdxMappingEntry>::const_iterator constEnd() const;

        QList<AtomIdxMappingEntry>::const_iterator find(AtomIdx atom) const;
        QList<AtomIdxMappingEntry>::const_iterator find(CGAtomIdx atom) const;

        void append(const AtomIdxMappingEntry &entry);
        void append(const AtomIdxMapping &other);

        void clear();

        bool isEmpty() const;

        int size() const;
        int count() const;

        const AtomIdxMappingEntry &operator[](int i) const;

        AtomIdxMappingEntry take(int i);
        AtomIdxMappingEntry take(const AtomIdx &atom);
        AtomIdxMappingEntry take(const CGAtomIdx &atom);

        void remove(int i);
        void remove(const AtomIdx &atom);
        void remove(const CGAtomIdx &atom);

        QList<AtomIdx> unmappedIn0() const;
        QList<AtomIdx> unmappedIn1() const;

        QList<AtomIdx> mappedIn0() const;
        QList<AtomIdx> mappedIn1() const;

        bool isUnmappedIn0(const AtomIdx &atom) const;
        bool isUnmappedIn1(const AtomIdx &atom) const;

        QHash<AtomIdx, AtomIdx> map0to1(bool include_unmapped = false) const;
        QHash<AtomIdx, AtomIdx> map1to0(bool include_unmapped = false) const;

    private:
        void assertSane() const;
        void rebuild();

        /** The list of atom entries */
        QList<AtomIdxMappingEntry> entries;

        /** The set of atoms that are unmapped in the reference state */
        QSet<AtomIdx> unmapped0_set;

        /** The set of atoms that are unmapped in the perturbed state */
        QSet<AtomIdx> unmapped1_set;

        /** The list of atoms that are unmapped in the reference state */
        QList<AtomIdx> unmapped0_list;

        /** The list of atoms that are unmapped in the perturbed state */
        QList<AtomIdx> unmapped1_list;

        /** The list of atoms that are mapped in the reference state */
        QList<AtomIdx> mapped0_list;

        /** The list of atoms that are mapped in the perturbed state */
        QList<AtomIdx> mapped1_list;

        /** The mapping from reference to perturbed atoms,
         *  which includes the unmapped atoms
         */
        QHash<AtomIdx, AtomIdx> map0_to_1_inc;

        /** The mapping from perturbed to reference atoms,
         *  which includes the unmapped atoms
         */
        QHash<AtomIdx, AtomIdx> map1_to_0_inc;

        /** The mapping from reference to perturbed atoms,
         *  which excludes the unmapped atoms
         */
        QHash<AtomIdx, AtomIdx> map0_to_1_exc;

        /** The mapping from perturbed to reference atoms,
         *  which excludes the unmapped atoms
         */
        QHash<AtomIdx, AtomIdx> map1_to_0_exc;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

    /** Return whether or not the passed atom is unmapped in the
     *  reference state - the atom index should be for the
     *  merged molecule
     */
    inline bool AtomIdxMapping::isUnmappedIn0(const AtomIdx &atom) const
    {
        return unmapped0_set.contains(atom);
    }

    /** Return whether or not the passed atom is unmapped in the
     *  perturbed state - the atom index should be for the
     *  merged molecule
     */
    inline bool AtomIdxMapping::isUnmappedIn1(const AtomIdx &atom) const
    {
        return unmapped1_set.contains(atom);
    }

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::AtomIdxMappingEntry)
Q_DECLARE_METATYPE(SireMol::AtomIdxMapping)

SIRE_EXPOSE_CLASS(SireMol::AtomIdxMappingEntry)
SIRE_EXPOSE_CLASS(SireMol::AtomIdxMapping)

SIRE_END_HEADER

#endif
