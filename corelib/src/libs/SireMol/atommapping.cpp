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

#include "atommapping.h"

#include "SireMaths/align.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<AtomMapping> r_mapping;

QDataStream &operator<<(QDataStream &ds, const AtomMapping &mapping)
{
    writeHeader(ds, r_mapping, 2);

    SharedDataStream sds(ds);

    sds << mapping.atms0 << mapping.atms1
        << mapping.orig_atms0 << mapping.orig_atms1
        << mapping.unmapped_atms0 << mapping.unmapped_atms1
        << mapping.map0 << mapping.map1
        << static_cast<const Property &>(mapping);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AtomMapping &mapping)
{
    VersionID v = readHeader(ds, r_mapping);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> mapping.atms0 >> mapping.atms1 >> mapping.orig_atms0 >> mapping.orig_atms1 >> mapping.unmapped_atms0 >> mapping.unmapped_atms1 >> mapping.map0 >> mapping.map1 >> static_cast<Property &>(mapping);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mapping.atms0 >> mapping.atms1 >> static_cast<Property &>(mapping);

        mapping.map0 = PropertyMap();
        mapping.map1 = PropertyMap();

        mapping.orig_atms0 = mapping.atms0;
        mapping.orig_atms1 = mapping.atms1;

        mapping.unmapped_atms0.clear();
        mapping.unmapped_atms1.clear();
    }
    else
        throw version_error(v, "1,2", r_mapping, CODELOC);

    return ds;
}

AtomMapping::AtomMapping() : ConcreteProperty<AtomMapping, Property>()
{
}

AtomMapping::AtomMapping(const SelectorM<Atom> &atoms0,
                         const SelectorM<Atom> &atoms1,
                         const PropertyMap &m0,
                         const PropertyMap &m1)
    : ConcreteProperty<AtomMapping, Property>(),
      atms0(atoms0), atms1(atoms1),
      orig_atms0(atoms0), orig_atms1(atoms1),
      map0(m0), map1(m1)
{
    if (atms0.count() != atms1.count())
    {
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of atoms in 'atoms0' (%1) is not equal to the "
                                                "number of atoms in 'atoms1' (%2). You can only create a "
                                                "mapping if the number of atoms is the same.")
                                                .arg(atms0.count())
                                                .arg(atms1.count()),
                                            CODELOC);
    }
}

AtomMapping::AtomMapping(const SelectorM<Atom> &atoms0,
                         const SelectorM<Atom> &atoms1,
                         const PropertyMap &map)
    : ConcreteProperty<AtomMapping, Property>()
{
    this->operator=(AtomMapping(atoms0, atoms1, map, map));
}

AtomMapping::AtomMapping(const MoleculeView &mol0, const MoleculeView &mol1,
                         const PropertyMap &map)
    : ConcreteProperty<AtomMapping, Property>()
{
    this->operator=(AtomMapping(SelectorM<Atom>(mol0.atoms()),
                                SelectorM<Atom>(mol1.atoms()),
                                map));
}

AtomMapping::AtomMapping(const MoleculeView &mol0, const MoleculeView &mol1,
                         const PropertyMap &m0, const PropertyMap &m1)
    : ConcreteProperty<AtomMapping, Property>()
{
    this->operator=(AtomMapping(SelectorM<Atom>(mol0.atoms()),
                                SelectorM<Atom>(mol1.atoms()),
                                m0, m1));
}

AtomMapping::AtomMapping(const SelectorMol &mols0, const SelectorMol &mols1,
                         const PropertyMap &map)
    : ConcreteProperty<AtomMapping, Property>()
{
    this->operator=(AtomMapping(mols0.atoms(), mols1.atoms(), map));
}

AtomMapping::AtomMapping(const SelectorMol &mols0, const SelectorMol &mols1,
                         const PropertyMap &m0, const PropertyMap &m1)
    : ConcreteProperty<AtomMapping, Property>()
{
    this->operator=(AtomMapping(mols0.atoms(), mols1.atoms(), m0, m1));
}

AtomMapping::AtomMapping(const SelectorM<Atom> &atoms0,
                         const SelectorM<Atom> &atoms1,
                         const SelectorM<Atom> &matched_atoms0,
                         const SelectorM<Atom> &matched_atoms1,
                         const SireBase::PropertyMap &m0,
                         const SireBase::PropertyMap &m1)
    : ConcreteProperty<AtomMapping, Property>(),
      atms0(matched_atoms0), atms1(matched_atoms1),
      orig_atms0(atoms0), orig_atms1(atoms1),
      map0(m0), map1(m1)
{
    if (atms0.count() != atms1.count())
    {
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of atoms in 'atoms0' (%1) is not equal to the "
                                                "number of atoms in 'atoms1' (%2). You can only create a "
                                                "mapping if the number of atoms is the same.")
                                                .arg(atms0.count())
                                                .arg(atms1.count()),
                                            CODELOC);
    }

    // find the indexes of missing atoms
    for (int i = 0; i < atoms0.count(); ++i)
    {
        if (not atms0.contains(atoms0[i]))
            this->unmapped_atms0.append(i);
    }

    for (int i = 0; i < atoms1.count(); ++i)
    {
        if (not atms1.contains(atoms1[i]))
            this->unmapped_atms1.append(i);
    }
}

AtomMapping::AtomMapping(const SelectorM<Atom> &atoms0,
                         const SelectorM<Atom> &atoms1,
                         const SelectorM<Atom> &matched_atoms0,
                         const SelectorM<Atom> &matched_atoms1,
                         const SireBase::PropertyMap &map)
    : ConcreteProperty<AtomMapping, Property>()
{
    this->operator=(AtomMapping(atoms0, atoms1,
                                matched_atoms0, matched_atoms1,
                                map, map));
}

AtomMapping::AtomMapping(const AtomMapping &other)
    : ConcreteProperty<AtomMapping, Property>(other),
      atms0(other.atms0), atms1(other.atms1),
      orig_atms0(other.orig_atms0), orig_atms1(other.orig_atms1),
      unmapped_atms0(other.unmapped_atms0), unmapped_atms1(other.unmapped_atms1),
      map0(other.map0), map1(other.map1)
{
}

AtomMapping::~AtomMapping()
{
}

AtomMapping &AtomMapping::operator=(const AtomMapping &other)
{
    if (this != &other)
    {
        atms0 = other.atms0;
        atms1 = other.atms1;
        orig_atms0 = other.orig_atms0;
        orig_atms1 = other.orig_atms1;
        unmapped_atms0 = other.unmapped_atms0;
        unmapped_atms1 = other.unmapped_atms1;
        map0 = other.map0;
        map1 = other.map1;

        Property::operator=(other);
    }

    return *this;
}

bool AtomMapping::operator==(const AtomMapping &other) const
{
    return orig_atms0 == other.orig_atms0 and orig_atms1 == other.orig_atms1 and
           atms0 == other.atms0 and atms1 == other.atms1 and
           map0 == other.map0 and map1 == other.map1;
}

bool AtomMapping::operator!=(const AtomMapping &other) const
{
    return not this->operator==(other);
}

Atom AtomMapping::operator[](int i) const
{
    return this->atms1[i];
}

Atom AtomMapping::operator[](const Atom &atom) const
{
    return this->map(atom, false);
}

SelectorM<Atom> AtomMapping::operator[](const SireBase::Slice &slice) const
{
    return this->atms1[slice];
}

SelectorM<Atom> AtomMapping::operator[](const QList<qint64> &idxs) const
{
    return this->atms1[idxs];
}

SelectorM<Atom> AtomMapping::operator[](const Selector<Atom> &atoms) const
{
    return this->operator[](SelectorM<Atom>(atoms));
}

SelectorM<Atom> AtomMapping::operator[](const SelectorM<Atom> &atoms) const
{
    return this->map(atoms, false);
}

AtomMapping *AtomMapping::clone() const
{
    return new AtomMapping(*this);
}

const char *AtomMapping::what() const
{
    return AtomMapping::typeName();
}

const char *AtomMapping::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AtomMapping>());
}

int AtomMapping::size() const
{
    return this->atms0.count();
}

int AtomMapping::count() const
{
    return this->size();
}

bool AtomMapping::isEmpty() const
{
    return this->atms0.isEmpty();
}

QString AtomMapping::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("AtomMapping::empty");
    }
    else
    {
        QStringList parts;

        const auto n = this->count();

        auto atom_to_string = [](const Atom &atom)
        {
            return QString("%1 Atom( %2:%3 )").arg(atom.data().number().toString()).arg(atom.name().value()).arg(atom.number().value());
        };

        if (n <= 10)
        {
            for (int i = 0; i < n; ++i)
            {
                const auto atom0 = this->atms0[i];
                const auto atom1 = this->atms1[i];

                parts.append(QString("%1: %2 <=> %3").arg(i).arg(atom_to_string(atom0)).arg(atom_to_string(atom1)));
            }
        }
        else
        {
            for (int i = 0; i < 5; ++i)
            {
                const auto atom0 = this->atms0[i];
                const auto atom1 = this->atms1[i];

                parts.append(QString("%1: %2 <=> %3").arg(i).arg(atom_to_string(atom0)).arg(atom_to_string(atom1)));
            }

            parts.append("...");

            for (int i = n - 5; i < n; ++i)
            {
                const auto atom0 = this->atms0[i];
                const auto atom1 = this->atms1[i];

                parts.append(QString("%1: %2 <=> %3").arg(i).arg(atom_to_string(atom0)).arg(atom_to_string(atom1)));
            }
        }

        return QObject::tr("AtomMapping( size=%1, unmapped0=%2, unmapped1=%3\n%4\n)")
            .arg(n)
            .arg(this->unmapped_atms0.count())
            .arg(this->unmapped_atms1.count())
            .arg(parts.join("\n"));
    }
}

/** Return the original reference atoms. This is the collection of both
 *  mapped and unmapped reference atoms
 */
const SelectorM<Atom> &AtomMapping::atoms0() const
{
    return this->orig_atms0;
}

/** Return the original mapped atoms. This is the collection of both
 *  mapped and unmapped mapped atoms
 */
const SelectorM<Atom> &AtomMapping::atoms1() const
{
    return this->orig_atms1;
}

/** Return an AtomMapping that swaps the reference and mapped atoms */
AtomMapping AtomMapping::swap() const
{
    AtomMapping ret;

    ret.orig_atms0 = this->orig_atms1;
    ret.orig_atms1 = this->orig_atms0;

    ret.atms0 = this->atms1;
    ret.atms1 = this->atms0;

    ret.map0 = this->map1;
    ret.map1 = this->map0;

    ret.unmapped_atms0 = this->unmapped_atms1;
    ret.unmapped_atms1 = this->unmapped_atms0;

    return ret;
}

/** Return whether or not the forward mapping contains the
 *  passed atom - this returns true if the atom is contained
 *  in the original reference atoms, i.e. it doesn't guarantee
 *  that the atom is mapped. Use the 'isMapped' method to
 *  check if the atom is mapped.
 */
bool AtomMapping::contains(const Atom &atom) const
{
    return this->orig_atms0.contains(atom);
}

/** Map from 'atom' (which must be in the reference atoms) to
 *  the corresponding atom in the mapped atoms */
Atom AtomMapping::map(const Atom &atom, bool find_all) const
{
    const auto idxs = this->atms0.find(atom);

    if (idxs.isEmpty())
    {
        if (find_all)
            throw SireMol::missing_atom(QObject::tr(
                                            "Cannot find '%1' in the mapping.")
                                            .arg(atom.toString()),
                                        CODELOC);

        return Atom();
    }

    return this->operator[](idxs[0]);
}

/** Map from the passed 'atoms' (which must all be in the reference
 *  atoms) to the corresponding atoms in the mapped atoms. The
 *  mapped atoms will be returned in the same order as the
 *  reference atoms appeared in 'atoms'. If 'find_all` is false
 *  then this will use null atoms in the map when the mapped
 *  atom cannot be found */
SelectorM<Atom> AtomMapping::map(const MoleculeView &atoms,
                                 bool find_all) const
{
    return this->map(SelectorM<Atom>(atoms.atoms()), find_all);
}

/** Map from the passed 'atoms' (which must all be in the reference
 *  atoms) to the corresponding atoms in the mapped atoms. The
 *  mapped atoms will be returned in the same order as the
 *  reference atoms appeared in 'atoms'. If 'find_all` is false
 *  then this will use null atoms in the map when the mapped
 *  atom cannot be found */
SelectorM<Atom> AtomMapping::map(const SelectorM<Atom> &atoms,
                                 bool find_all) const
{
    SelectorM<Atom> ret;

    for (int i = 0; i < atoms.count(); ++i)
    {
        const auto atom_idxs = this->atms0.find(atoms[i]);

        if (atom_idxs.isEmpty())
        {
            if (find_all)
                throw SireMol::missing_atom(QObject::tr(
                                                "Cannot find '%1' in the mapping.")
                                                .arg(atoms[i].toString()),
                                            CODELOC);

            ret += Atom();
        }
        else
        {
            ret += this->atms1[atom_idxs[0]];
        }
    }

    return ret;
}

/** Find and return the equivalent of 'atom' in the passed container.
 *  This maps 'atom' from the reference to the mapped atom, and then
 *  locates and returns the mapped atom from the container.
 */
Atom AtomMapping::find(const Atom &atom, const MoleculeView &container,
                       bool find_all) const
{
    return this->find(atom, SelectorM<Atom>(container.atoms()), find_all);
}

/** Find and return the equivalent of 'atom' in the passed container.
 *  This maps 'atom' from the reference to the mapped atom, and then
 *  locates and returns the mapped atom from the container.
 */
Atom AtomMapping::find(const Atom &atom,
                       const SelectorM<Atom> &container,
                       bool find_all) const
{
    const auto mapped = this->map(atom, find_all);

    if (mapped.isEmpty() and not find_all)
        return Atom();

    auto filtered = container.filter(mapped);

    if (filtered.isEmpty())
    {
        throw SireMol::missing_atom(QObject::tr(
                                        "There is no atom in %1 that matches atom %2, which was mapped "
                                        "to atom %3.")
                                        .arg(container.toString())
                                        .arg(atom.toString())
                                        .arg(mapped.toString()),
                                    CODELOC);
    }
    else if (filtered.count() != 1)
    {
        throw SireError::program_bug(QObject::tr(
                                         "There should only be a single match for the mapping: %1")
                                         .arg(filtered.toString()),
                                     CODELOC);
    }

    return filtered(0);
}

/** Find and return the equivalent of 'atom' in the passed container.
 *  This maps 'atom' from the reference to the mapped atom, and then
 *  locates and returns the mapped atom from the container.
 */
Selector<Atom> AtomMapping::find(const MoleculeView &atoms,
                                 const MoleculeView &container,
                                 bool find_all) const
{
    auto found = this->find(SelectorM<Atom>(atoms.atoms()),
                            container.atoms(), find_all);

    if (not found.isEmpty())
    {
        // this is definitely a single-molecule set
        return Selector<Atom>(found(0).data(), found.IDs());
    }
    else
        return Selector<Atom>();
}

/** Find and return the equivalent of 'atoms' in the passed container.
 *  This maps 'atoms' from the reference to the mapped atoms, and then
 *  locates and returns the mapped atoms from the container. Note that
 *  all atoms must be found, and they will be returned in the same
 *  order as 'atoms'
 */
SelectorM<Atom> AtomMapping::find(const SelectorM<Atom> &atoms,
                                  const MoleculeView &container,
                                  bool find_all) const
{
    return this->find(atoms, SelectorM<Atom>(container.atoms()), find_all);
}

/** Find and return the equivalent of 'atoms' in the passed container.
 *  This maps 'atoms' from the reference to the mapped atoms, and then
 *  locates and returns the mapped atoms from the container. Note that
 *  all atoms must be found, and they will be returned in the same
 *  order as 'atoms'
 */
SelectorM<Atom> AtomMapping::find(const MoleculeView &atoms,
                                  const SelectorM<Atom> &container,
                                  bool find_all) const
{
    return this->find(SelectorM<Atom>(atoms.atoms()), container, find_all);
}

/** Find and return the equivalent of 'atoms' in the passed container.
 *  This maps 'atoms' from the reference to the mapped atoms, and then
 *  locates and returns the mapped atoms from the container. Note that
 *  all atoms must be found, and they will be returned in the same
 *  order as 'atoms'
 */
SelectorM<Atom> AtomMapping::find(const SelectorM<Atom> &atoms,
                                  const SelectorM<Atom> &container,
                                  bool find_all) const
{
    if (atoms.isEmpty())
        return SelectorM<Atom>();

    else if (container.isEmpty() and not find_all)
        return SelectorM<Atom>();

    const auto mapped = this->map(atoms, find_all);

    // need this way round as we have to have the same
    // ordering as mapped
    auto filtered = mapped.filter(container);

    if (find_all and (filtered.count() < mapped.count()))
    {
        throw SireMol::missing_atom(QObject::tr(
                                        "There were some atoms in %1 that could not be located "
                                        "in the container %2, when mapped to %3.")
                                        .arg(atoms.toString())
                                        .arg(container.toString())
                                        .arg(mapped.toString()),
                                    CODELOC);
    }
    else if (filtered.count() > mapped.count())
    {
        throw SireError::program_bug(QObject::tr(
                                         "There should only be a single match per atom for the mapping: %1")
                                         .arg(filtered.toString()),
                                     CODELOC);
    }

    // need to update so we have the same versions
    // as in the container
    filtered.update(container.molecules());

    return filtered;
}

/** Return all of the reference atoms that have been mapped,
 *  in the same order as the mapped atoms they match with
 */
SelectorM<Atom> AtomMapping::mappedAtoms0() const
{
    return this->atms0;
}

/** Return all of the mapped atoms that have been mapped,
 *  in the same order as the reference atoms they match with
 */
SelectorM<Atom> AtomMapping::mappedAtoms1() const
{
    return this->atms1;
}

/** Return all of the reference atoms that haven't been mapped,
 *  in the same order as they appear in the original reference
 */
SelectorM<Atom> AtomMapping::unmappedAtoms0() const
{
    if (this->unmapped_atms0.isEmpty())
    {
        return SelectorM<Atom>();
    }
    else
    {
        return this->orig_atms0[this->unmapped_atms0];
    }
}

/** Return all of the mapped atoms that haven't been mapped,
 *  in the same order as they appear in the original mapped atoms
 */
SelectorM<Atom> AtomMapping::unmappedAtoms1() const
{
    if (this->unmapped_atms1.isEmpty())
    {
        return SelectorM<Atom>();
    }
    else
    {
        return this->orig_atms1[this->unmapped_atms1];
    }
}

/** Return whether or not the passed reference atom has been
 *  mapped to a mapped atom
 */
bool AtomMapping::isMapped(const Atom &atom) const
{
    return this->atms0.contains(atom);
}

/** Return whether or not this mapping refers to only a single molecule */
bool AtomMapping::isSingleMolecule() const
{
    return this->orig_atms0.isSingleMolecule() and this->orig_atms1.isSingleMolecule();
}

/** Assert that this mapping refers only to a single molecule */
void AtomMapping::assertSingleMolecule() const
{
    if (not this->isSingleMolecule())
    {
        throw SireError::incompatible_error(QObject::tr(
                                                "This mapping refers to more than one molecule, and so "
                                                "cannot be used in this context."),
                                            CODELOC);
    }
}

/** Return the property map used to find properties of the
 *  reference molecule
 */
const PropertyMap &AtomMapping::propertyMap0() const
{
    return this->map0;
}

/** Return the property map used to find properties of the
 *  mapped molecule
 */
const PropertyMap &AtomMapping::propertyMap1() const
{
    return this->map1;
}

/** Return the mapping where the perturbed state (1) has been
 *  aligned against the reference state (0).
 */
AtomMapping AtomMapping::align() const
{
    return this->alignTo0();
}

/** Return the mapping where the perturbed state (1) has been
 *  aligned against the reference state (0).
 */
AtomMapping AtomMapping::alignTo0() const
{
    if (this->isEmpty())
        return *this;
    else if (this->atms0.isEmpty())
        return *this;

    QVector<Vector> coords0 = this->atms0.property<Vector>(map0["coordinates"]).toVector();
    QVector<Vector> coords1 = this->atms1.property<Vector>(map1["coordinates"]).toVector();

    AtomMapping ret(*this);

    if ((this->count() == 1) and ((coords0.count() == 1) or (coords1.count() == 1)))
    {
        // if we've only mapped a single atom and one molecule is a monatomic ion,
        // then simply replace the coordinates of the mapped atom.

        auto atom0 = this->atms0[0];
        auto atom1 = this->atms1[0];

        auto mol = this->orig_atms1.molecules()[0];

        mol = mol.edit().atom(atom1.index())
                 .setProperty(map1["coordinates"].source(), coords0[0])
                 .molecule().commit();

        ret.atms1.update(mol);
        ret.orig_atms1.update(mol);
    }
    else
    {
        // calculate the transform to do a RMSD aligment of the two sets of coordinates
        auto transform = SireMaths::getAlignment(coords0, coords1, true);

        auto mols1 = this->orig_atms1.molecules();

        for (int i = 0; i < mols1.count(); ++i)
        {
            auto mol = mols1[i].move().transform(transform, map1).commit();

            ret.atms1.update(mol);
            ret.orig_atms1.update(mol);
        }
    }

    return ret;
}

/** Return the mapping where the perturbed state (1) has been
 *  aligned against the reference state (0).
 */
AtomMapping AtomMapping::alignTo1() const
{
    if (this->isEmpty())
        return *this;
    else if (this->atms0.isEmpty())
        return *this;

    QVector<Vector> coords0 = this->atms0.property<Vector>(map0["coordinates"]).toVector();
    QVector<Vector> coords1 = this->atms1.property<Vector>(map1["coordinates"]).toVector();

    AtomMapping ret(*this);

    if ((this->count() == 1) and ((coords0.count() == 1) or (coords1.count() == 1)))
    {
        // if we've only mapped a single atom and one molecule is a monatomic ion,
        // then simply replace the coordinates of the mapped atom.

        auto atom0 = this->atms0[0];
        auto atom1 = this->atms1[0];

        auto mol = this->orig_atms0.molecules()[0];

        mol = mol.edit().atom(atom0.index())
                 .setProperty(map0["coordinates"].source(), coords0[0])
                 .molecule().commit();

        ret.atms0.update(mol);
        ret.orig_atms0.update(mol);
    }
    else
    {
        // calculate the transform to do a RMSD aligment of the two sets of coordinates
        auto transform = SireMaths::getAlignment(coords1, coords0, true);

        auto mols0 = this->orig_atms0.molecules();

        for (int i = 0; i < mols0.count(); ++i)
        {
            auto mol = mols0[i].move().transform(transform, map0).commit();

            ret.atms0.update(mol);
            ret.orig_atms0.update(mol);
        }
    }

    return ret;
}
