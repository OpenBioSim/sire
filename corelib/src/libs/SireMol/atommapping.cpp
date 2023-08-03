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

#include "atommapping.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<AtomMapping> r_mapping;

QDataStream &operator<<(QDataStream &ds, const AtomMapping &mapping)
{
    writeHeader(ds, r_mapping, 1);

    SharedDataStream sds(ds);

    sds << mapping.atms0 << mapping.atms1
        << static_cast<const Property &>(mapping);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AtomMapping &mapping)
{
    VersionID v = readHeader(ds, r_mapping);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mapping.atms0 >> mapping.atms1 >> static_cast<Property &>(mapping);
    }
    else
        throw version_error(v, "1", r_mapping, CODELOC);

    return ds;
}

AtomMapping::AtomMapping() : ConcreteProperty<AtomMapping, Property>()
{
}

AtomMapping::AtomMapping(const SelectorM<Atom> &atoms0,
                         const SelectorM<Atom> &atoms1)
    : ConcreteProperty<AtomMapping, Property>(),
      atms0(atoms0), atms1(atoms1)
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

AtomMapping::AtomMapping(const AtomMapping &other)
    : ConcreteProperty<AtomMapping, Property>(other),
      atms0(other.atms0), atms1(other.atms1)
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
        Property::operator=(other);
    }

    return *this;
}

bool AtomMapping::operator==(const AtomMapping &other) const
{
    return atms0 == other.atms0 and atms1 == other.atms1;
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

        return QObject::tr("AtomMapping( size=%1\n%2\n)").arg(n).arg(parts.join("\n"));
    }
}

/** Return the reference atoms. We map from these atom to
 *  the mapped atoms (atoms1) */
const SelectorM<Atom> &AtomMapping::atoms0() const
{
    return this->atms0;
}

/** Return the mapped atoms. We map from the reference atoms (atoms0)
 *  to these atoms. */
const SelectorM<Atom> &AtomMapping::atoms1() const
{
    return this->atms1;
}

/** Return an AtomMapping that swaps the reference and mapped atoms */
AtomMapping AtomMapping::swap() const
{
    return AtomMapping(this->atms1, this->atms0);
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
SelectorM<Atom> AtomMapping::map(const Selector<Atom> &atoms,
                                 bool find_all) const
{
    return this->map(SelectorM<Atom>(atoms), find_all);
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
Selector<Atom> AtomMapping::find(const Selector<Atom> &atoms,
                                 const MoleculeView &container,
                                 bool find_all) const
{
    auto found = this->find(SelectorM<Atom>(atoms),
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
SelectorM<Atom> AtomMapping::find(const Selector<Atom> &atoms,
                                  const SelectorM<Atom> &container,
                                  bool find_all) const
{
    return this->find(SelectorM<Atom>(atoms), container, find_all);
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
