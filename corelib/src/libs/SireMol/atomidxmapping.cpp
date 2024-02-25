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

#include "atomidxmapping.h"

#include "moleculeinfodata.h"

#include "SireID/index.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

////////
//////// Implementation of AtomIdxMappingEntry
////////

RegisterMetaType<AtomIdxMappingEntry> r_entry(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const AtomIdxMappingEntry &entry)
{
    writeHeader(ds, r_entry, 1);

    SharedDataStream sds(ds);

    sds << entry.atomidx0 << entry.atomidx1
        << entry.cgatomidx0 << entry.cgatomidx1
        << entry.unmapped0;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AtomIdxMappingEntry &entry)
{
    auto v = readHeader(ds, r_entry);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> entry.atomidx0 >> entry.atomidx1 >> entry.cgatomidx0 >> entry.cgatomidx1 >> entry.unmapped0;
    }
    else
        throw SireStream::version_error(v, "1", r_entry, CODELOC);

    return ds;
}

/** Null constructor */
AtomIdxMappingEntry::AtomIdxMappingEntry()
    : atomidx0(), atomidx1(), cgatomidx0(), cgatomidx1(), unmapped0(false)
{
}

/** Construct to map from 'index0' in 'molinfo0' to 'index1' in 'molinfo1'.
 *  There is no mapping in the perturbed state if 'index1' is null. There
 *  is no mapping in the reference state if 'is_unmapped_in_reference' is true.
 *  It is an error to have 'index0' null, or if you cannot look up
 *  'index0' in 'molinfo0' or 'index1' in 'molinfo1'.
 */
AtomIdxMappingEntry::AtomIdxMappingEntry(const AtomIdx &index0, const AtomIdx &index1,
                                         const MoleculeInfoData &molinfo0,
                                         const MoleculeInfoData &molinfo1,
                                         bool is_unmapped_in_reference)
    : atomidx0(index0), atomidx1(index1), unmapped0(is_unmapped_in_reference)
{
    if (atomidx0.isNull())
    {
        throw SireError::incompatible_error(QObject::tr("The atom index in the reference state is null!"), CODELOC);
    }

    cgatomidx0 = molinfo0.cgAtomIdx(atomidx0);

    if (not atomidx1.isNull())
    {
        cgatomidx1 = molinfo1.cgAtomIdx(atomidx1);
    }
}

/** Copy constructor */
AtomIdxMappingEntry::AtomIdxMappingEntry(const AtomIdxMappingEntry &other)
    : atomidx0(other.atomidx0), atomidx1(other.atomidx1),
      cgatomidx0(other.cgatomidx0), cgatomidx1(other.cgatomidx1),
      unmapped0(other.unmapped0)
{
}

/** Destructor */
AtomIdxMappingEntry::~AtomIdxMappingEntry()
{
}

/** Assignment operator */
AtomIdxMappingEntry &AtomIdxMappingEntry::operator=(const AtomIdxMappingEntry &other)
{
    if (this != &other)
    {
        atomidx0 = other.atomidx0;
        atomidx1 = other.atomidx1;
        cgatomidx0 = other.cgatomidx0;
        cgatomidx1 = other.cgatomidx1;
        unmapped0 = other.unmapped0;
    }

    return *this;
}

/** Equality operator */
bool AtomIdxMappingEntry::operator==(const AtomIdxMappingEntry &other) const
{
    return atomidx0 == other.atomidx0 and atomidx1 == other.atomidx1 and
           cgatomidx0 == other.cgatomidx0 and cgatomidx1 == other.cgatomidx1 and
           unmapped0 == other.unmapped0;
}

/** Inequality operator */
bool AtomIdxMappingEntry::operator!=(const AtomIdxMappingEntry &other) const
{
    return not this->operator==(other);
}

/** Clone this object */
AtomIdxMappingEntry *AtomIdxMappingEntry::clone() const
{
    return new AtomIdxMappingEntry(*this);
}

/** Return whether or not this is a null entry */
bool AtomIdxMappingEntry::isNull() const
{
    return atomidx0.isNull();
}

/** Convert this object to a string */
QString AtomIdxMappingEntry::toString() const
{
    if (this->isNull())
    {
        return QString("AtomIdxMappingEntry: null");
    }
    else if (this->isUnmappedIn0())
    {
        return QString("%1 unmapped -> %2").arg(this->atomIdx0().value()).arg(this->atomIdx1().value());
    }
    else if (this->isUnmappedIn1())
    {
        return QString("%1 -> unmapped").arg(this->atomIdx0().value());
    }
    else
    {
        return QString("%1 -> %2").arg(this->atomIdx0().value()).arg(this->atomIdx1().value());
    }
}

const char *AtomIdxMappingEntry::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AtomIdxMappingEntry>());
}

const char *AtomIdxMappingEntry::what() const
{
    return AtomIdxMappingEntry::typeName();
}

/** Return whether or not this atom is unmapped in the reference state */
bool AtomIdxMappingEntry::isUnmappedIn0() const
{
    return unmapped0;
}

/** Return whether or not this atom is unmapped in the perturbed state */
bool AtomIdxMappingEntry::isUnmappedIn1() const
{
    return atomidx1.isNull();
}

/** Return the atom index in the reference state. This will always have
 *  a value, even if the atom is unmapped in the reference state
 *  (this signals that any parameter with this index should be zero)
 */
AtomIdx AtomIdxMappingEntry::atomIdx0() const
{
    return atomidx0;
}

/** Return the atom index in the perturbed state, or a null index if
 *  the atom is unmapped in the perturbed state */
AtomIdx AtomIdxMappingEntry::atomIdx1() const
{
    return atomidx1;
}

/** Return the atom index in the reference state. This will always have
 *  a value, even if the atom is unmapped in the reference state
 *  (this signals that any parameter with this index should be zero)
 */
CGAtomIdx AtomIdxMappingEntry::cgAtomIdx0() const
{
    return cgatomidx0;
}

/** Return the atom index in the perturbed state, or a null index if
 *  the atom is unmapped in the perturbed state */
CGAtomIdx AtomIdxMappingEntry::cgAtomIdx1() const
{
    return cgatomidx1;
}

////////
//////// Implementation of AtomIdxMapping
////////

RegisterMetaType<AtomIdxMapping> r_mapping;

QDataStream &operator<<(QDataStream &ds, const AtomIdxMapping &mapping)
{
    writeHeader(ds, r_mapping, 1);

    SharedDataStream sds(ds);

    sds << mapping.entries << static_cast<const Property &>(mapping);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AtomIdxMapping &mapping)
{
    auto v = readHeader(ds, r_mapping);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> mapping.entries >> static_cast<Property &>(mapping);
    }
    else
        throw SireStream::version_error(v, "1", r_mapping, CODELOC);

    return ds;
}

/** Null constructor */
AtomIdxMapping::AtomIdxMapping()
    : ConcreteProperty<AtomIdxMapping, Property>()
{
}

/** Construct from a single entry */
AtomIdxMapping::AtomIdxMapping(const AtomIdxMappingEntry &entry)
    : ConcreteProperty<AtomIdxMapping, Property>()
{
    if (not entry.isNull())
        entries.append(entry);
}

/** Assert that this object is sane */
void AtomIdxMapping::assertSane() const
{
    QHash<qint64, AtomIdxMappingEntry> seen;
    seen.reserve(entries.size());

    for (const auto &entry : entries)
    {
        if (entry.isNull())
        {
            throw SireError::incompatible_error(QObject::tr("The AtomIdxMapping contains a null entry!"), CODELOC);
        }
        else if (seen.contains(entry.atomIdx0().value()))
        {
            throw SireError::incompatible_error(QObject::tr("The AtomIdxMapping contains a duplicate entry (%1 vs %2)!")
                                                    .arg(entry.toString())
                                                    .arg(seen[entry.atomIdx0().value()].toString()),
                                                CODELOC);
        }
        else
        {
            seen.insert(entry.atomIdx0().value(), entry);
        }
    }
}

/** Construct from a list of entries */
AtomIdxMapping::AtomIdxMapping(const QList<AtomIdxMappingEntry> &e)
    : ConcreteProperty<AtomIdxMapping, Property>()
{
    bool any_null = false;

    for (const auto &entry : e)
    {
        if (entry.isNull())
        {
            any_null = true;
            break;
        }
    }

    if (not any_null)
    {
        entries = e;
    }
    else
    {
        for (const auto &entry : e)
        {
            if (not entry.isNull())
            {
                entries.append(entry);
            }
        }
    }

    this->assertSane();
}

/** Copy constructor */
AtomIdxMapping::AtomIdxMapping(const AtomIdxMapping &other)
    : ConcreteProperty<AtomIdxMapping, Property>(other)
{
}

/** Destructor */
AtomIdxMapping::~AtomIdxMapping()
{
}

/** Assignment operator */
AtomIdxMapping &AtomIdxMapping::operator=(const AtomIdxMapping &other)
{
    if (this != &other)
    {
        entries = other.entries;
        Property::operator=(other);
    }

    return *this;
}

/** Equality operator */
bool AtomIdxMapping::operator==(const AtomIdxMapping &other) const
{
    return entries == other.entries and Property::operator==(other);
}

/** Inequality operator */
bool AtomIdxMapping::operator!=(const AtomIdxMapping &other) const
{
    return not this->operator==(other);
}

/** Clone this object */
AtomIdxMapping *AtomIdxMapping::clone() const
{
    return new AtomIdxMapping(*this);
}

/** Convert this object to a string */
QString AtomIdxMapping::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("AtomIdxMapping::empty");
    }
    else
    {
        QStringList parts;

        const auto n = this->count();

        if (n <= 10)
        {
            for (int i = 0; i < n; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->operator[](i).toString()));
            }
        }
        else
        {
            for (int i = 0; i < 5; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->operator[](i).toString()));
            }

            parts.append("...");

            for (int i = n - 5; i < n; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->operator[](i).toString()));
            }
        }

        return QObject::tr("AtomIdxMapping( size=%2\n%3\n)").arg(n).arg(parts.join("\n"));
    }
}

const char *AtomIdxMapping::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AtomIdxMapping>());
}

const char *AtomIdxMapping::what() const
{
    return AtomIdxMapping::typeName();
}

/** Add another mapping to this one */
AtomIdxMapping &AtomIdxMapping::operator+=(const AtomIdxMapping &other)
{
    this->append(other);
    return *this;
}

/** Add another mapping to this one */
AtomIdxMapping &AtomIdxMapping::operator+=(const AtomIdxMappingEntry &entry)
{
    this->append(entry);
    return *this;
}

/** Add another mapping to this one */
AtomIdxMapping AtomIdxMapping::operator+(const AtomIdxMapping &other) const
{
    AtomIdxMapping result(*this);
    result += other;
    return result;
}

/** Add another mapping to this one */
AtomIdxMapping AtomIdxMapping::operator+(const AtomIdxMappingEntry &entry) const
{
    AtomIdxMapping result(*this);
    result += entry;
    return result;
}

/** Return an iterator to the beginning of the list */
QList<AtomIdxMappingEntry>::const_iterator AtomIdxMapping::begin() const
{
    return entries.begin();
}

/** Return an iterator to the beginning of the list */
QList<AtomIdxMappingEntry>::const_iterator AtomIdxMapping::constBegin() const
{
    return entries.constBegin();
}

/** Return an iterator to the end of the list */
QList<AtomIdxMappingEntry>::const_iterator AtomIdxMapping::end() const
{
    return entries.end();
}

/** Return an iterator to the end of the list */
QList<AtomIdxMappingEntry>::const_iterator AtomIdxMapping::constEnd() const
{
    return entries.constEnd();
}

/** Find an entry in the list */
QList<AtomIdxMappingEntry>::const_iterator AtomIdxMapping::find(AtomIdx atom) const
{
    return std::find_if(entries.constBegin(), entries.constEnd(), [atom](const AtomIdxMappingEntry &entry)
                        { return entry.atomIdx0() == atom; });
}

/** Find an entry in the list */
QList<AtomIdxMappingEntry>::const_iterator AtomIdxMapping::find(CGAtomIdx atom) const
{
    return std::find_if(entries.constBegin(), entries.constEnd(), [atom](const AtomIdxMappingEntry &entry)
                        { return entry.cgAtomIdx0() == atom; });
}

/** Append an entry to the list */
void AtomIdxMapping::append(const AtomIdxMappingEntry &new_entry)
{
    if (new_entry.isNull())
        return;

    for (const auto &entry : entries)
    {
        if (entry.atomIdx0() == new_entry.atomIdx0())
        {
            throw SireError::incompatible_error(QObject::tr("The AtomIdxMapping already contains an entry for atom %1!")
                                                    .arg(new_entry.atomIdx0().value()),
                                                CODELOC);
        }
    }

    entries.append(new_entry);
}

/** Append all of the passed entries of other onto this list */
void AtomIdxMapping::append(const AtomIdxMapping &other)
{
    AtomIdxMapping ret(*this);

    for (const auto &entry : other.entries)
    {
        ret.append(entry);
    }

    *this = ret;
}

/** Clear the list */
void AtomIdxMapping::clear()
{
    entries.clear();
}

/** Return whether or not the list is empty */
bool AtomIdxMapping::isEmpty() const
{
    return entries.isEmpty();
}

/** Return the size of the list */
int AtomIdxMapping::size() const
{
    return entries.size();
}

/** Return the count of the list */
int AtomIdxMapping::count() const
{
    return entries.count();
}

/** Return the entry at index 'i' */
const AtomIdxMappingEntry &AtomIdxMapping::operator[](int i) const
{
    i = SireID::Index(i).map(entries.size());
    return entries[i];
}

/** Take the entry at index 'i' */
AtomIdxMappingEntry AtomIdxMapping::take(int i)
{
    i = SireID::Index(i).map(entries.size());
    return entries.takeAt(i);
}

/** Take the entry for atom 'atom' */
AtomIdxMappingEntry AtomIdxMapping::take(const AtomIdx &atom)
{
    auto it = this->find(atom);

    if (it != this->constEnd())
    {
        return this->take(std::distance(this->constBegin(), it));
    }
    else
    {
        throw SireMol::missing_atom(QObject::tr(
                                        "The AtomIdxMapping does not contain an entry for atom %1!")
                                        .arg(atom.toString()),
                                    CODELOC);

        return AtomIdxMappingEntry();
    }
}

/** Take the entry for atom 'atom' */
AtomIdxMappingEntry AtomIdxMapping::take(const CGAtomIdx &atom)
{
    auto it = this->find(atom);

    if (it != this->constEnd())
    {
        return this->take(std::distance(this->constBegin(), it));
    }
    else
    {
        throw SireMol::missing_atom(QObject::tr(
                                        "The AtomIdxMapping does not contain an entry for atom %1!")
                                        .arg(atom.toString()),
                                    CODELOC);

        return AtomIdxMappingEntry();
    }
}

/** Remove the entry at index 'i' */
void AtomIdxMapping::remove(int i)
{
    i = SireID::Index(i).map(entries.size());
    entries.removeAt(i);
}

/** Remove the entry for atom 'atom' */
void AtomIdxMapping::remove(const AtomIdx &atom)
{
    auto it = this->find(atom);

    if (it != this->constEnd())
    {
        this->remove(std::distance(this->constBegin(), it));
    }
}

/** Remove the entry for atom 'atom' */
void AtomIdxMapping::remove(const CGAtomIdx &atom)
{
    auto it = this->find(atom);

    if (it != this->constEnd())
    {
        this->remove(std::distance(this->constBegin(), it));
    }
}
