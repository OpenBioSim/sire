/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "sireglobal.h"

#include "twoatomfunctions.h"

#include "SireBase/console.h"

#include "SireCAS/symbols.h"

#include "SireMol/atommatcher.h"
#include "SireMol/atomselection.h"
#include "SireMol/moleculeinfodata.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMM::detail;

using namespace SireBase;
using namespace SireCAS;
using namespace SireMol;
using namespace SireStream;

//////
////// Implementation of TwoAtomFunction
//////

QDataStream &operator<<(QDataStream &ds, const TwoAtomFunction &twoatomfunc)
{
    ds << twoatomfunc.atm0 << twoatomfunc.atm1 << static_cast<const AtomFunction &>(twoatomfunc);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, TwoAtomFunction &twoatomfunc)
{
    ds >> twoatomfunc.atm0 >> twoatomfunc.atm1 >> static_cast<AtomFunction &>(twoatomfunc);

    return ds;
}

/** Constructor */
TwoAtomFunction::TwoAtomFunction() : AtomFunction()
{
}

/** Construct for the specified pair of atoms with the specified function */
TwoAtomFunction::TwoAtomFunction(const CGAtomIdx &atom0, const CGAtomIdx &atom1, const SireCAS::Expression &function)
    : AtomFunction(function), atm0(atom0), atm1(atom1)
{
}

/** Copy constructor */
TwoAtomFunction::TwoAtomFunction(const TwoAtomFunction &other) : AtomFunction(other), atm0(other.atm0), atm1(other.atm1)
{
}

/** Destructor */
TwoAtomFunction::~TwoAtomFunction()
{
}

/** Copy assignment operator */
TwoAtomFunction &TwoAtomFunction::operator=(const TwoAtomFunction &other)
{
    AtomFunction::operator=(other);
    atm0 = other.atm0;
    atm1 = other.atm1;

    return *this;
}

/** Comparison operator */
bool TwoAtomFunction::operator==(const TwoAtomFunction &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and AtomFunction::operator==(other);
}

/** Comparison operator */
bool TwoAtomFunction::operator!=(const TwoAtomFunction &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1 or AtomFunction::operator!=(other);
}

/** Return a string representation */
QString TwoAtomFunction::toString() const
{
    return QObject::tr("TwoAtomFunction( %1 <-> %2 : %3 )")
        .arg(atm0.toString(), atm1.toString(), this->function().toString());
}

//////
////// Implementation of detail::IDPair
//////

SIRE_ALWAYS_INLINE QDataStream &operator<<(QDataStream &ds, const IDPair &idpair)
{
    ds << idpair.atom0 << idpair.atom1;
    return ds;
}

SIRE_ALWAYS_INLINE QDataStream &operator>>(QDataStream &ds, IDPair &idpair)
{
    ds >> idpair.atom0 >> idpair.atom1;
    return ds;
}

IDPair::IDPair(quint32 atm0, quint32 atm1) : atom0(atm0), atom1(atm1)
{
    if (atm0 > atm1)
    {
        qSwap(atom0, atom1);
    }
}

IDPair::IDPair(const IDPair &other) : atom0(other.atom0), atom1(other.atom1)
{
}

IDPair::~IDPair()
{
}

IDPair &IDPair::operator=(const IDPair &other)
{
    atom0 = other.atom0;
    atom1 = other.atom1;
    return *this;
}

bool IDPair::operator==(const IDPair &other) const
{
    return atom0 == other.atom0 and atom1 == other.atom1;
}

bool IDPair::operator!=(const IDPair &other) const
{
    return atom0 != other.atom0 or atom1 != other.atom1;
}

bool IDPair::operator<(const IDPair &other) const
{
    return atom0 < other.atom0 or (atom0 == other.atom0 and atom1 < other.atom1);
}

bool IDPair::operator<=(const IDPair &other) const
{
    return atom0 < other.atom0 or (atom0 == other.atom0 and atom1 <= other.atom1);
}

bool IDPair::operator>(const IDPair &other) const
{
    return atom0 > other.atom0 or (atom0 == other.atom0 and atom1 > other.atom1);
}

bool IDPair::operator>=(const IDPair &other) const
{
    return atom0 > other.atom0 or (atom0 == other.atom0 and atom1 >= other.atom1);
}

//////
////// Implementation of TwoAtomFunctions
//////

static const RegisterMetaType<TwoAtomFunctions> r_twoatomfuncs;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const TwoAtomFunctions &twoatomfuncs)
{
    writeHeader(ds, r_twoatomfuncs, 1);

    SharedDataStream sds(ds);

    sds << twoatomfuncs.potentials_by_atoms << static_cast<const AtomFunctions &>(twoatomfuncs);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, TwoAtomFunctions &twoatomfuncs)
{
    VersionID v = readHeader(ds, r_twoatomfuncs);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> twoatomfuncs.potentials_by_atoms >> static_cast<AtomFunctions &>(twoatomfuncs);
    }
    else
        throw version_error(v, "1", r_twoatomfuncs, CODELOC);

    return ds;
}

/** Constructor */
TwoAtomFunctions::TwoAtomFunctions() : ConcreteProperty<TwoAtomFunctions, AtomFunctions>()
{
}

/** Construct the container to hold the set of two-atom functions
    for the molecule whose data is in 'moldata' */
TwoAtomFunctions::TwoAtomFunctions(const MoleculeData &moldata)
    : ConcreteProperty<TwoAtomFunctions, AtomFunctions>(moldata)
{
}

/** Construct the container to hold the set of two-atom functions
    for the molecule whose layout is in 'molinfo' */
TwoAtomFunctions::TwoAtomFunctions(const MoleculeInfoData &molinfo)
    : ConcreteProperty<TwoAtomFunctions, AtomFunctions>(molinfo)
{
}

/** Copy constructor */
TwoAtomFunctions::TwoAtomFunctions(const TwoAtomFunctions &other)
    : ConcreteProperty<TwoAtomFunctions, AtomFunctions>(other), potentials_by_atoms(other.potentials_by_atoms)
{
}

/** Destructor */
TwoAtomFunctions::~TwoAtomFunctions()
{
}

/** Copy assignment operator */
TwoAtomFunctions &TwoAtomFunctions::operator=(const TwoAtomFunctions &other)
{
    AtomFunctions::operator=(other);
    potentials_by_atoms = other.potentials_by_atoms;

    return *this;
}

/** Comparison operator */
bool TwoAtomFunctions::operator==(const TwoAtomFunctions &other) const
{
    return AtomFunctions::operator==(other) and potentials_by_atoms == other.potentials_by_atoms;
}

/** Comparison operator */
bool TwoAtomFunctions::operator!=(const TwoAtomFunctions &other) const
{
    return AtomFunctions::operator!=(other) or potentials_by_atoms != other.potentials_by_atoms;
}

inline QString _id_string(const MoleculeInfoData &info, int atom)
{
    return QString("%1:%2").arg(info.name(AtomIdx(atom))).arg(info.number(AtomIdx(atom)));
}

inline QString _pretty_string(const MoleculeInfoData &info, const IDPair &pair, const Expression &func)
{
    QString id = QString("%1-%2").arg(_id_string(info, pair.atom0), 7).arg(_id_string(info, pair.atom1), -7);

    return QString("%1 : %2").arg(id, -15).arg(func.toString());
}

/** Return a string representation */
QString TwoAtomFunctions::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("TwoAtomFunctions::empty");
    }
    else
    {
        QStringList parts;

        auto keys = potentials_by_atoms.keys();
        const int n = keys.count();
        std::sort(keys.begin(), keys.end());

        if (n <= 10)
        {
            for (int i = 0; i < n; ++i)
            {
                parts.append(
                    QObject::tr("%1: %2").arg(i).arg(_pretty_string(info(), keys[i], potentials_by_atoms[keys[i]])));
            }
        }
        else
        {
            for (int i = 0; i < 5; ++i)
            {
                parts.append(
                    QObject::tr("%1: %2").arg(i).arg(_pretty_string(info(), keys[i], potentials_by_atoms[keys[i]])));
            }

            parts.append("...");

            for (int i = n - 5; i < n; ++i)
            {
                parts.append(
                    QObject::tr("%1: %2").arg(i).arg(_pretty_string(info(), keys[i], potentials_by_atoms[keys[i]])));
            }
        }

        return QObject::tr("TwoAtomFunctions( size=%1\n%2\n)").arg(n).arg(parts.join("\n"));
    }
}

/** Set the potential energy function used by atoms 'atom0' and 'atom1'
    to be equal to 'expression' - this replaces any existing expression

    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
void TwoAtomFunctions::set(AtomIdx atom0, AtomIdx atom1, const Expression &expression)
{
    quint32 atm0 = atom0.map(info().nAtoms());
    quint32 atm1 = atom1.map(info().nAtoms());

    if (atm0 == atm1)
        throw SireMol::duplicate_atom(
            QObject::tr("You cannot add a function that acts between the same atom! (%1)").arg(atm0), CODELOC);

    potentials_by_atoms.insert(IDPair(atm0, atm1), expression);
    AtomFunctions::addSymbols(expression.symbols());
}

/** Set the potential energy function used by atoms 'atom0' and 'atom1'
    to be equal to 'expression' - this replaces any existing expression

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void TwoAtomFunctions::set(const AtomID &atom0, const AtomID &atom1, const Expression &expression)
{
    this->set(info().atomIdx(atom0), info().atomIdx(atom1), expression);
}

/** Set the potential energy function used for the bond identified by 'bondid'
    to be equal to 'expression' - this replaces any existing expression

    Note that this replaces both 1-2 and 2-1

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void TwoAtomFunctions::set(const BondID &bondid, const Expression &expression)
{
    AtomIdx atom0 = info().atomIdx(bondid.atom0());
    AtomIdx atom1 = info().atomIdx(bondid.atom1());

    this->clear(atom0, atom1);
    this->clear(atom1, atom0);

    this->set(atom0, atom1, expression);
}

/** Check if any of the symbols in 'symbols' need to be removed... */
void TwoAtomFunctions::removeSymbols(QSet<Symbol> symbols)
{
    for (QHash<IDPair, Expression>::const_iterator it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd(); ++it)
    {
        if (symbols.isEmpty())
            return;

        symbols.subtract(it.value().symbols());
    }

    // the only remaining symbols are ones that no longer exist
    // in this set
    AtomFunctions::removeSymbols(symbols);
}

/** Clear any function that acts between the atoms 'atom0' and 'atom1'

    \throw SireError::invalid_index
*/
void TwoAtomFunctions::clear(AtomIdx atom0, AtomIdx atom1)
{
    quint32 atm0 = atom0.map(info().nAtoms());
    quint32 atm1 = atom1.map(info().nAtoms());

    TwoAtomFunctions::removeSymbols(potentials_by_atoms.take(IDPair(atm0, atm1)).symbols());
}

/** Clear all functions that involve the atom 'atom'

    \throw SireError::invalid_index
*/
void TwoAtomFunctions::clear(AtomIdx atom)
{
    quint32 atm = atom.map(info().nAtoms());

    QList<IDPair> keys = potentials_by_atoms.keys();

    foreach (const IDPair &key, keys)
    {
        if (key.atom0 == atm or key.atom1 == atm)
        {
            TwoAtomFunctions::removeSymbols(potentials_by_atoms.take(key).symbols());
        }
    }
}

/** Clear all functions that invole any of the atoms in 'atoms'
 *  - if 'exclusive' is true, then this only removes functions
 *  that exclusively involve these atoms - if false, then
 *  if removes functions that involve any of these atoms
 */
void TwoAtomFunctions::clear(const QList<AtomIdx> &atoms, bool exclusive)
{
    QSet<quint32> atms;
    atms.reserve(atoms.count());

    for (const auto &atom : atoms)
    {
        atms.insert(atom.map(info().nAtoms()));
    }

    QList<IDPair> keys = potentials_by_atoms.keys();

    if (exclusive)
    {
        for (const auto &key : keys)
        {
            if (atms.contains(key.atom0) and atms.contains(key.atom1))
            {
                TwoAtomFunctions::removeSymbols(potentials_by_atoms.take(key).symbols());
            }
        }
    }
    else
    {
        for (const auto &key : keys)
        {
            if (atms.contains(key.atom0) or atms.contains(key.atom1))
            {
                TwoAtomFunctions::removeSymbols(potentials_by_atoms.take(key).symbols());
            }
        }
    }
}

/** Clear any function that acts on the atoms identified by 'atom'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void TwoAtomFunctions::clear(const AtomID &atom)
{
    QList<AtomIdx> atomidxs = atom.map(info());

    foreach (AtomIdx atomidx, atomidxs)
    {
        this->clear(atomidx);
    }
}

/** Clear any function that acts between the atoms 'atom0' and 'atom1'

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void TwoAtomFunctions::clear(const AtomID &atom0, const AtomID &atom1)
{
    QList<AtomIdx> atoms0 = atom0.map(info());
    QList<AtomIdx> atoms1 = atom1.map(info());

    foreach (AtomIdx atm0, atoms0)
    {
        foreach (AtomIdx atm1, atoms1)
        {
            this->clear(atm0, atm1);
        }
    }
}

/** Clear the potential that acts over the bond identified by 'bondid'

    Note that this removes both 1-2 and 2-1

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void TwoAtomFunctions::clear(const BondID &bondid)
{
    this->clear(bondid.atom0(), bondid.atom1());
    this->clear(bondid.atom1(), bondid.atom0());
}

/** Completely clear all of the functions from this set */
void TwoAtomFunctions::clear()
{
    potentials_by_atoms.clear();
    AtomFunctions::removeSymbols();
}

/** Perform the substitutions contained in 'identities' in all of
    the expressions in this set. This could be useful if you have
    defined these expressions with respect to a lambda parameter,
    and now want to set that value of lambda */
void TwoAtomFunctions::substitute(const Identities &identities)
{
    AtomFunctions::removeSymbols();

    for (QHash<IDPair, Expression>::iterator it = potentials_by_atoms.begin(); it != potentials_by_atoms.end(); ++it)
    {
        it.value() = it.value().substitute(identities);
        AtomFunctions::addSymbols(it.value().symbols());
    }
}

/** Return whether or not this is empty (has no potentials for any internals) */
bool TwoAtomFunctions::isEmpty() const
{
    return potentials_by_atoms.isEmpty();
}

/** Return the function acting between the atoms 'atom0' and 'atom1'.
    This returns an empty expression if there is no expression between
    these atoms

    \throw SireError::invalid_index
    \throw SireMol::duplicate_atom
*/
Expression TwoAtomFunctions::potential(AtomIdx atom0, AtomIdx atom1) const
{
    quint32 atm0 = atom0.map(info().nAtoms());
    quint32 atm1 = atom1.map(info().nAtoms());

    if (atm0 == atm1)
        throw SireMol::duplicate_atom(
            QObject::tr("There is no potential that acts between the same atom! (%1)").arg(atm0), CODELOC);

    return potentials_by_atoms.value(IDPair(atm0, atm1));
}

/** Return the function acting between the atoms 'atom0' and 'atom1'.
    This returns an empty expression if there is no expression between
    these atoms

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression TwoAtomFunctions::potential(const AtomID &atom0, const AtomID &atom1) const
{
    return this->potential(info().atomIdx(atom0), info().atomIdx(atom1));
}

/** Return the function acting on the bond identified by 'bondid'.
    This returns an empty expression if there is no expression on
    this bond

    This searches first for the function 1-2, and if that is not
    found then it returns the function for 2-1

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression TwoAtomFunctions::potential(const BondID &bondid) const
{
    AtomIdx atom0 = info().atomIdx(bondid.atom0());
    AtomIdx atom1 = info().atomIdx(bondid.atom1());

    quint32 atm0 = atom0.map(info().nAtoms());
    quint32 atm1 = atom1.map(info().nAtoms());

    if (potentials_by_atoms.contains(IDPair(atm0, atm1)))
    {
        return potentials_by_atoms.value(IDPair(atm0, atm1));
    }
    else
        return potentials_by_atoms.value(IDPair(atm1, atm0));
}

/** Return the force (derivative of the potential with respect to 'symbol')
    between the atoms 'atom0' and 'atom1'

    \throw SireError::invalid_index
*/
Expression TwoAtomFunctions::force(AtomIdx atom0, AtomIdx atom1, const Symbol &symbol) const
{
    return -(this->potential(atom0, atom1).differentiate(symbol));
}

/** Return the force (derivative of the potential with respect to 'symbol')
    between the atoms 'atom0' and 'atom1'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression TwoAtomFunctions::force(const AtomID &atom0, const AtomID &atom1, const Symbol &symbol) const
{
    return -(this->potential(atom0, atom1).differentiate(symbol));
}

/** Return the force (derivative of the potential with respect to 'symbol')
    on the bond identified by 'bondid'

    This searches first for the function 1-2, and if that is not
    found then it returns the function for 2-1

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
Expression TwoAtomFunctions::force(const BondID &bondid, const Symbol &symbol) const
{
    return -(this->potential(bondid).differentiate(symbol));
}

/** Return the potential energy functions acting between the identified
    pairs of atoms - if exclusive is true then only return potentials where
    both atoms are in the bond
*/
QVector<TwoAtomFunction> TwoAtomFunctions::potentials(const QList<AtomIdx> &atms, bool exclusive) const
{
    QVector<TwoAtomFunction> funcs;
    funcs.reserve(potentials_by_atoms.count());

    QSet<AtomIdx> atoms(atms.begin(), atms.end());

    for (QHash<IDPair, Expression>::const_iterator it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd(); ++it)
    {
        if (exclusive)
        {
            if (atoms.contains(AtomIdx(it.key().atom0)) and atoms.contains(AtomIdx(it.key().atom1)))
            {
                funcs.append(TwoAtomFunction(info().cgAtomIdx(AtomIdx(it.key().atom0)),
                                             info().cgAtomIdx(AtomIdx(it.key().atom1)), it.value()));
            }
        }
        else
        {
            if (atoms.contains(AtomIdx(it.key().atom0)) or atoms.contains(AtomIdx(it.key().atom1)))
            {
                funcs.append(TwoAtomFunction(info().cgAtomIdx(AtomIdx(it.key().atom0)),
                                             info().cgAtomIdx(AtomIdx(it.key().atom1)), it.value()));
            }
        }
    }

    return funcs;
}

/** Return the potential energy functions acting between the identified
    pairs of atoms */
QVector<TwoAtomFunction> TwoAtomFunctions::potentials() const
{
    QVector<TwoAtomFunction> funcs(potentials_by_atoms.count());

    TwoAtomFunction *funcs_array = funcs.data();

    int i = 0;

    for (QHash<IDPair, Expression>::const_iterator it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd(); ++it)
    {
        funcs_array[i] = TwoAtomFunction(info().cgAtomIdx(AtomIdx(it.key().atom0)),
                                         info().cgAtomIdx(AtomIdx(it.key().atom1)), it.value());

        ++i;
    }

    return funcs;
}

/** Return the force functions acting between the identified
    pairs of atoms, for the given symbol */
QVector<TwoAtomFunction> TwoAtomFunctions::forces(const Symbol &symbol) const
{
    QVector<TwoAtomFunction> forces;
    forces.reserve(potentials_by_atoms.count());

    for (QHash<IDPair, Expression>::const_iterator it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd(); ++it)
    {
        Expression force = it.value().differentiate(symbol);

        if (not force.isZero())
        {
            forces.append(TwoAtomFunction(info().cgAtomIdx(AtomIdx(it.key().atom0)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom1)), -force));
        }
    }

    return forces;
}

/** Return the set of functions where only functions that involve the
    atoms in 'selected_atoms' are included. If 'isstrict' is true, then
    only include functions where all of the atoms are in 'selected_atoms',
    while if 'isstrict' is false, include functions where at least one
    atom is in 'selected_atoms' */
TwoAtomFunctions TwoAtomFunctions::includeOnly(const AtomSelection &selected_atoms, bool isstrict) const
{
    TwoAtomFunctions ret(*this);

    QMutableHashIterator<IDPair, Expression> it(ret.potentials_by_atoms);

    if (isstrict)
    {
        while (it.hasNext())
        {
            it.next();

            if (not(selected_atoms.selected(AtomIdx(it.key().atom0)) and
                    selected_atoms.selected(AtomIdx(it.key().atom1))))
            {
                it.remove();
            }
        }
    }
    else
    {
        while (it.hasNext())
        {
            it.next();

            if (not(selected_atoms.selected(AtomIdx(it.key().atom0)) or
                    selected_atoms.selected(AtomIdx(it.key().atom1))))
            {
                it.remove();
            }
        }
    }

    return ret;
}

/** This returns the total number of functions in this set */
int TwoAtomFunctions::nFunctions() const
{
    return potentials_by_atoms.count();
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom matching
    functions in 'atommatcher' to match atoms from the current molecule
    to the atoms in the molecule whose layout is in 'molinfo'

    This only copies the TwoAtomFunction for pairs of atoms that
    are successfully matched - it does not copy functions for atoms
    that are not matched. Use TwoAtomFunctions::nFunctions() to check
    if the number of functions in the returned set is the same as
    the number in this set, if you want to ensure that all of the
    functions have been copied.

    \throw SireError::incompatible_error
*/
PropertyPtr TwoAtomFunctions::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                      const AtomMatcher &atommatcher) const
{
    if (not atommatcher.changesOrder(this->info(), molinfo))
    {
        // the order of the atoms remains the same - this means that the
        // AtomIdx indicies are still valid
        TwoAtomFunctions ret(molinfo);
        ret.potentials_by_atoms = this->potentials_by_atoms;
        return ret;
    }

    QHash<AtomIdx, AtomIdx> matched_atoms = atommatcher.match(this->info(), molinfo);

    return this->_pvt_makeCompatibleWith(molinfo, matched_atoms);
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom mapping
    in 'map' to match atoms from the current molecule to the atoms in
    the molecule whose layout is in 'molinfo'

    This only copies the TwoAtomFunction for pairs of atoms that
    are successfully matched - it does not copy functions for atoms
    that are not matched. Use TwoAtomFunctions::nFunctions() to check
    if the number of functions in the returned set is the same as
    the number in this set, if you want to ensure that all of the
    functions have been copied.

    \throw SireError::incompatible_error
*/
PropertyPtr TwoAtomFunctions::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                      const QHash<AtomIdx, AtomIdx> &map) const
{
    TwoAtomFunctions ret(molinfo);

    for (QHash<IDPair, Expression>::const_iterator it = potentials_by_atoms.constBegin();
         it != potentials_by_atoms.constEnd(); ++it)
    {
        AtomIdx new_atom0 = map.value(AtomIdx(it.key().atom0), AtomIdx(-1));
        AtomIdx new_atom1 = map.value(AtomIdx(it.key().atom1), AtomIdx(-1));

        if (new_atom0 == -1 or new_atom1 == -1)
            continue;

        ret.set(new_atom0, new_atom1, it.value());
    }

    return ret;
}

const char *TwoAtomFunctions::typeName()
{
    return QMetaType::typeName(qMetaTypeId<TwoAtomFunctions>());
}

template <class T>
QSet<T> _to_set(const QList<T> &vals)
{
    QSet<T> ret;
    ret.reserve(vals.count());

    for (const auto &val : vals)
    {
        ret.insert(val);
    }

    return ret;
}

/** Merge this property with another property */
PropertyList TwoAtomFunctions::merge(const MolViewProperty &other,
                                     const AtomIdxMapping &mapping,
                                     const QString &ghost,
                                     const SireBase::PropertyMap &map) const
{
    if (not other.isA<TwoAtomFunctions>())
    {
        throw SireError::incompatible_error(QObject::tr("Cannot merge %1 with %2 as they are different types.")
                                                .arg(this->what())
                                                .arg(other.what()),
                                            CODELOC);
    }

    if (not ghost.isEmpty())
    {
        Console::warning(QObject::tr("The ghost parameter '%1' for bond parameters is ignored").arg(ghost));
    }

    const TwoAtomFunctions &ref = *this;
    const TwoAtomFunctions &pert = other.asA<TwoAtomFunctions>();

    TwoAtomFunctions prop0 = ref;
    TwoAtomFunctions prop1 = ref;

    // the prop1 properties are made by finding all of the atoms that
    // are involved in bonds in 'pert' and removing any bonds involving
    // only those atoms from 'prop1', and then adding back the matching
    // bonds from 'pert'. Use 'true' to only remove bonds where both
    // atoms are in the mapping
    prop1.clear(mapping.mappedIn1(), true);

    // get the mapping from the perturbed to reference states, including
    // atoms that don't exist in the reference state. In all cases,
    // the values are the indexes in the merged molecule
    auto map1to0 = mapping.map1to0(true);

    // now find all of the bonds in 'pert' where both atoms in the
    // bond are in map1to0.keys() - i.e. exist and are mapped from
    // the perturbed state
    const auto pert_bonds = pert.potentials(map1to0.keys(), true);

    for (const auto &pert_bond : pert_bonds)
    {
        const auto atom0 = map1to0.value(info().atomIdx(pert_bond.atom0()));
        const auto atom1 = map1to0.value(info().atomIdx(pert_bond.atom1()));

        prop1.set(atom0, atom1, pert_bond.function());

        if (mapping.isUnmappedIn0(atom0) or mapping.isUnmappedIn0(atom1))
        {
            // the prop0 properties are nearly correct - we just need to add
            // in bonds from 'pert' that involve the atoms that are not mapped
            // in the reference state - this way, those added atoms are held
            // by a constant bond potential, so won't fly away in the
            // simulation of the reference state
            prop0.set(atom0, atom1, pert_bond.function());
        }
    }

    // now add in the bonds to the perturbed state from the reference
    // state for any atoms that aren't mapped to the perturbed state.
    // This way, the removed atoms are held by a constant bond potential,
    // so won't fly away in the simulation of the perturbed state
    auto map0to1 = mapping.map0to1(true);

    const auto ref_bonds = prop0.potentials(map0to1.keys(), true);

    for (const auto &ref_bond : ref_bonds)
    {
        const auto atom0 = info().atomIdx(ref_bond.atom0());
        const auto atom1 = info().atomIdx(ref_bond.atom1());

        if (mapping.isUnmappedIn1(atom0) or mapping.isUnmappedIn1(atom1))
        {
            prop1.set(atom0, atom1, ref_bond.function());
        }
    }

    // check if we are allowed to change the size of a ring or break rings
    bool allow_ring_breaking = true;
    bool allow_ring_size_change = true;

    if (map.specified("allow_ring_breaking"))
    {
        allow_ring_breaking = map["allow_ring_breaking"].value().asABoolean();
    }

    if (map.specified("allow_ring_size_change"))
    {
        allow_ring_size_change = map["allow_ring_size_change"].value().asABoolean();
    }

    if (not(allow_ring_breaking or allow_ring_size_change))
    {
        if (prop0.nFunctions() != prop1.nFunctions())
        {
            // number of bond functions has changed - this indicates
            // (but not necessarily proves) that a ring has been broken
            // or a ring size has changed
            throw SireError::incompatible_error(
                QObject::tr("The number of bonds in the reference (%1) and "
                            "perturbed (%2) states is different, indicating that a ring has been broken or a ring size has changed. "
                            "If you want to allow this perturbation, set 'allow_ring_breaking' or 'allow_ring_size_change' to true.")
                    .arg(prop0.nFunctions())
                    .arg(prop1.nFunctions()),
                CODELOC);
        }
    }

    SireBase::PropertyList ret;

    ret.append(prop0);
    ret.append(prop1);

    return ret;
}
