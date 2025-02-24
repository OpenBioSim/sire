/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2025  Christopher Woods
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

#include "sireglobal.h"

#include "cmapfunctions.h"

#include "SireBase/console.h"

#include "SireMol/atommatcher.h"
#include "SireMol/atomselection.h"
#include "SireMol/moleculeinfodata.h"
#include "SireMol/moleculedata.h"
#include "SireMol/atommapping.h"

#include "SireError/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMM::detail;

using namespace SireBase;
using namespace SireMol;
using namespace SireStream;

//////
////// Implementation of CMAPFunction
//////

static const RegisterMetaType<CMAPFunction> r_cmapfunc(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CMAPFunction &func)
{
    writeHeader(ds, r_cmapfunc, 1);

    SharedDataStream sds(ds);

    sds << func.atm0 << func.atm1 << func.atm2 << func.atm3 << func.atm4
        << func.param;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CMAPFunction &func)
{
    VersionID v = readHeader(ds, r_cmapfunc);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> func.atm0 >> func.atm1 >> func.atm2 >> func.atm3 >> func.atm4 >> func.param;
    }
    else
        throw version_error(v, "1", r_cmapfunc, CODELOC);

    return ds;
}

CMAPFunction::CMAPFunction()
{
}

CMAPFunction::CMAPFunction(const CGAtomIdx &a0, const CGAtomIdx &a1, const CGAtomIdx &a2,
                           const CGAtomIdx &a3, const CGAtomIdx &a4, const CMAPParameter &p)
    : atm0(a0), atm1(a1), atm2(a2), atm3(a3), atm4(a4), param(p)
{
}

CMAPFunction::CMAPFunction(const CMAPFunction &other)
    : atm0(other.atm0), atm1(other.atm1), atm2(other.atm2), atm3(other.atm3), atm4(other.atm4), param(other.param)
{
}

CMAPFunction::~CMAPFunction()
{
}

CMAPFunction &CMAPFunction::operator=(const CMAPFunction &other)
{
    atm0 = other.atm0;
    atm1 = other.atm1;
    atm2 = other.atm2;
    atm3 = other.atm3;
    atm4 = other.atm4;
    param = other.param;

    return *this;
}

bool CMAPFunction::operator==(const CMAPFunction &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and atm2 == other.atm2 and atm3 == other.atm3 and atm4 == other.atm4 and param == other.param;
}

bool CMAPFunction::operator!=(const CMAPFunction &other) const
{
    return not CMAPFunction::operator==(other);
}

const char *CMAPFunction::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CMAPFunction>());
}

const char *CMAPFunction::what() const
{
    return CMAPFunction::typeName();
}

CMAPFunction *CMAPFunction::clone() const
{
    return new CMAPFunction(*this);
}

QString CMAPFunction::toString() const
{
    return QString("CMAPFunction(%1-%2-%3-%4-%5\n%6)")
        .arg(atm0.toString())
        .arg(atm1.toString())
        .arg(atm2.toString())
        .arg(atm3.toString())
        .arg(atm4.toString())
        .arg(param.toString());
}

//////
////// Implementation of detail::IDQuint
//////

QDataStream &operator<<(QDataStream &ds, const IDQuint &idquint)
{
    ds << idquint.atom0 << idquint.atom1 << idquint.atom2 << idquint.atom3 << idquint.atom4;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, IDQuint &idquint)
{
    ds >> idquint.atom0 >> idquint.atom1 >> idquint.atom2 >> idquint.atom3 >> idquint.atom4;

    return ds;
}

IDQuint::IDQuint(quint32 atm0, quint32 atm1, quint32 atm2, quint32 atm3, quint32 atm4)
    : atom0(atm0), atom1(atm1), atom2(atm2), atom3(atm3), atom4(atm4)
{
    // this is reversible - reverse to a canonical form
    if (atm0 > atm4)
    {
        qSwap(atom0, atom4);
        qSwap(atom1, atom3);
    }
    else if (atm0 == atm4 and atm1 > atm3)
    {
        qSwap(atom1, atom3);
    }
}

IDQuint::IDQuint(const IDQuint &other)
    : atom0(other.atom0), atom1(other.atom1), atom2(other.atom2), atom3(other.atom3), atom4(other.atom4)
{
}

IDQuint::~IDQuint()
{
}

IDQuint &IDQuint::operator=(const IDQuint &other)
{
    atom0 = other.atom0;
    atom1 = other.atom1;
    atom2 = other.atom2;
    atom3 = other.atom3;
    atom4 = other.atom4;

    return *this;
}

bool IDQuint::operator==(const IDQuint &other) const
{
    return atom0 == other.atom0 and atom1 == other.atom1 and atom2 == other.atom2 and atom3 == other.atom3 and atom4 == other.atom4;
}

bool IDQuint::operator!=(const IDQuint &other) const
{
    return not IDQuint::operator==(other);
}

bool IDQuint::operator>(const IDQuint &other) const
{
    return atom0 > other.atom0 or
           (atom0 == other.atom0 and
            (atom1 > other.atom1 or
             (atom1 == other.atom1 and
              (atom2 > other.atom2 or
               (atom2 == other.atom2 and
                (atom3 > other.atom3 or
                 (atom3 == other.atom3 and atom4 > other.atom4)))))));
}

bool IDQuint::operator>=(const IDQuint &other) const
{
    return IDQuint::operator>(other) or IDQuint::operator==(other);
}

bool IDQuint::operator<(const IDQuint &other) const
{
    return not IDQuint::operator>=(other);
}

bool IDQuint::operator<=(const IDQuint &other) const
{
    return not IDQuint::operator>(other);
}

//////
////// Implementation of CMAPFunctions
//////

static const RegisterMetaType<CMAPFunctions> r_cmapfuncs;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CMAPFunctions &cmapfuncs)
{
    writeHeader(ds, r_cmapfuncs, 1);

    SharedDataStream sds(ds);

    sds << cmapfuncs.molinfo << cmapfuncs.parameters_by_atoms << static_cast<const MoleculeProperty &>(cmapfuncs);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CMAPFunctions &cmapfuncs)
{
    VersionID v = readHeader(ds, r_cmapfuncs);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> cmapfuncs.molinfo >>
            cmapfuncs.parameters_by_atoms >>
            static_cast<MoleculeProperty &>(cmapfuncs);
    }
    else
        throw version_error(v, "1", r_cmapfuncs, CODELOC);

    return ds;
}

CMAPFunctions::CMAPFunctions()
{
}

CMAPFunctions::CMAPFunctions(const MoleculeData &moldata)
    : ConcreteProperty<CMAPFunctions, MoleculeProperty>(),
      molinfo(moldata.info())
{
}

CMAPFunctions::CMAPFunctions(const MoleculeInfoData &molinfo)
    : ConcreteProperty<CMAPFunctions, MoleculeProperty>(),
      molinfo(molinfo)
{
}

CMAPFunctions::CMAPFunctions(const CMAPFunctions &other)
    : ConcreteProperty<CMAPFunctions, MoleculeProperty>(other),
      molinfo(other.molinfo),
      parameters_by_atoms(other.parameters_by_atoms)
{
}

CMAPFunctions::~CMAPFunctions()
{
}

CMAPFunctions &CMAPFunctions::operator=(const CMAPFunctions &other)
{
    MoleculeProperty::operator=(other);
    parameters_by_atoms = other.parameters_by_atoms;

    return *this;
}

bool CMAPFunctions::operator==(const CMAPFunctions &other) const
{
    return MoleculeProperty::operator==(other) and parameters_by_atoms == other.parameters_by_atoms;
}

bool CMAPFunctions::operator!=(const CMAPFunctions &other) const
{
    return MoleculeProperty::operator!=(other) or parameters_by_atoms != other.parameters_by_atoms;
}

QString CMAPFunctions::toString() const
{
    return QString("CMAPFunctions(%1)").arg(parameters_by_atoms.size());
}

CMAPFunctions *CMAPFunctions::clone() const
{
    return new CMAPFunctions(*this);
}

const char *CMAPFunctions::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CMAPFunctions>());
}

const char *CMAPFunctions::what() const
{
    return CMAPFunctions::typeName();
}

const MoleculeInfoData &CMAPFunctions::info() const
{
    return *molinfo;
}

void CMAPFunctions::set(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2,
                        AtomIdx atom3, AtomIdx atom4, const CMAPParameter &parameter)
{
    quint32 atm0 = atom0.map(info().nAtoms());
    quint32 atm1 = atom1.map(info().nAtoms());
    quint32 atm2 = atom2.map(info().nAtoms());
    quint32 atm3 = atom3.map(info().nAtoms());
    quint32 atm4 = atom4.map(info().nAtoms());

    if (atm0 == atm1 or atm0 == atm2 or atm0 == atm3 or atm0 == atom4 or
        atm1 == atm2 or atm1 == atm3 or atm1 == atm4 or
        atm2 == atm3 or atm2 == atm4 or atm3 == atm4)
        throw SireMol::duplicate_atom(QObject::tr("You cannot add a function that acts between the same atoms! "
                                                  "(%1-%2-%3-%4-%5)")
                                          .arg(atm0)
                                          .arg(atm1)
                                          .arg(atm2)
                                          .arg(atm3)
                                          .arg(atm4),
                                      CODELOC);

    if (parameter.isNull())
    {
        parameters_by_atoms.remove(IDQuint(atm0, atm1, atm2, atm3, atm4));
    }
    else
    {
        parameters_by_atoms.insert(IDQuint(atm0, atm1, atm2, atm3, atm4), parameter);
    }
}

void CMAPFunctions::set(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2,
                        const AtomID &atom3, const AtomID &atom4, const CMAPParameter &parameter)
{
    this->set(info().atomIdx(atom0), info().atomIdx(atom1), info().atomIdx(atom2),
              info().atomIdx(atom3), info().atomIdx(atom4), parameter);
}

void CMAPFunctions::set(const CMAPFunction &param)
{
    this->set(info().atomIdx(param.atom()), info().atomIdx(param.atom1()), info().atomIdx(param.atom2()),
              info().atomIdx(param.atom3()), info().atomIdx(param.atom4()), param.parameter());
}

void CMAPFunctions::clear()
{
    parameters_by_atoms.clear();
}

void CMAPFunctions::clear(AtomIdx atom)
{
    quint32 atm = atom.map(info().nAtoms());

    QList<IDQuint> keys = parameters_by_atoms.keys();

    foreach (const IDQuint &key, keys)
    {
        if (key.atom0 == atm or key.atom1 == atm or key.atom2 == atm or key.atom3 == atm or key.atom4 == atm)
        {
            parameters_by_atoms.remove(key);
        }
    }
}

void CMAPFunctions::clear(const AtomID &atom)
{
    QList<AtomIdx> atomidxs = atom.map(info());

    foreach (AtomIdx atomidx, atomidxs)
    {
        this->clear(atomidx);
    }
}

/** Clear all functions that involve any of the atoms in 'atoms'
 *  - if 'exclusive' is true, then this only removes functions
 *  that exclusively involve these atoms - if false, then
 *  if removes functions that involve any of these atoms
 */
void CMAPFunctions::clear(const QList<AtomIdx> &atoms, bool exclusive)
{
    QSet<quint32> atms;
    atms.reserve(atoms.count());

    for (const auto &atom : atoms)
    {
        atms.insert(atom.map(info().nAtoms()));
    }

    QList<IDQuint> keys = parameters_by_atoms.keys();

    if (exclusive)
    {
        for (const auto &key : keys)
        {
            if (atms.contains(key.atom0) and atms.contains(key.atom1) and atms.contains(key.atom2) and atms.contains(key.atom3) and atms.contains(key.atom4))
            {
                parameters_by_atoms.remove(key);
            }
        }
    }
    else
    {
        for (const auto &key : keys)
        {
            if (atms.contains(key.atom0) or atms.contains(key.atom1) or atms.contains(key.atom2) or atms.contains(key.atom3) or atms.contains(key.atom4))
            {
                parameters_by_atoms.remove(key);
            }
        }
    }
}

void CMAPFunctions::clear(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3, AtomIdx atom4)
{
    quint32 atm0 = atom0.map(info().nAtoms());
    quint32 atm1 = atom1.map(info().nAtoms());
    quint32 atm2 = atom2.map(info().nAtoms());
    quint32 atm3 = atom3.map(info().nAtoms());
    quint32 atm4 = atom4.map(info().nAtoms());

    parameters_by_atoms.remove(IDQuint(atm0, atm1, atm2, atm3, atm4));
}

void CMAPFunctions::clear(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2, const AtomID &atom3, const AtomID &atom4)
{
    QList<AtomIdx> atoms0 = atom0.map(info());
    QList<AtomIdx> atoms1 = atom1.map(info());
    QList<AtomIdx> atoms2 = atom2.map(info());
    QList<AtomIdx> atoms3 = atom3.map(info());
    QList<AtomIdx> atoms4 = atom4.map(info());

    foreach (AtomIdx atm0, atoms0)
    {
        foreach (AtomIdx atm1, atoms1)
        {
            foreach (AtomIdx atm2, atoms2)
            {
                foreach (AtomIdx atm3, atoms3)
                {
                    foreach (AtomIdx atm4, atoms4)
                    {
                        this->clear(atm0, atm1, atm2, atm3, atm4);
                    }
                }
            }
        }
    }
}

bool CMAPFunctions::isEmpty() const
{
    return parameters_by_atoms.isEmpty();
}

int CMAPFunctions::nFunctions() const
{
    return parameters_by_atoms.size();
}

CMAPParameter CMAPFunctions::parameter(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3, AtomIdx atom4) const
{
    quint32 atm0 = atom0.map(info().nAtoms());
    quint32 atm1 = atom1.map(info().nAtoms());
    quint32 atm2 = atom2.map(info().nAtoms());
    quint32 atm3 = atom3.map(info().nAtoms());
    quint32 atm4 = atom4.map(info().nAtoms());

    return parameters_by_atoms.value(IDQuint(atm0, atm1, atm2, atm3, atm4));
}

CMAPParameter CMAPFunctions::parameter(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2, const AtomID &atom3, const AtomID &atom4) const
{
    return this->parameter(info().atomIdx(atom0), info().atomIdx(atom1), info().atomIdx(atom2),
                           info().atomIdx(atom3), info().atomIdx(atom4));
}

QVector<CMAPFunction> CMAPFunctions::parameters() const
{
    QVector<CMAPFunction> funcs(parameters_by_atoms.count());

    CMAPFunction *funcs_array = funcs.data();

    int i = 0;

    for (QHash<IDQuint, CMAPParameter>::const_iterator it = parameters_by_atoms.constBegin();
         it != parameters_by_atoms.constEnd(); ++it)
    {
        funcs_array[i] = CMAPFunction(
            info().cgAtomIdx(AtomIdx(it.key().atom0)), info().cgAtomIdx(AtomIdx(it.key().atom1)),
            info().cgAtomIdx(AtomIdx(it.key().atom2)), info().cgAtomIdx(AtomIdx(it.key().atom3)),
            info().cgAtomIdx(AtomIdx(it.key().atom4)), it.value());

        ++i;
    }

    return funcs;
}

QVector<CMAPFunction> CMAPFunctions::parameters(const QList<AtomIdx> &atms, bool exclusive) const
{
    QVector<CMAPFunction> funcs;
    funcs.reserve(parameters_by_atoms.count());

    QSet<AtomIdx> atoms(atms.begin(), atms.end());

    for (QHash<IDQuint, CMAPParameter>::const_iterator it = parameters_by_atoms.constBegin();
         it != parameters_by_atoms.constEnd(); ++it)
    {
        if (exclusive)
        {
            if (atoms.contains(AtomIdx(it.key().atom0)) and atoms.contains(AtomIdx(it.key().atom1)) and atoms.contains(AtomIdx(it.key().atom2)) and atoms.contains(AtomIdx(it.key().atom3)) and atoms.contains(AtomIdx(it.key().atom4)))
            {
                funcs.append(CMAPFunction(info().cgAtomIdx(AtomIdx(it.key().atom0)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom1)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom2)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom3)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom4)),
                                          it.value()));
            }
        }
        else
        {
            if (atoms.contains(AtomIdx(it.key().atom0)) or atoms.contains(AtomIdx(it.key().atom1)) or atoms.contains(AtomIdx(it.key().atom2)) or atoms.contains(AtomIdx(it.key().atom3)) or atoms.contains(AtomIdx(it.key().atom4)))
            {
                funcs.append(CMAPFunction(info().cgAtomIdx(AtomIdx(it.key().atom0)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom1)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom2)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom3)),
                                          info().cgAtomIdx(AtomIdx(it.key().atom4)),
                                          it.value()));
            }
        }
    }

    return funcs;
}

/** Return the set of functions where only functions that involve the
    atoms in 'selected_atoms' are included. If 'isstrict' is true, then
    only include functions where all of the atoms are in 'selected_atoms',
    while if 'isstrict' is false, include functions where at least one
    atom is in 'selected_atoms' */
CMAPFunctions CMAPFunctions::includeOnly(const AtomSelection &selected_atoms, bool isstrict) const
{
    CMAPFunctions ret(*this);

    QMutableHashIterator<IDQuint, CMAPParameter> it(ret.parameters_by_atoms);

    if (isstrict)
    {
        while (it.hasNext())
        {
            it.next();

            if (not(selected_atoms.selected(AtomIdx(it.key().atom0)) and
                    selected_atoms.selected(AtomIdx(it.key().atom1)) and
                    selected_atoms.selected(AtomIdx(it.key().atom2)) and
                    selected_atoms.selected(AtomIdx(it.key().atom3)) and
                    selected_atoms.selected(AtomIdx(it.key().atom4))))
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
                    selected_atoms.selected(AtomIdx(it.key().atom1)) or
                    selected_atoms.selected(AtomIdx(it.key().atom2)) or
                    selected_atoms.selected(AtomIdx(it.key().atom3)) or
                    selected_atoms.selected(AtomIdx(it.key().atom4))))
            {
                it.remove();
            }
        }
    }

    return ret;
}

/** Return whether or not this property is compatible with the molecule
    whose layout information is in 'molinfo' */
bool CMAPFunctions::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return info() == molinfo;
}

/** Merge this property with another property */
PropertyList CMAPFunctions::merge(const MolViewProperty &other,
                                  const AtomIdxMapping &mapping,
                                  const QString &ghost,
                                  const SireBase::PropertyMap &map) const
{
    if (not other.isA<CMAPFunctions>())
    {
        throw SireError::incompatible_error(QObject::tr("Cannot merge %1 with %2 as they are different types.")
                                                .arg(this->what())
                                                .arg(other.what()),
                                            CODELOC);
    }

    if (not ghost.isEmpty())
    {
        Console::warning(QObject::tr("The ghost parameter '%1' for CMAP parameters is ignored").arg(ghost));
    }

    const CMAPFunctions &ref = *this;
    const CMAPFunctions &pert = other.asA<CMAPFunctions>();

    CMAPFunctions prop0 = ref;
    CMAPFunctions prop1 = ref;

    // the prop1 properties are made by finding all of the atoms that
    // are involved in cmaps in 'pert' and removing any involving
    // only those atoms from 'prop1', and then adding back the matching
    // cmap from 'pert'. Use 'true' to only remove cmaps where all
    // atoms are in the mapping
    prop1.clear(mapping.mappedIn1(), true);

    // get the mapping from the perturbed to reference states, including
    // atoms that don't exist in the reference state. In all cases,
    // the values are the indexes in the merged molecule
    auto map1to0 = mapping.map1to0(true);

    // now find all of the cmaps in 'pert' where all atoms in the
    // cmap are in map1to0.keys() - i.e. exist and are mapped from
    // the perturbed state
    const auto pert_cmaps = pert.parameters(map1to0.keys(), true);

    for (const auto &pert_cmap : pert_cmaps)
    {
        const auto atom0 = map1to0.value(info().atomIdx(pert_cmap.atom0()));
        const auto atom1 = map1to0.value(info().atomIdx(pert_cmap.atom1()));
        const auto atom2 = map1to0.value(info().atomIdx(pert_cmap.atom2()));
        const auto atom3 = map1to0.value(info().atomIdx(pert_cmap.atom3()));
        const auto atom4 = map1to0.value(info().atomIdx(pert_cmap.atom4()));

        prop1.set(atom0, atom1, atom2, atom3, atom4, pert_cmap.parameter());

        if (mapping.isUnmappedIn0(atom0) or mapping.isUnmappedIn0(atom1) or mapping.isUnmappedIn0(atom2) or mapping.isUnmappedIn0(atom3) or mapping.isUnmappedIn0(atom4))
        {
            // the prop0 properties are nearly correct - we just need to add
            // in cmaps from 'pert' that involve the atoms that are not mapped
            // in the reference state - this way, those added atoms are held
            // by a constant cmap potential, so won't fly away in the
            // simulation of the reference state
            prop0.set(atom0, atom1, atom2, atom3, atom4, pert_cmap.parameter());
        }
    }

    // now add in the cmaps to the perturbed state from the reference
    // state for any atoms that aren't mapped to the perturbed state.
    // This way, the removed atoms are held by a constant potential,
    // so won't fly away in the simulation of the perturbed state
    auto map0to1 = mapping.map0to1(true);

    const auto ref_cmaps = prop0.parameters(map0to1.keys(), true);

    for (const auto &ref_cmap : ref_cmaps)
    {
        const auto atom0 = info().atomIdx(ref_cmap.atom0());
        const auto atom1 = info().atomIdx(ref_cmap.atom1());
        const auto atom2 = info().atomIdx(ref_cmap.atom2());
        const auto atom3 = info().atomIdx(ref_cmap.atom3());
        const auto atom4 = info().atomIdx(ref_cmap.atom4());

        if (mapping.isUnmappedIn1(atom0) or mapping.isUnmappedIn1(atom1) or mapping.isUnmappedIn1(atom2) or mapping.isUnmappedIn1(atom3) or mapping.isUnmappedIn1(atom4))
        {
            prop1.set(atom0, atom1, atom2, atom3, atom4, ref_cmap.parameter());
        }
    }

    SireBase::PropertyList ret;

    ret.append(prop0);
    ret.append(prop1);

    return ret;
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom matching
    functions in 'atommatcher' to match atoms from the current molecule
    to the atoms in the molecule whose layout is in 'molinfo'

    This only copies the CMAPFunction for pairs of atoms that
    are successfully matched - it does not copy functions for atoms
    that are not matched. Use CMAPFunctions::nFunctions() to check
    if the number of functions in the returned set is the same as
    the number in this set, if you want to ensure that all of the
    functions have been copied.

    \throw SireError::incompatible_error
*/
PropertyPtr CMAPFunctions::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                   const AtomMatcher &atommatcher) const
{
    if (not atommatcher.changesOrder(this->info(), molinfo))
    {
        // the order of the atoms remains the same - this means that the
        // AtomIdx indicies are still valid
        CMAPFunctions ret(molinfo);
        ret.parameters_by_atoms = this->parameters_by_atoms;
        return ret;
    }

    QHash<AtomIdx, AtomIdx> matched_atoms = atommatcher.match(this->info(), molinfo);

    return this->_pvt_makeCompatibleWith(molinfo, matched_atoms);
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom mapping
    in 'map' to match atoms from the current molecule to the atoms in
    the molecule whose layout is in 'molinfo'

    This only copies the CMAPFunction for pairs of atoms that
    are successfully matched - it does not copy functions for atoms
    that are not matched. Use CMAPFunctions::nFunctions() to check
    if the number of functions in the returned set is the same as
    the number in this set, if you want to ensure that all of the
    functions have been copied.

    \throw SireError::incompatible_error
*/
PropertyPtr CMAPFunctions::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                   const QHash<AtomIdx, AtomIdx> &map) const
{
    CMAPFunctions ret(molinfo);

    for (QHash<IDQuint, CMAPParameter>::const_iterator it = parameters_by_atoms.constBegin();
         it != parameters_by_atoms.constEnd(); ++it)
    {
        AtomIdx new_atom0 = map.value(AtomIdx(it.key().atom0), AtomIdx(-1));
        AtomIdx new_atom1 = map.value(AtomIdx(it.key().atom1), AtomIdx(-1));
        AtomIdx new_atom2 = map.value(AtomIdx(it.key().atom2), AtomIdx(-1));
        AtomIdx new_atom3 = map.value(AtomIdx(it.key().atom3), AtomIdx(-1));
        AtomIdx new_atom4 = map.value(AtomIdx(it.key().atom4), AtomIdx(-1));

        if (new_atom0 == -1 or new_atom1 == -1 or new_atom2 == -1 or new_atom3 == -1 or new_atom4 == -1)
            continue;

        ret.set(new_atom0, new_atom1, new_atom2, new_atom3, new_atom4, it.value());
    }

    return ret;
}
