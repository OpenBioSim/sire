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

#include "selectormol.h"
#include "evaluatorm.h"
#include "moleculegroup.h"
#include "moleculegroups.h"
#include "molecules.h"
#include "selectorm.hpp"
#include "selector.hpp"

#include "atomid.h"
#include "cgid.h"
#include "chainid.h"
#include "molid.h"
#include "molidx.h"
#include "molname.h"
#include "molnum.h"
#include "resid.h"
#include "segid.h"

#include "SireID/index.h"

#include "SireBase/lazyevaluator.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireMol;
using namespace SireID;

RegisterMetaType<SelectorMol> r_smol;
RegisterMetaType<SelectorM<Atom>> r_smatm;
RegisterMetaType<SelectorM<Residue>> r_smres;
RegisterMetaType<SelectorM<Chain>> r_smchn;
RegisterMetaType<SelectorM<Segment>> r_smseg;
RegisterMetaType<SelectorM<CutGroup>> r_smcg;
RegisterMetaType<Selector<Atom>> r_satm;
RegisterMetaType<Selector<Residue>> r_sres;
RegisterMetaType<Selector<Chain>> r_schn;
RegisterMetaType<Selector<Segment>> r_sseg;
RegisterMetaType<Selector<CutGroup>> r_scg;

/** Serialise to a binary datastream */
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorMol &mols)
{
    writeHeader(ds, r_smol, 1);

    SharedDataStream sds(ds);

    sds << mols.mols << static_cast<const Property &>(mols);

    return ds;
}

/** Extract from a binary datastream */
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorMol &mols)
{
    VersionID v = readHeader(ds, r_smol);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> mols.mols >> static_cast<Property &>(mols);
    }
    else
        throw version_error(v, "1", r_smol, CODELOC);

    return ds;
}

/** Serialise to a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorM<Atom> &mols)
{
    writeHeader(ds, r_satm, 1);

    SharedDataStream sds(ds);

    sds << mols.vws << static_cast<const Property &>(mols);

    return ds;
}

/** Extract from a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorM<Atom> &mols)
{
    VersionID v = readHeader(ds, r_satm);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> mols.vws >> static_cast<Property &>(mols);
    }
    else
        throw version_error(v, "1", r_satm, CODELOC);

    return ds;
}
/** Serialise to a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorM<Residue> &mols)
{
    writeHeader(ds, r_sres, 1);

    SharedDataStream sds(ds);

    sds << mols.vws << static_cast<const Property &>(mols);

    return ds;
}

/** Extract from a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorM<Residue> &mols)
{
    VersionID v = readHeader(ds, r_sres);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> mols.vws >> static_cast<Property &>(mols);
    }
    else
        throw version_error(v, "1", r_sres, CODELOC);

    return ds;
}
/** Serialise to a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorM<Chain> &mols)
{
    writeHeader(ds, r_schn, 1);

    SharedDataStream sds(ds);

    sds << mols.vws << static_cast<const Property &>(mols);

    return ds;
}

/** Extract from a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorM<Chain> &mols)
{
    VersionID v = readHeader(ds, r_schn);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> mols.vws >> static_cast<Property &>(mols);
    }
    else
        throw version_error(v, "1", r_schn, CODELOC);

    return ds;
}
/** Serialise to a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorM<Segment> &mols)
{
    writeHeader(ds, r_sseg, 1);

    SharedDataStream sds(ds);

    sds << mols.vws << static_cast<const Property &>(mols);

    return ds;
}

/** Extract from a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorM<Segment> &mols)
{
    VersionID v = readHeader(ds, r_sseg);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> mols.vws >> static_cast<Property &>(mols);
    }
    else
        throw version_error(v, "1", r_sseg, CODELOC);

    return ds;
}
/** Serialise to a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorM<CutGroup> &mols)
{
    writeHeader(ds, r_scg, 1);

    SharedDataStream sds(ds);

    sds << mols.vws << static_cast<const Property &>(mols);

    return ds;
}

/** Extract from a binary datastream */
template <>
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorM<CutGroup> &mols)
{
    VersionID v = readHeader(ds, r_scg);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> mols.vws >> static_cast<Property &>(mols);
    }
    else
        throw version_error(v, "1", r_scg, CODELOC);

    return ds;
}

SelectorMol::SelectorMol() : ConcreteProperty<SelectorMol, Property>()
{
}

SelectorMol::SelectorMol(const MoleculeView &mol) : ConcreteProperty<SelectorMol, Property>()
{
    this->mols.append(mol.molecule());
}

SelectorMol::SelectorMol(const Molecules &molecules) : ConcreteProperty<SelectorMol, Property>()
{
    if (not molecules.isEmpty())
    {
        auto toList = [](const QSet<MolNum> &molnums)
        { return molnums.values(); };

        auto molnums = toList(molecules.molNums());

        // sort them, as this is also likely the order the molecules
        // were read in from a file, and so more likely to be the
        // order the user would expect
        std::sort(molnums.begin(), molnums.end());

        this->mols.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            this->mols.append(molecules.at(molnum).molecule());
        }
    }
}

SelectorMol::SelectorMol(const QList<Molecule> &molecules) : ConcreteProperty<SelectorMol, Property>()
{
    mols = molecules;
}

SelectorMol::SelectorMol(const QVector<Molecule> &molecules) : ConcreteProperty<SelectorMol, Property>()
{
    mols = QList<Molecule>(molecules.constBegin(), molecules.constEnd());
}

SelectorMol::SelectorMol(const MoleculeGroup &molecules) : ConcreteProperty<SelectorMol, Property>()
{
    if (not molecules.isEmpty())
    {
        const auto molnums = molecules.molNums();
        this->mols.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            this->mols.append(molecules.at(molnum).molecule());
        }
    }
}

SelectorMol::SelectorMol(const MolGroupsBase &molecules) : ConcreteProperty<SelectorMol, Property>()
{
    if (not molecules.isEmpty())
    {
        const auto molnums = molecules.molNums();
        this->mols.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            this->mols.append(molecules.at(molnum).molecule());
        }
    }
}

SelectorMol::SelectorMol(const SelectResult &molecules) : ConcreteProperty<SelectorMol, Property>()
{
    if (not molecules.isEmpty())
    {
        this->mols.reserve(molecules.count());

        for (const auto &mol : molecules)
        {
            this->mols.append(mol->molecule());
        }
    }
}

SelectorMol::SelectorMol(const QList<SelectorMol> &other) : ConcreteProperty<SelectorMol, Property>()
{
    for (const auto &o : other)
    {
        this->mols += o.mols;
    }
}

SelectorMol::SelectorMol(const QVector<SelectorMol> &other) : ConcreteProperty<SelectorMol, Property>()
{
    this->operator=(SelectorMol(QList<SelectorMol>(other.constBegin(), other.constEnd())));
}

SelectorMol::SelectorMol(const SelectorMol &other) : ConcreteProperty<SelectorMol, Property>(), mols(other.mols)
{
}

SelectorMol::~SelectorMol()
{
}

const char *SelectorMol::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorMol>());
}

SelectorMol &SelectorMol::operator=(const SelectorMol &other)
{
    if (this != &other)
    {
        mols = other.mols;
        Property::operator=(other);
    }

    return *this;
}

bool SelectorMol::operator==(const SelectorMol &other) const
{
    return mols == other.mols;
}

bool SelectorMol::operator!=(const SelectorMol &other) const
{
    return not SelectorMol::operator==(other);
}

Molecule SelectorMol::operator[](int i) const
{
    return this->molecule(i);
}

SelectorMol SelectorMol::operator[](const Slice &slice) const
{
    return this->molecules(slice);
}

SelectorMol SelectorMol::operator[](const QList<qint64> &idxs) const
{
    return this->molecules(idxs);
}

Molecule SelectorMol::operator[](const QString &name) const
{
    return this->molecule(name);
}

Molecule SelectorMol::operator[](const MolIdx &molid) const
{
    return this->molecule(molid);
}

Molecule SelectorMol::operator[](const MolName &molid) const
{
    return this->molecule(molid);
}

Molecule SelectorMol::operator[](const MolNum &molid) const
{
    return this->molecule(molid);
}

Molecule SelectorMol::operator[](const MolID &molid) const
{
    return this->molecule(molid);
}

SelectorMol SelectorMol::operator+(const Molecule &mol) const
{
    SelectorMol ret(*this);
    ret += mol;
    return ret;
}

SelectorMol SelectorMol::operator+(const SelectorMol &mols) const
{
    SelectorMol ret(*this);
    ret += mols;
    return ret;
}

SelectorMol &SelectorMol::operator+=(const Molecule &mol)
{
    this->append(mol);
    return *this;
}

SelectorMol &SelectorMol::operator+=(const SelectorMol &mols)
{
    this->append(mols);
    return *this;
}

void SelectorMol::append(const Molecule &mol)
{
    this->mols.append(mol);
}

void SelectorMol::append(const SelectorMol &other)
{
    this->mols.append(other.mols);
}

QList<MolViewPtr> SelectorMol::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(mols.count());

    for (const auto &mol : mols)
    {
        l.append(MolViewPtr(mol.clone()));
    }

    return l;
}

Molecules SelectorMol::toMolecules() const
{
    return Molecules(mols);
}

SelectResult SelectorMol::toSelectResult() const
{
    return SelectResult(this->mols);
}

SelectorMol SelectorMol::extract() const
{
    return *this;
}

bool SelectorMol::isSelector() const
{
    return true;
}

/** Return the index of 'view' in this list of molecules.
 *  Returns an empty list if the view is not in this collection
 */
QList<qint64> SelectorMol::find(const Molecule &view) const
{
    QList<qint64> matches;

    for (int i = 0; i < mols.count(); ++i)
    {
        if (mols.at(i).isSameMolecule(view))
        {
            matches.append(i);
            break;
        }
    }

    return matches;
}

/** Return the index of 'views' in this list of molecules.
 *  Returns an empty list if none of the views are in this collection
 */
QList<qint64> SelectorMol::find(const SelectorMol &views) const
{
    QList<qint64> matches;

    for (int i = 0; i < views.count(); ++i)
    {
        const auto &view = views[i];

        for (int j = 0; j < mols.count(); ++j)
        {
            if (mols.at(j).isSameMolecule(view))
            {
                matches.append(j);
                break;
            }
        }
    }

    return matches;
}

void SelectorMol::update(const Molecules &molecules)
{
    // better to create a map from MolNum to index here
    QMultiHash<MolNum, int> molnum_to_idx;
    molnum_to_idx.reserve(mols.count());

    int i = 0;

    for (const auto &mol : mols)
    {
        molnum_to_idx.insert(mol.number(), i);
        i += 1;
    }

    for (const auto &mol : molecules)
    {
        const auto molnum = mol.data().number();

        auto it = molnum_to_idx.constFind(molnum);

        while (it != molnum_to_idx.constEnd() && it.key() == molnum)
        {
            mols[it.value()].update(mol.data());
            ++it;
        }
    }
}

void SelectorMol::update(const MoleculeView &molview)
{
    this->update(molview.data());
}

void SelectorMol::update(const MoleculeData &moldata)
{
    const auto molnum = moldata.number();

    QList<int> idx;

    int i = 0;
    for (const auto &mol : mols)
    {
        if (mol.number() == molnum)
            idx.append(i);

        i += 1;
    }

    for (auto i : idx)
    {
        mols[i].update(moldata);
    }
}

void SelectorMol::update(const SelectorMol &molecules)
{
    this->update(molecules.toMolecules());
}

EvaluatorM SelectorMol::evaluate() const
{
    return EvaluatorM(*this);
}

SelectResult SelectorMol::search(const QString &search_string) const
{
    Select search(search_string);
    return search(this->toSelectResult());
}

SelectResult SelectorMol::search(const QString &search_string, const PropertyMap &map) const
{
    Select search(search_string);
    return search(this->toSelectResult(), map);
}

int SelectorMol::count() const
{
    return this->mols.count();
}

int SelectorMol::size() const
{
    return this->mols.count();
}

bool SelectorMol::contains(const MolNum &molnum) const
{
    for (const auto &mol : this->mols)
    {
        if (mol.number() == molnum)
            return true;
    }

    return false;
}

Molecule SelectorMol::molecule(int i, const PropertyMap &) const
{
    return this->mols.at(Index(i).map(this->mols.count()));
}

Molecule SelectorMol::molecule(const QString &name, const PropertyMap &map) const
{
    auto m = this->molecules(name, map);

    if (m.count() > 1)
        throw SireMol::duplicate_molecule(
            QObject::tr("More than one molecule matches %1. Number of matches is %2.").arg(name).arg(m.count()),
            CODELOC);

    BOOST_ASSERT(m.count() != 0);

    return m[0];
}

Molecule SelectorMol::molecule(const MolIdx &molid, const PropertyMap &map) const
{
    auto m = this->molecules(molid, map);

    if (m.count() > 1)
        throw SireMol::duplicate_molecule(QObject::tr("More than one molecule matches %1. Number of matches is %2.")
                                              .arg(molid.toString())
                                              .arg(m.count()),
                                          CODELOC);

    BOOST_ASSERT(m.count() != 0);

    return m[0];
}

Molecule SelectorMol::molecule(const MolName &molid, const PropertyMap &map) const
{
    auto m = this->molecules(molid, map);

    if (m.count() > 1)
        throw SireMol::duplicate_molecule(QObject::tr("More than one molecule matches %1. Number of matches is %2.")
                                              .arg(molid.toString())
                                              .arg(m.count()),
                                          CODELOC);

    BOOST_ASSERT(m.count() != 0);

    return m[0];
}

Molecule SelectorMol::molecule(const MolNum &molid, const PropertyMap &map) const
{
    auto m = this->molecules(molid, map);

    if (m.count() > 1)
        throw SireMol::duplicate_molecule(QObject::tr("More than one molecule matches %1. Number of matches is %2.")
                                              .arg(molid.toString())
                                              .arg(m.count()),
                                          CODELOC);

    BOOST_ASSERT(m.count() != 0);

    return m[0];
}

Molecule SelectorMol::molecule(const MolID &molid, const PropertyMap &map) const
{
    auto m = this->molecules(molid, map);

    if (m.count() > 1)
        throw SireMol::duplicate_molecule(QObject::tr("More than one molecule matches %1. Number of matches is %2.")
                                              .arg(molid.toString())
                                              .arg(m.count()),
                                          CODELOC);

    BOOST_ASSERT(m.count() != 0);

    return m[0];
}

SelectorMol SelectorMol::molecules() const
{
    return *this;
}

SelectorMol SelectorMol::molecules(int i, const PropertyMap &map) const
{
    return SelectorMol(this->molecule(i, map));
}

SelectorMol SelectorMol::molecules(const SireBase::Slice &slice, const PropertyMap &) const
{
    SelectorMol m;

    for (auto it = slice.begin(this->mols.count()); not it.atEnd(); it.next())
    {
        m.mols.append(this->mols.at(it.value()));
    }

    return m;
}

SelectorMol SelectorMol::molecules(const QList<qint64> &idxs, const PropertyMap &) const
{
    SelectorMol m;

    for (const auto &idx : idxs)
    {
        m.mols.append(this->molecule(idx));
    }

    return m;
}

SelectorMol SelectorMol::molecules(const QString &name, const PropertyMap &map) const
{
    try
    {
        return this->molecules(MolName(name));
    }
    catch (const SireError::exception &e)
    {
        try
        {
            SelectorMol ret(this->search(name, map));

            if (ret.count() == 0)
                throw SireMol::missing_molecule(QObject::tr("No molecule matches '%1'").arg(name), CODELOC);

            return ret;
        }
        catch (...)
        {
            if (name.length() < 5)
                // likely a name error
                e.throwSelf();
            else
                // likely a syntax error
                throw;
        }
    }

    return SelectorMol();
}

SelectorMol SelectorMol::molecules(const MolIdx &molid, const PropertyMap &map) const
{
    return this->molecules(molid.value(), map);
}

SelectorMol SelectorMol::molecules(const MolName &molname, const PropertyMap &) const
{
    SelectorMol m;

    for (const auto &mol : this->mols)
    {
        if (mol.name() == molname)
        {
            m.mols.append(mol);
        }
    }

    if (m.isEmpty())
        throw SireMol::missing_molecule(QObject::tr("There is no molecule called %1 in this container.").arg(molname),
                                        CODELOC);

    return m;
}

SelectorMol SelectorMol::molecules(const MolNum &molnum, const PropertyMap &) const
{
    SelectorMol m;

    for (const auto &mol : this->mols)
    {
        if (mol.number() == molnum)
        {
            m.mols.append(mol);
            // there should only be one matching molecule
            break;
        }
    }

    if (m.isEmpty())
        throw SireMol::missing_molecule(
            QObject::tr("There is no molecule with number %1 in this container.").arg(molnum), CODELOC);

    return m;
}

MoleculeGroup SelectorMol::toMoleculeGroup() const
{
    MoleculeGroup grp("all");

    for (const auto &mol : this->mols)
    {
        grp.add(mol);
    }

    return grp;
}

SelectorMol SelectorMol::molecules(const MolID &molid, const PropertyMap &) const
{
    auto molnums = molid.map(this->toMoleculeGroup());

    QHash<MolNum, Molecule> m;

    // QSet<MolNum> molnums_set(molnums.constBegin(), molnums.constEnd());
    QSet<MolNum> molnums_set = _list_to_set(molnums);

    for (const auto &mol : this->mols)
    {
        if (molnums_set.contains(mol.number()))
        {
            if (not m.contains(mol.number()))
            {
                m.insert(mol.number(), mol);

                if (m.count() == molnums_set.count())
                    break;
            }
        }
    }

    SelectorMol ret;

    for (const auto &molnum : molnums)
    {
        ret.mols.append(m[molnum]);
    }

    return ret;
}

Atom SelectorMol::atom(int i, const PropertyMap &) const
{
    i = Index(i).map(this->nAtoms());

    for (const auto &mol : this->mols)
    {
        if (i < mol.nAtoms())
            return mol.atom(AtomIdx(i));
        else
            i -= mol.nAtoms();
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Atom();
}

Atom SelectorMol::atom(const QString &name, const PropertyMap &map) const
{
    auto atms = this->atoms(name, map);

    if (atms.count() > 1)
        throw SireMol::duplicate_atom(QObject::tr("More than one atom matches the name '%1'. Number of matches is %2.")
                                          .arg(name)
                                          .arg(atms.count()),
                                      CODELOC);

    BOOST_ASSERT(not atms.isEmpty());

    return atms[0];
}

Atom SelectorMol::atom(const AtomID &atomid, const PropertyMap &map) const
{
    auto atms = this->atoms(atomid, map);

    if (atms.count() > 1)
        throw SireMol::duplicate_atom(QObject::tr("More than one atom matches '%1'. Number of matches is %2.")
                                          .arg(atomid.toString())
                                          .arg(atms.count()),
                                      CODELOC);

    BOOST_ASSERT(not atms.isEmpty());

    return atms[0];
}

Residue SelectorMol::residue(int i, const PropertyMap &) const
{
    i = Index(i).map(this->nResidues());

    for (const auto &mol : this->mols)
    {
        if (i < mol.nResidues())
            return mol.residue(ResIdx(i));
        else
            i -= mol.nResidues();
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Residue();
}

Residue SelectorMol::residue(const QString &name, const PropertyMap &map) const
{
    auto res = this->residues(name, map);

    if (res.count() > 1)
        throw SireMol::duplicate_residue(
            QObject::tr("More than one residue matches the name '%1'. Number of matches is %2.")
                .arg(name)
                .arg(res.count()),
            CODELOC);

    BOOST_ASSERT(not res.isEmpty());

    return res[0];
}

Residue SelectorMol::residue(const ResID &resid, const PropertyMap &map) const
{
    auto res = this->residues(resid, map);

    if (res.count() > 1)
        throw SireMol::duplicate_residue(QObject::tr("More than one residue matches '%1'. Number of matches is %2.")
                                             .arg(resid.toString())
                                             .arg(res.count()),
                                         CODELOC);

    BOOST_ASSERT(not res.isEmpty());

    return res[0];
}

Chain SelectorMol::chain(int i, const PropertyMap &) const
{
    i = Index(i).map(this->nChains());

    for (const auto &mol : this->mols)
    {
        if (i < mol.nChains())
            return mol.chain(ChainIdx(i));
        else
            i -= mol.nChains();
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Chain();
}

Chain SelectorMol::chain(const QString &name, const PropertyMap &map) const
{
    auto cs = this->chains(name, map);

    if (cs.count() > 1)
        throw SireMol::duplicate_chain(
            QObject::tr("More than one chain matches the name '%1'. Number of matches is %2.")
                .arg(name)
                .arg(cs.count()),
            CODELOC);

    BOOST_ASSERT(not cs.isEmpty());

    return cs[0];
}

Chain SelectorMol::chain(const ChainID &chainid, const PropertyMap &map) const
{
    auto cs = this->chains(chainid, map);

    if (cs.count() > 1)
        throw SireMol::duplicate_chain(QObject::tr("More than one chain matches '%1'. Number of matches is %2.")
                                           .arg(chainid.toString())
                                           .arg(cs.count()),
                                       CODELOC);

    BOOST_ASSERT(not cs.isEmpty());

    return cs[0];
}

Segment SelectorMol::segment(int i, const PropertyMap &) const
{
    i = Index(i).map(this->nSegments());

    for (const auto &mol : this->mols)
    {
        if (i < mol.nSegments())
            return mol.segment(SegIdx(i));
        else
            i -= mol.nSegments();
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Segment();
}

Segment SelectorMol::segment(const QString &name, const PropertyMap &map) const
{
    auto segs = this->segments(name, map);

    if (segs.count() > 1)
        throw SireMol::duplicate_segment(
            QObject::tr("More than one segment matches the name '%1'. Number of matches is %2.")
                .arg(name)
                .arg(segs.count()),
            CODELOC);

    BOOST_ASSERT(not segs.isEmpty());

    return segs[0];
}

Segment SelectorMol::segment(const SegID &segid, const PropertyMap &map) const
{
    auto segs = this->segments(segid, map);

    if (segs.count() > 1)
        throw SireMol::duplicate_segment(QObject::tr("More than one segment matches '%1'. Number of matches is %2.")
                                             .arg(segid.toString())
                                             .arg(segs.count()),
                                         CODELOC);

    BOOST_ASSERT(not segs.isEmpty());

    return segs[0];
}

CutGroup SelectorMol::cutGroup(int i, const PropertyMap &) const
{
    i = Index(i).map(this->nCutGroups());

    for (const auto &mol : this->mols)
    {
        if (i < mol.nCutGroups())
            return mol.cutGroup(CGIdx(i));
        else
            i -= mol.nCutGroups();
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return CutGroup();
}

CutGroup SelectorMol::cutGroup(const QString &name, const PropertyMap &map) const
{
    auto cgs = this->cutGroups(name, map);

    if (cgs.count() > 1)
        throw SireMol::duplicate_cutgroup(
            QObject::tr("More than one CutGroup matches the name '%1'. Number of matches is %2.")
                .arg(name)
                .arg(cgs.count()),
            CODELOC);

    BOOST_ASSERT(not cgs.isEmpty());

    return cgs[0];
}

CutGroup SelectorMol::cutGroup(const CGID &cgid, const PropertyMap &map) const
{
    auto cgs = this->cutGroups(cgid, map);

    if (cgs.count() > 1)
        throw SireMol::duplicate_cutgroup(QObject::tr("More than one CutGroup matches '%1'. Number of matches is %2.")
                                              .arg(cgid.toString())
                                              .arg(cgs.count()),
                                          CODELOC);

    BOOST_ASSERT(not cgs.isEmpty());

    return cgs[0];
}

SelectorM<Atom> SelectorMol::atoms() const
{
    return SelectorM<Atom>(*this);
}

SelectorM<Atom> SelectorMol::atoms(int i, const PropertyMap &map) const
{
    return SelectorM<Atom>(this->atom(i, map));
}

SelectorM<Atom> SelectorMol::atoms(const SireBase::Slice &slice, const PropertyMap &) const
{
    return SelectorM<Atom>(*this, slice);
}

SelectorM<Atom> SelectorMol::atoms(const QList<qint64> &idxs, const PropertyMap &) const
{
    return SelectorM<Atom>(*this, idxs);
}

SelectorM<Atom> SelectorMol::atoms(const QString &name, const PropertyMap &map) const
{
    return SelectorM<Atom>(*this, name, map);
}

SelectorM<Atom> SelectorMol::atoms(const AtomID &atomid, const PropertyMap &) const
{
    return SelectorM<Atom>(*this, atomid);
}

SelectorM<Residue> SelectorMol::residues() const
{
    return SelectorM<Residue>(*this);
}

SelectorM<Residue> SelectorMol::residues(int i, const PropertyMap &map) const
{
    return SelectorM<Residue>(this->residue(i, map));
}

SelectorM<Residue> SelectorMol::residues(const SireBase::Slice &slice, const PropertyMap &) const
{
    return SelectorM<Residue>(*this, slice);
}

SelectorM<Residue> SelectorMol::residues(const QList<qint64> &idxs, const PropertyMap &) const
{
    return SelectorM<Residue>(*this, idxs);
}

SelectorM<Residue> SelectorMol::residues(const QString &name, const PropertyMap &map) const
{
    return SelectorM<Residue>(*this, name, map);
}

SelectorM<Residue> SelectorMol::residues(const ResID &resid, const PropertyMap &) const
{
    return SelectorM<Residue>(*this, resid);
}

SelectorM<Chain> SelectorMol::chains() const
{
    return SelectorM<Chain>(*this);
}

SelectorM<Chain> SelectorMol::chains(int i, const PropertyMap &map) const
{
    return SelectorM<Chain>(this->chain(i, map));
}

SelectorM<Chain> SelectorMol::chains(const SireBase::Slice &slice, const PropertyMap &) const
{
    return SelectorM<Chain>(*this, slice);
}

SelectorM<Chain> SelectorMol::chains(const QList<qint64> &idxs, const PropertyMap &) const
{
    return SelectorM<Chain>(*this, idxs);
}

SelectorM<Chain> SelectorMol::chains(const QString &name, const PropertyMap &map) const
{
    return SelectorM<Chain>(*this, name, map);
}

SelectorM<Chain> SelectorMol::chains(const ChainID &chainid, const PropertyMap &) const
{
    return SelectorM<Chain>(*this, chainid);
}

SelectorM<Segment> SelectorMol::segments() const
{
    return SelectorM<Segment>(*this);
}

SelectorM<Segment> SelectorMol::segments(int i, const PropertyMap &map) const
{
    return SelectorM<Segment>(this->segment(i, map));
}

SelectorM<Segment> SelectorMol::segments(const SireBase::Slice &slice, const PropertyMap &) const
{
    return SelectorM<Segment>(*this, slice);
}

SelectorM<Segment> SelectorMol::segments(const QList<qint64> &idxs, const PropertyMap &) const
{
    return SelectorM<Segment>(*this, idxs);
}

SelectorM<Segment> SelectorMol::segments(const QString &name, const PropertyMap &map) const
{
    return SelectorM<Segment>(*this, name, map);
}

SelectorM<Segment> SelectorMol::segments(const SegID &segid, const PropertyMap &) const
{
    return SelectorM<Segment>(*this, segid);
}

SelectorM<CutGroup> SelectorMol::cutGroups() const
{
    return SelectorM<CutGroup>(*this);
}

SelectorM<CutGroup> SelectorMol::cutGroups(int i, const PropertyMap &map) const
{
    return SelectorM<CutGroup>(this->cutGroup(i, map));
}

SelectorM<CutGroup> SelectorMol::cutGroups(const SireBase::Slice &slice, const PropertyMap &) const
{
    return SelectorM<CutGroup>(*this, slice);
}

SelectorM<CutGroup> SelectorMol::cutGroups(const QList<qint64> &idxs, const PropertyMap &) const
{
    return SelectorM<CutGroup>(*this, idxs);
}

SelectorM<CutGroup> SelectorMol::cutGroups(const QString &name, const PropertyMap &map) const
{
    return SelectorM<CutGroup>(*this, name, map);
}

SelectorM<CutGroup> SelectorMol::cutGroups(const CGID &cgid, const PropertyMap &) const
{
    return SelectorM<CutGroup>(*this, cgid);
}

QVector<MolNum> SelectorMol::molNums() const
{
    if (this->isEmpty())
        return QVector<MolNum>();

    QVector<MolNum> molnums;
    molnums.reserve(this->mols.count());

    for (const auto &mol : this->mols)
    {
        molnums.append(mol.number());
    }

    return molnums;
}

QList<MolNum> SelectorMol::IDs() const
{
    return this->numbers();
}

QList<MolIdx> SelectorMol::indexes() const
{
    if (this->isEmpty())
        return QList<MolIdx>();

    QList<MolIdx> idxs;
    idxs.reserve(this->mols.count());

    for (int i = 0; i < this->mols.count(); ++i)
    {
        idxs.append(MolIdx(i));
    }

    return idxs;
}

QList<MolNum> SelectorMol::numbers() const
{
    if (this->isEmpty())
        return QList<MolNum>();

    QList<MolNum> nums;
    nums.reserve(this->mols.count());

    for (const auto &mol : this->mols)
    {
        nums.append(mol.number());
    }

    return nums;
}

QList<MolName> SelectorMol::names() const
{
    if (this->isEmpty())
        return QList<MolName>();

    QList<MolName> names;
    names.reserve(this->mols.count());

    for (const auto &mol : this->mols)
    {
        names.append(mol.name());
    }

    return names;
}

int SelectorMol::nAtoms() const
{
    int nats = 0;

    for (const auto &mol : this->mols)
    {
        nats += mol.nAtoms();
    }

    return nats;
}

int SelectorMol::nResidues() const
{
    int nres = 0;

    for (const auto &mol : this->mols)
    {
        nres += mol.nResidues();
    }

    return nres;
}

int SelectorMol::nChains() const
{
    int nc = 0;

    for (const auto &mol : this->mols)
    {
        nc += mol.nChains();
    }

    return nc;
}

int SelectorMol::nSegments() const
{
    int ns = 0;

    for (const auto &mol : this->mols)
    {
        ns += mol.nSegments();
    }

    return ns;
}

int SelectorMol::nCutGroups() const
{
    int nc = 0;

    for (const auto &mol : this->mols)
    {
        nc += mol.nCutGroups();
    }

    return nc;
}

int SelectorMol::nMolecules() const
{
    return this->mols.count();
}

int SelectorMol::nFrames() const
{
    return this->nFrames(PropertyMap());
}

int SelectorMol::nFrames(const SireBase::PropertyMap &map) const
{
    return SireMol::detail::_nFrames(this->mols, map);
}

void SelectorMol::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

void SelectorMol::loadFrame(int frame, const LazyEvaluator &evaluator)
{
    this->loadFrame(frame, evaluator, PropertyMap());
}

void SelectorMol::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

void SelectorMol::saveFrame()
{
    this->saveFrame(PropertyMap());
}

void SelectorMol::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

void SelectorMol::deleteAllFrames()
{
    this->deleteAllFrames(PropertyMap());
}

void SelectorMol::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    LazyEvaluator evaluator;
    this->loadFrame(frame, evaluator, map);
}

void SelectorMol::loadFrame(int frame, const LazyEvaluator &evaluator,
                            const SireBase::PropertyMap &map)
{
    SireMol::detail::_loadFrame(this->mols, frame, evaluator, map);
}

void SelectorMol::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->mols, frame, map);
}

void SelectorMol::saveFrame(const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->mols, map);
}

void SelectorMol::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_deleteFrame(this->mols, frame, map);
}

void SelectorMol::deleteAllFrames(const SireBase::PropertyMap &map)
{
    SireMol::detail::_deleteAllFrames(this->mols, map);
}

bool SelectorMol::isEmpty() const
{
    return this->mols.isEmpty();
}

SelectorMol::const_iterator SelectorMol::begin() const
{
    return this->mols.begin();
}

SelectorMol::const_iterator SelectorMol::end() const
{
    return this->mols.end();
}

SelectorMol::const_iterator SelectorMol::constBegin() const
{
    return this->mols.constBegin();
}

SelectorMol::const_iterator SelectorMol::constEnd() const
{
    return this->mols.constEnd();
}

QString SelectorMol::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("SelectorMol::empty");
    }
    else
    {
        QStringList parts;

        const auto n = this->count();

        if (n <= 10)
        {
            for (int i = 0; i < n; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->mols[i].toString()));
            }
        }
        else
        {
            for (int i = 0; i < 5; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->mols[i].toString()));
            }

            parts.append("...");

            for (int i = n - 5; i < n; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->mols[i].toString()));
            }
        }

        return QObject::tr("SelectorMol( size=%1\n%2\n)").arg(n).arg(parts.join("\n"));
    }
}

namespace SireMol
{
    template class SelectorM<Atom>;
    template class SelectorM<Residue>;
    template class SelectorM<Chain>;
    template class SelectorM<Segment>;
    template class SelectorM<CutGroup>;
} // namespace SireMol
