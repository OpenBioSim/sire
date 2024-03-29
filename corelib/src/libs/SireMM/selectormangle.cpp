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

#include "selectormangle.h"

#include "SireID/index.h"

#include "SireCAS/expression.h"

#include "SireBase/parallel.h"

#include "SireBase/errors.h"
#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireMol;
using namespace SireMM;
using namespace SireID;

RegisterMetaType<SelectorMAngle> r_sang;

/** Serialise to a binary datastream */
SIREMM_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorMAngle &angs)
{
    writeHeader(ds, r_sang, 1);

    SharedDataStream sds(ds);

    sds << angs.angs << static_cast<const Property &>(angs);

    return ds;
}

/** Extract from a binary datastream */
SIREMM_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorMAngle &angs)
{
    VersionID v = readHeader(ds, r_sang);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> angs.angs >> static_cast<Property &>(angs);
    }
    else
        throw version_error(v, "1", r_sang, CODELOC);

    return ds;
}

SelectorMAngle::SelectorMAngle() : ConcreteProperty<SelectorMAngle, Property>()
{
}

SelectorMAngle::SelectorMAngle(const Angle &view) : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not view.isEmpty())
        angs.append(SelectorAngle(view));
}

SelectorMAngle::SelectorMAngle(const Molecules &mols, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not mols.isEmpty())
    {
        auto toList = [](const QSet<MolNum> &molnums)
        { return molnums.values(); };

        auto molnums = toList(mols.molNums());

        // sort them, as this is also likely the order the molecules
        // were read in from a file, and so more likely to be the
        // order the user would expect
        std::sort(molnums.begin(), molnums.end());

        this->angs.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorAngle a(mols.at(molnum), map);

            if (not a.isEmpty())
                this->angs.append(a);
        }
    }
}

SelectorMAngle::SelectorMAngle(const MoleculeGroup &mols, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->angs.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorAngle a(mols.at(molnum), map);

            if (not a.isEmpty())
                this->angs.append(a);
        }
    }
}

SelectorMAngle::SelectorMAngle(const MolGroupsBase &mols, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->angs.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorAngle a(mols.at(molnum), map);

            if (not a.isEmpty())
                this->angs.append(a);
        }
    }
}

SelectorMAngle::SelectorMAngle(const SelectResult &mols, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not mols.isEmpty())
    {
        this->angs.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorAngle a;

            if (mol->isA<SelectorAngle>())
                a = mol->asA<SelectorAngle>();
            else
                a = SelectorAngle(*mol, map);

            if (not a.isEmpty())
                this->angs.append(a);
        }
    }
}

SelectorMAngle::SelectorMAngle(const SelectResult &mols, const AngleID &angle, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not mols.isEmpty())
    {
        this->angs.reserve(mols.count());

        for (const auto &mol : mols)
        {
            try
            {
                auto a = SelectorAngle(*mol, angle, map);

                if (not a.isEmpty())
                    this->angs.append(a);
            }
            catch (...)
            {
            }
        }
    }
}

SelectorMAngle::SelectorMAngle(const SelectorAngle &angles) : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not angles.isEmpty())
        angs.append(angles);
}

SelectorMAngle::SelectorMAngle(const SelectorMol &mols, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    if (not mols.isEmpty())
    {
        this->angs.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorAngle a(mol, map);

            if (not a.isEmpty())
                angs.append(a);
        }
    }
}

void SelectorMAngle::_append(const SelectorAngle &angles)
{
    if (angles.isEmpty())
        return;

    if (this->angs.isEmpty())
        this->angs.append(angles);
    else
    {
        for (int i = 0; i < angles.count(); ++i)
        {
            this->_append(angles(i));
        }
    }
}

void SelectorMAngle::_append(const Angle &angle)
{
    if (this->angs.isEmpty())
    {
        this->angs.append(SelectorAngle(angle));
    }
    else if (this->angs.last().data().number() != angle.data().number())
    {
        // new molecule
        this->angs.append(SelectorAngle(angle));
    }
    else
    {
        // a new view in the current molecule
        this->angs.last() = this->angs.last().add(angle);
    }
}

SelectorMAngle::SelectorMAngle(const SelectorMAngle &angles, const SireBase::Slice &slice)
    : SireBase::ConcreteProperty<SelectorMAngle, Property>()
{
    for (auto it = slice.begin(angles.count()); not it.atEnd(); it.next())
    {
        this->_append(angles[it.value()]);
    }
}

SelectorMAngle::SelectorMAngle(const SelectorMAngle &angles, const QList<qint64> &idxs)
    : SireBase::ConcreteProperty<SelectorMAngle, Property>()
{
    for (const auto &idx : idxs)
    {
        this->_append(angles[idx]);
    }
}

SelectorMAngle::SelectorMAngle(const SelectorM<Atom> &atoms, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    for (const auto &mol_atoms : atoms)
    {
        const auto angles = SelectorAngle(mol_atoms, map);
        this->_append(angles);
    }
}

SelectorMAngle::SelectorMAngle(const SelectorM<Atom> &atoms0, const SelectorM<Atom> &atoms1, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    for (const auto &mol_atoms0 : atoms0)
    {
        for (const auto &mol_atoms1 : atoms1)
        {
            if (mol_atoms0.isSameMolecule(mol_atoms1))
            {
                const auto angles = SelectorAngle(mol_atoms0, mol_atoms1, map);
                this->_append(angles);
            }
        }
    }
}

SelectorMAngle::SelectorMAngle(const SelectorM<Atom> &atoms0, const SelectorM<Atom> &atoms1,
                               const SelectorM<Atom> &atoms2, const PropertyMap &map)
    : ConcreteProperty<SelectorMAngle, Property>()
{
    for (const auto &mol_atoms0 : atoms0)
    {
        for (const auto &mol_atoms1 : atoms1)
        {
            if (mol_atoms0.isSameMolecule(mol_atoms1))
            {
                for (const auto &mol_atoms2 : atoms2)
                {
                    if (mol_atoms0.isSameMolecule(mol_atoms2))
                    {
                        const auto angles = SelectorAngle(mol_atoms0, mol_atoms1, mol_atoms2, map);

                        this->_append(angles);
                    }
                }
            }
        }
    }
}

SelectorMAngle::SelectorMAngle(const SelectorMAngle &other)
    : ConcreteProperty<SelectorMAngle, Property>(), angs(other.angs)
{
}

SelectorMAngle::~SelectorMAngle()
{
}

const char *SelectorMAngle::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorMAngle>());
}

SelectorMAngle &SelectorMAngle::operator=(const SelectorMAngle &other)
{
    if (this != &other)
    {
        angs = other.angs;
        Property::operator=(other);
    }

    return *this;
}

bool SelectorMAngle::operator==(const SelectorMAngle &other) const
{
    return angs == other.angs;
}

bool SelectorMAngle::operator!=(const SelectorMAngle &other) const
{
    return not operator==(other);
}

Angle SelectorMAngle::operator[](int i) const
{
    i = SireID::Index(i).map(this->count());

    for (const auto &a : angs)
    {
        if (i < a.count())
        {
            return a(i);
        }
        else
        {
            i -= a.count();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Angle();
}

SelectorMAngle SelectorMAngle::operator[](const SireBase::Slice &slice) const
{
    return SelectorMAngle(*this, slice);
}

SelectorMAngle SelectorMAngle::operator[](const QList<qint64> &idxs) const
{
    return SelectorMAngle(*this, idxs);
}

SelectorMAngle SelectorMAngle::operator[](const AngleID &id) const
{
    SelectorMAngle ret;

    for (const auto &a : angs)
    {
        try
        {
            auto r = a(id);

            if (not r.isEmpty())
            {
                ret.angs.append(r);
            }
        }
        catch (...)
        {
        }
    }

    return ret;
}

Angle SelectorMAngle::operator()(int i) const
{
    return this->operator[](i);
}

SelectorMAngle SelectorMAngle::operator()(const SireBase::Slice &slice) const
{
    return this->operator[](slice);
}

SelectorMAngle SelectorMAngle::operator()(const QList<qint64> &idxs) const
{
    return this->operator[](idxs);
}

SelectorMAngle SelectorMAngle::operator()(const AngleID &id) const
{
    return this->operator[](id);
}

SelectorMol SelectorMAngle::extract() const
{
    const int nmols = this->count();

    const bool uses_parallel = nmols < 16;

    QVector<Molecule> mols(nmols);
    Molecule *mols_data = mols.data();

    if (uses_parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                          {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                mols_data[i] = this->operator()(i).extract();
            } });
    }
    else
    {
        for (int i = 0; i < nmols; ++i)
        {
            mols_data[i] = this->operator()(i).extract();
        }
    }

    return SelectorMol(mols);
}

bool SelectorMAngle::isSelector() const
{
    return true;
}

QList<MolViewPtr> SelectorMAngle::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(angs.count());

    for (const auto &ang : angs)
    {
        l.append(MolViewPtr(ang.clone()));
    }

    return l;
}

Molecules SelectorMAngle::toMolecules() const
{
    return Molecules(this->angs);
}

void SelectorMAngle::update(const Molecules &molecules)
{
    // better to create a map from MolNum to index here
    QMultiHash<MolNum, int> molnum_to_idx;
    molnum_to_idx.reserve(this->angs.count());

    int i = 0;

    for (const auto &mol : this->angs)
    {
        molnum_to_idx.insert(mol.data().number(), i);
        i += 1;
    }

    for (const auto &mol : molecules)
    {
        const auto molnum = mol.data().number();

        auto it = molnum_to_idx.constFind(molnum);

        while (it != molnum_to_idx.constEnd() && it.key() == molnum)
        {
            this->angs[it.value()].update(mol.data());
            ++it;
        }
    }
}

void SelectorMAngle::update(const SelectorMol &molecules)
{
    this->update(molecules.toMolecules());
}

void SelectorMAngle::update(const MoleculeData &moldata)
{
    QList<int> idx;

    const auto molnum = moldata.number();

    int i = 0;
    for (const auto &mol : this->angs)
    {
        if (mol.data().number() == molnum)
        {
            idx.append(i);
        }
        i += 1;
    }

    for (auto i : idx)
    {
        this->angs[i].update(moldata);
    }
}

void SelectorMAngle::update(const MoleculeView &molview)
{
    this->update(molview.data());
}

int SelectorMAngle::count() const
{
    int n = 0;

    for (const auto &a : angs)
    {
        n += a.count();
    }

    return n;
}

int SelectorMAngle::size() const
{
    return this->count();
}

EvaluatorM SelectorMAngle::evaluate() const
{
    return EvaluatorM(this->atoms());
}

MoleculeGroup SelectorMAngle::toMoleculeGroup() const
{
    MoleculeGroup grp;

    for (const auto &a : this->angs)
    {
        grp.add(a);
    }

    return grp;
}

SelectResult SelectorMAngle::toSelectResult() const
{
    QList<MolViewPtr> r;

    for (const auto &a : angs)
    {
        r.append(a);
    }

    return SelectResult(r);
}

Molecule SelectorMAngle::molecule(int i, const PropertyMap &map) const
{
    return this->molecules().molecule(i);
}

Molecule SelectorMAngle::molecule(const QString &name, const PropertyMap &map) const
{
    return this->molecules().molecule(name);
}

Molecule SelectorMAngle::molecule(const MolID &molid, const PropertyMap &map)
{
    return this->molecules().molecule(molid);
}

SelectorMol SelectorMAngle::molecules() const
{
    QList<Molecule> mols;

    for (const auto &a : this->angs)
    {
        mols.append(a.molecule());
    }

    return SelectorMol(mols);
}

SelectorMol SelectorMAngle::molecules(int i, const PropertyMap &map) const
{
    return this->molecules().molecules(i);
}

SelectorMol SelectorMAngle::molecules(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->molecules().molecules(slice);
}

SelectorMol SelectorMAngle::molecules(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->molecules().molecules(idxs);
}

SelectorMol SelectorMAngle::molecules(const QString &name, const PropertyMap &map) const
{
    return this->molecules().molecules(name);
}

SelectorMol SelectorMAngle::molecules(const MolID &molid, const PropertyMap &map) const
{
    return this->molecules().molecules(molid);
}

Atom SelectorMAngle::atom(int i, const PropertyMap &map) const
{
    return this->atoms().atom(i, map);
}

Atom SelectorMAngle::atom(const QString &name, const PropertyMap &map) const
{
    return this->atoms().atom(name, map);
}

Atom SelectorMAngle::atom(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atoms().atom(atomid, map);
}

Residue SelectorMAngle::residue(int i, const PropertyMap &map) const
{
    return this->residues().residue(i, map);
}

Residue SelectorMAngle::residue(const QString &name, const PropertyMap &map) const
{
    return this->residues().residue(name, map);
}

Residue SelectorMAngle::residue(const ResID &resid, const PropertyMap &map) const
{
    return this->residues().residue(resid, map);
}

Chain SelectorMAngle::chain(int i, const PropertyMap &map) const
{
    return this->chains().chain(i, map);
}

Chain SelectorMAngle::chain(const QString &name, const PropertyMap &map) const
{
    return this->chains().chain(name, map);
}

Chain SelectorMAngle::chain(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chains().chain(chainid, map);
}

Segment SelectorMAngle::segment(int i, const PropertyMap &map) const
{
    return this->segments().segment(i, map);
}

Segment SelectorMAngle::segment(const QString &name, const PropertyMap &map) const
{
    return this->segments().segment(name, map);
}

Segment SelectorMAngle::segment(const SegID &segid, const PropertyMap &map) const
{
    return this->segments().segment(segid, map);
}

CutGroup SelectorMAngle::cutGroup(int i, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(i, map);
}

CutGroup SelectorMAngle::cutGroup(const QString &name, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(name, map);
}

CutGroup SelectorMAngle::cutGroup(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(cgid, map);
}

SelectorM<Atom> SelectorMAngle::atoms() const
{
    QList<Selector<Atom>> ret;

    for (const auto &a : this->angs)
    {
        ret.append(a.atoms());
    }

    return SelectorM<Atom>(ret);
}

SelectorM<Atom> SelectorMAngle::atoms(int i, const PropertyMap &map) const
{
    return this->atoms().atoms(i, map);
}

SelectorM<Atom> SelectorMAngle::atoms(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->atoms().atoms(slice, map);
}

SelectorM<Atom> SelectorMAngle::atoms(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->atoms().atoms(idxs, map);
}

SelectorM<Atom> SelectorMAngle::atoms(const QString &name, const PropertyMap &map) const
{
    return this->atoms().atoms(name, map);
}

SelectorM<Atom> SelectorMAngle::atoms(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atoms().atoms(atomid, map);
}

SelectorM<Residue> SelectorMAngle::residues() const
{
    QList<Selector<Residue>> ret;

    for (const auto &a : this->angs)
    {
        ret.append(a.residues());
    }

    return SelectorM<Residue>(ret);
}

SelectorM<Residue> SelectorMAngle::residues(int i, const PropertyMap &map) const
{
    return this->residues().residues(i, map);
}

SelectorM<Residue> SelectorMAngle::residues(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->residues().residues(slice, map);
}

SelectorM<Residue> SelectorMAngle::residues(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->residues().residues(idxs, map);
}

SelectorM<Residue> SelectorMAngle::residues(const QString &name, const PropertyMap &map) const
{
    return this->residues().residues(name, map);
}

SelectorM<Residue> SelectorMAngle::residues(const ResID &resid, const PropertyMap &map) const
{
    return this->residues().residues(resid, map);
}

SelectorM<Chain> SelectorMAngle::chains() const
{
    QList<Selector<Chain>> ret;

    for (const auto &a : this->angs)
    {
        ret.append(a.chains());
    }

    return SelectorM<Chain>(ret);
}

SelectorM<Chain> SelectorMAngle::chains(int i, const PropertyMap &map) const
{
    return this->chains().chains(i, map);
}

SelectorM<Chain> SelectorMAngle::chains(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->chains().chains(slice, map);
}

SelectorM<Chain> SelectorMAngle::chains(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->chains().chains(idxs, map);
}

SelectorM<Chain> SelectorMAngle::chains(const QString &name, const PropertyMap &map) const
{
    return this->chains().chains(name, map);
}

SelectorM<Chain> SelectorMAngle::chains(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chains().chains(chainid, map);
}

SelectorM<Segment> SelectorMAngle::segments() const
{
    QList<Selector<Segment>> ret;

    for (const auto &a : this->angs)
    {
        ret.append(a.segments());
    }

    return SelectorM<Segment>(ret);
}

SelectorM<Segment> SelectorMAngle::segments(int i, const PropertyMap &map) const
{
    return this->segments().segments(i, map);
}

SelectorM<Segment> SelectorMAngle::segments(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->segments().segments(slice, map);
}

SelectorM<Segment> SelectorMAngle::segments(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->segments().segments(idxs, map);
}

SelectorM<Segment> SelectorMAngle::segments(const QString &name, const PropertyMap &map) const
{
    return this->segments().segments(name, map);
}

SelectorM<Segment> SelectorMAngle::segments(const SegID &segid, const PropertyMap &map) const
{
    return this->segments().segments(segid, map);
}

SelectorM<CutGroup> SelectorMAngle::cutGroups() const
{
    QList<Selector<CutGroup>> ret;

    for (const auto &a : this->angs)
    {
        ret.append(a.cutGroups());
    }

    return SelectorM<CutGroup>(ret);
}

SelectorM<CutGroup> SelectorMAngle::cutGroups(int i, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(i, map);
}

SelectorM<CutGroup> SelectorMAngle::cutGroups(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(slice, map);
}

SelectorM<CutGroup> SelectorMAngle::cutGroups(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(idxs, map);
}

SelectorM<CutGroup> SelectorMAngle::cutGroups(const QString &name, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(name, map);
}

SelectorM<CutGroup> SelectorMAngle::cutGroups(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(cgid, map);
}

SelectResult SelectorMAngle::search(const QString &search_string) const
{
    Select search(search_string);
    return search(this->toSelectResult());
}

SelectResult SelectorMAngle::search(const QString &search_string, const PropertyMap &map) const
{
    Select search(search_string);
    return search(this->toSelectResult(), map);
}

QList<AngleID> SelectorMAngle::IDs() const
{
    QList<AngleID> ret;

    for (const auto &a : this->angs)
    {
        ret += a.IDs();
    }

    return ret;
}

int SelectorMAngle::nAtoms() const
{
    int n = 0;

    for (const auto &a : this->angs)
    {
        n += a.nAtoms();
    }

    return n;
}

int SelectorMAngle::nResidues() const
{
    int n = 0;

    for (const auto &a : this->angs)
    {
        n += a.nResidues();
    }

    return n;
}

int SelectorMAngle::nChains() const
{
    int n = 0;

    for (const auto &a : this->angs)
    {
        n += a.nChains();
    }

    return n;
}

int SelectorMAngle::nSegments() const
{
    int n = 0;

    for (const auto &a : this->angs)
    {
        n += a.nSegments();
    }

    return n;
}

int SelectorMAngle::nCutGroups() const
{
    int n = 0;

    for (const auto &a : this->angs)
    {
        n += a.nCutGroups();
    }

    return n;
}

int SelectorMAngle::nMolecules() const
{
    return this->angs.count();
}

bool SelectorMAngle::isEmpty() const
{
    return this->angs.isEmpty();
}

int SelectorMAngle::nFrames() const
{
    return this->nFrames(PropertyMap());
}

int SelectorMAngle::nFrames(const SireBase::PropertyMap &map) const
{
    return SireMol::detail::_nFrames(this->angs, map);
}

void SelectorMAngle::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

void SelectorMAngle::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

void SelectorMAngle::saveFrame()
{
    this->saveFrame(PropertyMap());
}

void SelectorMAngle::deleteAllFrames()
{
    this->deleteAllFrames(PropertyMap());
}

void SelectorMAngle::deleteAllFrames(const SireBase::PropertyMap &map)
{
    SireMol::detail::_deleteAllFrames(this->angs, map);
}

void SelectorMAngle::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

void SelectorMAngle::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_loadFrame(this->angs, frame, map);
}

void SelectorMAngle::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->angs, frame, map);
}

void SelectorMAngle::saveFrame(const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->angs, map);
}

void SelectorMAngle::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_deleteFrame(this->angs, frame, map);
}

SelectorMAngle::const_iterator SelectorMAngle::begin() const
{
    return this->angs.constBegin();
}

SelectorMAngle::const_iterator SelectorMAngle::end() const
{
    return this->angs.constEnd();
}

SelectorMAngle::const_iterator SelectorMAngle::constBegin() const
{
    return this->angs.constBegin();
}

SelectorMAngle::const_iterator SelectorMAngle::constEnd() const
{
    return this->angs.constEnd();
}

QString SelectorMAngle::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("SelectorMAngle::empty");
    }
    else
    {
        QStringList parts;

        const auto n = this->count();

        if (n <= 10)
        {
            for (int i = 0; i < n; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3").arg(i).arg(view.data().number().toString()).arg(view.toString()));
            }
        }
        else
        {
            for (int i = 0; i < 5; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3").arg(i).arg(view.data().number().toString()).arg(view.toString()));
            }

            parts.append("...");

            for (int i = n - 5; i < n; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3").arg(i).arg(view.data().number().toString()).arg(view.toString()));
            }
        }

        return QObject::tr("SelectorMAngle( size=%2\n%3\n)").arg(n).arg(parts.join("\n"));
    }
}

SelectorMAngle SelectorMAngle::add(const SelectorMAngle &other) const
{
    SelectorMAngle ret(*this);

    for (const auto &value : other)
    {
        if (ret.isEmpty())
        {
            ret.angs.append(value);
        }
        else if (ret.angs.last().isSameMolecule(value))
        {
            for (int i = 0; i < value.count(); ++i)
            {
                ret._append(value(i));
            }
        }
        else
        {
            ret.angs.append(value);
        }
    }

    return ret;
}

SelectorMAngle SelectorMAngle::intersection(const SelectorMAngle &other) const
{
    if (this->count() < other.count())
        return other.intersection(*this);

    SelectorMAngle ret;

    for (const auto &val : angs)
    {
        SelectorAngle intersect;

        for (const auto &other_val : other)
        {
            if (val.isSameMolecule(other_val))
            {
                if (intersect.isEmpty())
                    intersect = val.intersection(other_val);
                else
                    intersect = intersect.add(val.intersection(other_val));
            }
        }

        ret.angs.append(intersect);
    }

    return ret;
}

SelectorMAngle SelectorMAngle::invert(const SireBase::PropertyMap &map) const
{
    SelectorMAngle ret;

    for (const auto &val : angs)
    {
        ret.angs.append(val.invert(map));
    }

    return ret;
}

SelectorMAngle SelectorMAngle::invert() const
{
    return this->invert(PropertyMap());
}

bool SelectorMAngle::hasProperty(const SireBase::PropertyName &key) const
{
    for (const auto &val : angs)
    {
        if (val.hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorMAngle::hasMetadata(const SireBase::PropertyName &key) const
{
    for (const auto &val : angs)
    {
        if (val.hasMetadata(key))
            return true;
    }

    return false;
}

bool SelectorMAngle::hasMetadata(const SireBase::PropertyName &key, const SireBase::PropertyName &metakey) const
{
    for (const auto &val : angs)
    {
        if (val.hasMetadata(key, metakey))
            return true;
    }

    return false;
}

template <class T>
inline QSet<T> _to_set(const QList<T> &l)
{
#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
    return l.toSet();
#else
    return QSet<T>(l.constBegin(), l.constEnd());
#endif
}

QStringList SelectorMAngle::propertyKeys() const
{
    QSet<QString> keys;

    for (const auto &val : angs)
    {
        keys += _to_set(val.propertyKeys());
    }

    return keys.values();
}

QStringList SelectorMAngle::metadataKeys() const
{
    QSet<QString> keys;

    for (const auto &val : angs)
    {
        keys += _to_set(val.metadataKeys());
    }

    return keys.values();
}

QStringList SelectorMAngle::metadataKeys(const SireBase::PropertyName &key) const
{
    QSet<QString> keys;

    for (const auto &val : angs)
    {
        keys += _to_set(val.metadataKeys(key));
    }

    return keys.values();
}

QList<SireBase::Properties> SelectorMAngle::properties() const
{
    QList<SireBase::Properties> props;

    for (const auto &val : angs)
    {
        props += val.properties();
    }

    return props;
}

QList<SireBase::PropertyPtr> SelectorMAngle::property(const SireBase::PropertyName &key) const
{
    QList<SireBase::PropertyPtr> props;

    bool has_prop = false;

    for (const auto &val : angs)
    {
        try
        {
            props += val.property(key);
            has_prop = true;
        }
        catch (const SireError::exception &)
        {
            PropertyPtr null(new NullProperty());

            for (int i = 0; i < val.count(); ++i)
            {
                props.append(null);
            }
        }
    }

    if (not has_prop)
        throw SireBase::missing_property(
            QObject::tr("None of the angles in this container have a property called %1.").arg(key.source()), CODELOC);

    return props;
}

QList<SireBase::PropertyPtr> SelectorMAngle::property(const SireBase::PropertyName &key,
                                                      const Property &default_value) const
{
    QList<SireBase::PropertyPtr> props;

    for (const auto &val : angs)
    {
        props += val.property(key, default_value);
    }

    return props;
}

QList<SireUnits::Dimension::Angle> SelectorMAngle::sizes() const
{
    return this->sizes(PropertyMap());
}

QList<SireUnits::Dimension::Angle> SelectorMAngle::sizes(const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::Angle> a;

    for (const auto &val : angs)
    {
        a += val.sizes(map);
    }

    return a;
}

QList<SireUnits::Dimension::Angle> SelectorMAngle::measures() const
{
    return this->sizes();
}

QList<SireUnits::Dimension::Angle> SelectorMAngle::measures(const SireBase::PropertyMap &map) const
{
    return this->sizes(map);
}

QList<SireCAS::Expression> SelectorMAngle::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<SireCAS::Expression> SelectorMAngle::potentials(const SireBase::PropertyMap &map) const
{
    QList<SireCAS::Expression> e;

    for (const auto &val : angs)
    {
        e += val.potentials(map);
    }

    return e;
}

QList<SireUnits::Dimension::GeneralUnit> SelectorMAngle::energies() const
{
    return this->energies(PropertyMap());
}

QList<SireUnits::Dimension::GeneralUnit> SelectorMAngle::energies(const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::GeneralUnit> e;

    for (const auto &val : angs)
    {
        e += val.energies(map);
    }

    return e;
}

SireUnits::Dimension::GeneralUnit SelectorMAngle::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::GeneralUnit SelectorMAngle::energy(const SireBase::PropertyMap &map) const
{
    SireUnits::Dimension::GeneralUnit nrg(0);

    for (const auto &val : angs)
    {
        nrg += val.energy(map);
    }

    return nrg;
}
