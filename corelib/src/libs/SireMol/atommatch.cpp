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

#include "atommatch.h"

#include "SireMol/core.h"
#include "SireMol/mover_metaid.h"

#include "SireBase/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

////////
//////// Implementation of AtomMatch
////////

static const RegisterMetaType<AtomMatch> r_atommatch;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const AtomMatch &match)
{
    writeHeader(ds, r_atommatch, 1);

    SharedDataStream sds(ds);
    sds << match.matches << match.reference << static_cast<const PartialMolecule &>(match);

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, AtomMatch &match)
{
    auto v = readHeader(ds, r_atommatch);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> match.matches >> match.reference >> static_cast<PartialMolecule &>(match);
    }
    else
        throw SireStream::version_error(v, "1", r_atommatch, CODELOC);

    return ds;
}

AtomMatch::AtomMatch() : ConcreteProperty<AtomMatch, PartialMolecule>()
{
}

AtomMatch::AtomMatch(const MoleculeView &molview)
    : ConcreteProperty<AtomMatch, PartialMolecule>()
{
    if (molview.isA<AtomMatch>())
    {
        this->operator=(molview.asA<AtomMatch>());
    }
    else
    {
        auto atoms = molview.atoms();

        if (not atoms.isEmpty())
        {
            const qint64 natoms = atoms.count();

            QList<qint64> idxs;

            for (qint64 i = 0; i < natoms; ++i)
            {
                idxs.append(i);
            }

            this->operator=(AtomMatch(atoms, idxs));
        }
    }
}

AtomMatch::AtomMatch(const Selector<Atom> &molview, const QList<qint64> &m)
    : ConcreteProperty<AtomMatch, PartialMolecule>()
{
    if (not m.isEmpty())
        this->operator=(AtomMatch(molview, QList<QList<qint64>>({m})));
}

AtomMatch::AtomMatch(const Selector<Atom> &molview, const QList<QList<qint64>> &m)
    : ConcreteProperty<AtomMatch, PartialMolecule>()
{
    const int count = molview.count();

    QSet<qint64> idxs;

    for (const auto &match : m)
    {
        for (auto i : match)
        {
            if (i < -count or i >= count)
            {
                throw SireError::invalid_index(QObject::tr(
                                                   "You cannot use index '%1' if the number of atoms in the view is %2.")
                                                   .arg(i)
                                                   .arg(count),
                                               CODELOC);
            }
            else
            {
                idxs.insert(Index(i).map(count));
            }
        }

        if (not match.isEmpty())
            matches.append(match);
    }

    if (not idxs.isEmpty())
    {
        reference = molview;

        QList<qint64> atoms(idxs.constBegin(), idxs.constEnd());
        std::sort(atoms.begin(), atoms.end());

        PartialMolecule::operator=(reference(atoms));
    }
}

AtomMatch::AtomMatch(const AtomMatch &other)
    : ConcreteProperty<AtomMatch, PartialMolecule>(other),
      reference(other.reference), matches(other.matches)
{
}

AtomMatch::~AtomMatch()
{
}

const char *AtomMatch::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AtomMatch>());
}

AtomMatch &AtomMatch::operator=(const AtomMatch &other)
{
    if (this != &other)
    {
        reference = other.reference;
        matches = other.matches;
        PartialMolecule::operator=(other);
    }

    return *this;
}

bool AtomMatch::operator==(const AtomMatch &other) const
{
    return matches == other.matches and PartialMolecule::operator==(other);
}

bool AtomMatch::operator!=(const AtomMatch &other) const
{
    return not this->operator==(other);
}

QString AtomMatch::toString() const
{
    if (this->isEmpty())
        return QObject::tr("AtomMatch::empty");

    auto getLine = [](int i, const Selector<Atom> &atoms)
    {
        int natoms = atoms.nAtoms();
        bool dotdotdot = false;

        if (natoms > 5)
        {
            dotdotdot = true;
            natoms = 5;
        }

        QStringList a;

        for (int i = 0; i < natoms; ++i)
        {
            const auto atom = atoms(i);

            a.append(QString("%1:%2")
                         .arg(atom.name().value())
                         .arg(atom.number().value()));
        }

        QString atomstring = a.join(",");

        if (dotdotdot)
            atomstring += "...";

        return QString("%1: [%2] %5")
            .arg(i)
            .arg(atoms.nAtoms())
            .arg(atomstring);
    };

    const int size = this->nGroups();

    QStringList parts;

    if (size <= 10)
    {
        for (int i = 0; i < this->nGroups(); ++i)
        {
            parts.append(getLine(i, this->group(i)));
        }
    }
    else
    {
        for (int i = 0; i < 5; ++i)
        {
            parts.append(getLine(i, this->group(i)));
        }

        parts.append("...");

        for (int i = size - 5; i < size; ++i)
        {
            parts.append(getLine(i, this->group(i)));
        }
    }

    return QObject::tr("AtomMatch( size=%1\n%2\n)").arg(this->nGroups()).arg(parts.join("\n"));
}

int AtomMatch::nGroups() const
{
    return matches.count();
}

Selector<Atom> AtomMatch::group(int i) const
{
    i = Index(i).map(this->nGroups());

    auto atoms = reference(matches[i]);
    atoms.update(this->data());

    return atoms;
}

QList<Selector<Atom>> AtomMatch::groups() const
{
    QList<Selector<Atom>> g;

    for (const auto &m : matches)
    {
        auto atoms = reference(m);
        atoms.update(this->data());
        g.append(atoms);
    }

    return g;
}

AtomMatch::operator Selector<Atom>() const
{
    return this->atoms();
}

////////
//////// Implementation of AtomMatchM
////////

static const RegisterMetaType<AtomMatchM> r_atommatchm;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const AtomMatchM &match)
{
    writeHeader(ds, r_atommatchm, 1);

    SharedDataStream sds(ds);
    sds << match.matches << static_cast<const Property &>(match);

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, AtomMatchM &match)
{
    auto v = readHeader(ds, r_atommatchm);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> match.matches >> static_cast<Property &>(match);
    }
    else
        throw SireStream::version_error(v, "1", r_atommatchm, CODELOC);

    return ds;
}

AtomMatchM::AtomMatchM()
    : ConcreteProperty<AtomMatchM, Property>()
{
}

AtomMatchM::AtomMatchM(const AtomMatch &match)
    : ConcreteProperty<AtomMatchM, Property>()
{
    if (not match.isEmpty())
    {
        this->operator=(AtomMatchM(QList<AtomMatch>({match})));
    }
}

AtomMatchM::AtomMatchM(const QList<AtomMatch> &m)
    : ConcreteProperty<AtomMatchM, Property>()
{
    for (const auto &match : m)
    {
        if (not match.isEmpty())
        {
            matches.append(match);
        }
    }
}

AtomMatchM::AtomMatchM(const SelectResult &mols)
    : ConcreteProperty<AtomMatchM, Property>()
{
    QList<AtomMatch> atoms;

    for (const auto &mol : mols)
    {
        atoms.append(AtomMatch(mol.read()));
    }

    this->operator=(AtomMatchM(atoms));
}

AtomMatchM::AtomMatchM(const AtomMatchM &other)
    : ConcreteProperty<AtomMatchM, Property>(other),
      matches(other.matches)
{
}

AtomMatchM::~AtomMatchM()
{
}

const char *AtomMatchM::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AtomMatchM>());
}

AtomMatchM &AtomMatchM::operator=(const AtomMatchM &other)
{
    if (this != &other)
    {
        matches = other.matches;
        Property::operator=(other);
    }

    return *this;
}

bool AtomMatchM::operator==(const AtomMatchM &other) const
{
    return matches == other.matches and
           Property::operator==(other);
}

bool AtomMatchM::operator!=(const AtomMatchM &other) const
{
    return not this->operator==(other);
}

bool AtomMatchM::isEmpty() const
{
    return this->matches.isEmpty();
}

QString AtomMatchM::toString() const
{
    if (this->isEmpty())
        return QObject::tr("AtomMatchM::empty");

    QStringList parts;

    auto getLine = [](int i, const Selector<Atom> &atoms)
    {
        int natoms = atoms.nAtoms();
        bool dotdotdot = false;

        if (natoms > 5)
        {
            dotdotdot = true;
            natoms = 5;
        }

        QStringList a;

        for (int i = 0; i < natoms; ++i)
        {
            const auto atom = atoms(i);

            a.append(QString("%1:%2")
                         .arg(atom.name().value())
                         .arg(atom.number().value()));
        }

        QString atomstring = a.join(",");

        if (dotdotdot)
            atomstring += "...";

        return QString("%1: %3:%4 [%2] %5")
            .arg(i)
            .arg(atoms.nAtoms())
            .arg(atoms.molecule().name().value())
            .arg(atoms.molecule().number().value())
            .arg(atomstring);
    };

    const int size = this->nGroups();

    if (size <= 10)
    {
        for (int i = 0; i < this->nGroups(); ++i)
        {
            parts.append(getLine(i, this->group(i)));
        }
    }
    else
    {
        for (int i = 0; i < 5; ++i)
        {
            parts.append(getLine(i, this->group(i)));
        }

        parts.append("...");

        for (int i = size - 5; i < size; ++i)
        {
            parts.append(getLine(i, this->group(i)));
        }
    }

    return QObject::tr("AtomMatchM( size=%1\n%2\n)").arg(this->nGroups()).arg(parts.join("\n"));
}

int AtomMatchM::nGroups() const
{
    int n = 0;

    for (const auto &match : matches)
    {
        n += match.nGroups();
    }

    return n;
}

Selector<Atom> AtomMatchM::group(int i) const
{
    i = Index(i).map(this->nGroups());

    for (const auto &match : matches)
    {
        const int n = match.nGroups();

        if (i < n)
        {
            return match.group(i);
        }
        else
        {
            i -= n;
        }
    }

    throw SireError::program_bug(QObject::tr(
                                     "Should not get here! %1")
                                     .arg(i),
                                 CODELOC);

    return Selector<Atom>();
}

QList<Selector<Atom>> AtomMatchM::groups() const
{
    QList<Selector<Atom>> ret;

    for (auto match : matches)
    {
        ret += match.groups();
    }

    return ret;
}

SelectorM<Atom> AtomMatchM::_atoms() const
{
    QList<Selector<Atom>> views;

    for (const auto &match : matches)
    {
        views.append(match);
    }

    return SelectorM<Atom>(views);
}

Atom AtomMatchM::operator[](int i) const
{
    return this->_atoms().operator[](i);
}

SelectorM<Atom> AtomMatchM::operator[](const SireBase::Slice &slice) const
{
    return this->_atoms().operator[](slice);
}

SelectorM<Atom> AtomMatchM::operator[](const QList<qint64> &idxs) const
{
    return this->_atoms().operator[](idxs);
}

Atom AtomMatchM::operator[](const QString &name) const
{
    return this->_atoms().operator[](name);
}

Atom AtomMatchM::operator[](const typename Atom::ID &id) const
{
    return this->_atoms().operator[](id);
}

Atom AtomMatchM::operator()(int i) const
{
    return this->_atoms().operator()(i);
}

Atom AtomMatchM::operator()(const QString &name) const
{
    return this->_atoms().operator()(name);
}

Atom AtomMatchM::operator()(const typename Atom::ID &id) const
{
    return this->_atoms().operator()(id);
}

QList<MolViewPtr> AtomMatchM::toList() const
{
    return this->_atoms().toList();
}

Molecules AtomMatchM::toMolecules() const
{
    return this->_atoms().toMolecules();
}

int AtomMatchM::count() const
{
    return this->_atoms().count();
}

int AtomMatchM::size() const
{
    return this->_atoms().size();
}

void AtomMatchM::update(const MoleculeView &molview)
{
    this->update(molview.data());
}

void AtomMatchM::update(const MoleculeData &moldata)
{
    const auto molnum = moldata.number();

    for (auto &mol : matches)
    {
        if (mol.data().number() == molnum)
        {
            mol.update(moldata);
        }
    }
}

void AtomMatchM::update(const Molecules &molecules)
{
    for (const auto &mol : molecules)
    {
        this->update(mol.data());
    }
}

void AtomMatchM::update(const SelectorMol &molecules)
{
    for (const auto &mol : molecules)
    {
        this->update(mol.data());
    }
}

EvaluatorM AtomMatchM::evaluate() const
{
    return this->_atoms().evaluate();
}

MoleculeGroup AtomMatchM::toMoleculeGroup() const
{
    return this->_atoms().toMoleculeGroup();
}

SelectResult AtomMatchM::toSelectResult() const
{
    return this->_atoms().toSelectResult();
}

bool AtomMatchM::isSelector() const
{
    return true;
}

SelectorMol AtomMatchM::extract() const
{
    return this->_atoms().extract();
}

QList<qint64> AtomMatchM::find(const Atom &view) const
{
    return this->_atoms().find(view);
}

QList<qint64> AtomMatchM::find(const Selector<Atom> &views) const
{
    return this->_atoms().find(views);
}

QList<qint64> AtomMatchM::find(const SelectorM<Atom> &views) const
{
    return this->_atoms().find(views);
}

SelectorM<Atom> AtomMatchM::intersection(const SelectorM<Atom> &other) const
{
    return this->_atoms().intersection(other);
}

SelectorM<Atom> AtomMatchM::intersection(const Selector<Atom> &views) const
{
    return this->_atoms().intersection(views);
}

SelectorM<Atom> AtomMatchM::intersection(const Atom &view) const
{
    return this->_atoms().intersection(view);
}

SelectorM<Atom> AtomMatchM::invert() const
{
    return this->_atoms().invert();
}

bool AtomMatchM::intersects(const SelectorM<Atom> &other) const
{
    return this->_atoms().intersects(other);
}

bool AtomMatchM::intersects(const Selector<Atom> &view) const
{
    return this->_atoms().intersects(view);
}

bool AtomMatchM::intersects(const Atom &view) const
{
    return this->_atoms().intersects(view);
}

bool AtomMatchM::contains(const SelectorM<Atom> &other) const
{
    return this->_atoms().contains(other);
}

bool AtomMatchM::contains(const Selector<Atom> &view) const
{
    return this->_atoms().contains(view);
}

bool AtomMatchM::contains(const Atom &view) const
{
    return this->_atoms().contains(view);
}

Molecule AtomMatchM::molecule(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().molecule(i, map);
}

Molecule AtomMatchM::molecule(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().molecule(name, map);
}

Molecule AtomMatchM::molecule(const MolID &molid, const SireBase::PropertyMap &map)
{
    return this->_atoms().molecule(molid, map);
}

SelectorMol AtomMatchM::molecules() const
{
    return this->_atoms().molecules();
}

SelectorMol AtomMatchM::molecules(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().molecules(i, map);
}

SelectorMol AtomMatchM::molecules(const SireBase::Slice &slice,
                                  const SireBase::PropertyMap &map) const
{
    return this->_atoms().molecules(slice, map);
}

SelectorMol AtomMatchM::molecules(const QList<qint64> &idxs, const SireBase::PropertyMap &map) const
{
    return this->_atoms().molecules(idxs, map);
}

SelectorMol AtomMatchM::molecules(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().molecules(name, map);
}

SelectorMol AtomMatchM::molecules(const MolID &molid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().molecules(molid, map);
}

Atom AtomMatchM::atom(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().atom(i, map);
}

Atom AtomMatchM::atom(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().atom(name, map);
}

Atom AtomMatchM::atom(const AtomID &atomid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().atom(atomid, map);
}

Residue AtomMatchM::residue(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().residue(i, map);
}

Residue AtomMatchM::residue(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().residue(name, map);
}

Residue AtomMatchM::residue(const ResID &resid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().residue(resid, map);
}

Chain AtomMatchM::chain(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().chain(i, map);
}

Chain AtomMatchM::chain(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().chain(name, map);
}

Chain AtomMatchM::chain(const ChainID &chainid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().chain(chainid, map);
}

Segment AtomMatchM::segment(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().segment(i, map);
}

Segment AtomMatchM::segment(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().segment(name, map);
}

Segment AtomMatchM::segment(const SegID &segid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().segment(segid, map);
}

CutGroup AtomMatchM::cutGroup(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroup(i, map);
}

CutGroup AtomMatchM::cutGroup(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroup(name, map);
}

CutGroup AtomMatchM::cutGroup(const CGID &cgid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroup(cgid, map);
}

SelectorM<Atom> AtomMatchM::atoms() const
{
    return this->_atoms().atoms();
}

SelectorM<Atom> AtomMatchM::atoms(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().atoms(i, map);
}

SelectorM<Atom> AtomMatchM::atoms(const SireBase::Slice &slice,
                                  const SireBase::PropertyMap &map) const
{
    return this->_atoms().atoms(slice, map);
}

SelectorM<Atom> AtomMatchM::atoms(const QList<qint64> &idxs, const SireBase::PropertyMap &map) const
{
    return this->_atoms().atoms(idxs, map);
}

SelectorM<Atom> AtomMatchM::atoms(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().atoms(name, map);
}

SelectorM<Atom> AtomMatchM::atoms(const AtomID &atomid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().atoms(atomid, map);
}

SelectorM<Residue> AtomMatchM::residues() const
{
    return this->_atoms().residues();
}

SelectorM<Residue> AtomMatchM::residues(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().residues(i, map);
}

SelectorM<Residue> AtomMatchM::residues(const SireBase::Slice &slice,
                                        const SireBase::PropertyMap &map) const
{
    return this->_atoms().residues(slice, map);
}

SelectorM<Residue> AtomMatchM::residues(const QList<qint64> &idxs,
                                        const SireBase::PropertyMap &map) const
{
    return this->_atoms().residues(idxs, map);
}

SelectorM<Residue> AtomMatchM::residues(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().residues(name, map);
}

SelectorM<Residue> AtomMatchM::residues(const ResID &resid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().residues(resid, map);
}

SelectorM<Chain> AtomMatchM::chains() const
{
    return this->_atoms().chains();
}

SelectorM<Chain> AtomMatchM::chains(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().chains(i, map);
}

SelectorM<Chain> AtomMatchM::chains(const SireBase::Slice &slice,
                                    const SireBase::PropertyMap &map) const
{
    return this->_atoms().chains(slice, map);
}

SelectorM<Chain> AtomMatchM::chains(const QList<qint64> &idxs,
                                    const SireBase::PropertyMap &map) const
{
    return this->_atoms().chains(idxs, map);
}

SelectorM<Chain> AtomMatchM::chains(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().chains(name, map);
}

SelectorM<Chain> AtomMatchM::chains(const ChainID &chainid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().chains(chainid, map);
}

SelectorM<Segment> AtomMatchM::segments() const
{
    return this->_atoms().segments();
}

SelectorM<Segment> AtomMatchM::segments(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().segments(i, map);
}

SelectorM<Segment> AtomMatchM::segments(const SireBase::Slice &slice,
                                        const SireBase::PropertyMap &map) const
{
    return this->_atoms().segments(slice, map);
}

SelectorM<Segment> AtomMatchM::segments(const QList<qint64> &idxs,
                                        const SireBase::PropertyMap &map) const
{
    return this->_atoms().segments(idxs, map);
}

SelectorM<Segment> AtomMatchM::segments(const QString &name, const SireBase::PropertyMap &map) const
{
    return this->_atoms().segments(name, map);
}

SelectorM<Segment> AtomMatchM::segments(const SegID &segid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().segments(segid, map);
}

SelectorM<CutGroup> AtomMatchM::cutGroups() const
{
    return this->_atoms().cutGroups();
}

SelectorM<CutGroup> AtomMatchM::cutGroups(int i, const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroups(i, map);
}

SelectorM<CutGroup> AtomMatchM::cutGroups(const SireBase::Slice &slice,
                                          const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroups(slice, map);
}

SelectorM<CutGroup> AtomMatchM::cutGroups(const QList<qint64> &idxs,
                                          const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroups(idxs, map);
}

SelectorM<CutGroup> AtomMatchM::cutGroups(const QString &name,
                                          const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroups(name, map);
}

SelectorM<CutGroup> AtomMatchM::cutGroups(const CGID &cgid, const SireBase::PropertyMap &map) const
{
    return this->_atoms().cutGroups(cgid, map);
}

SelectResult AtomMatchM::search(const QString &search_string) const
{
    return this->_atoms().search(search_string);
}

SelectResult AtomMatchM::search(const QString &search_string, const SireBase::PropertyMap &map) const
{
    return this->_atoms().search(search_string, map);
}

QList<Atom::Index> AtomMatchM::IDs() const
{
    return this->_atoms().IDs();
}

QList<Atom::Index> AtomMatchM::indexes() const
{
    return this->_atoms().indexes();
}

QList<Atom::Number> AtomMatchM::numbers() const
{
    return this->_atoms().numbers();
}

QList<Atom::Name> AtomMatchM::names() const
{
    return this->_atoms().names();
}

int AtomMatchM::nAtoms() const
{
    return this->_atoms().nAtoms();
}

int AtomMatchM::nResidues() const
{
    return this->_atoms().nResidues();
}

int AtomMatchM::nChains() const
{
    return this->_atoms().nChains();
}

int AtomMatchM::nSegments() const
{
    return this->_atoms().nSegments();
}

int AtomMatchM::nCutGroups() const
{
    return this->_atoms().nCutGroups();
}

int AtomMatchM::nMolecules() const
{
    return this->_atoms().nMolecules();
}

int AtomMatchM::nFrames() const
{
    return this->_atoms().nFrames();
}

int AtomMatchM::nFrames(const SireBase::PropertyMap &map) const
{
    return this->_atoms().nFrames(map);
}

void AtomMatchM::loadFrame(int frame)
{
    for (auto &mol : matches)
    {
        mol.loadFrame(frame);
    }
}

void AtomMatchM::saveFrame(int frame)
{
    for (auto &mol : matches)
    {
        mol.saveFrame(frame);
    }
}

void AtomMatchM::saveFrame()
{
    for (auto &mol : matches)
    {
        mol.saveFrame();
    }
}

void AtomMatchM::deleteFrame(int frame)
{
    for (auto &mol : matches)
    {
        mol.deleteFrame(frame);
    }
}

void AtomMatchM::deleteAllFrames()
{
    for (auto &mol : matches)
    {
        mol.deleteAllFrames();
    }
}

void AtomMatchM::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    for (auto &mol : matches)
    {
        mol.loadFrame(frame, map);
    }
}

void AtomMatchM::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    for (auto &mol : matches)
    {
        mol.saveFrame(frame, map);
    }
}

void AtomMatchM::saveFrame(const SireBase::PropertyMap &map)
{
    for (auto &mol : matches)
    {
        mol.saveFrame(map);
    }
}

void AtomMatchM::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    for (auto &mol : matches)
    {
        mol.deleteFrame(frame, map);
    }
}

void AtomMatchM::deleteAllFrames(const SireBase::PropertyMap &map)
{
    for (auto &mol : matches)
    {
        mol.deleteAllFrames(map);
    }
}

AtomMatchM::const_iterator AtomMatchM::begin() const
{
    return matches.begin();
}

AtomMatchM::const_iterator AtomMatchM::end() const
{
    return matches.end();
}

AtomMatchM::const_iterator AtomMatchM::constBegin() const
{
    return matches.constBegin();
}

AtomMatchM::const_iterator AtomMatchM::constEnd() const
{
    return matches.constEnd();
}

bool AtomMatchM::hasProperty(const PropertyName &key) const
{
    return this->_atoms().hasProperty(key);
}

bool AtomMatchM::hasMetadata(const PropertyName &metakey) const
{
    return this->_atoms().hasMetadata(metakey);
}

bool AtomMatchM::hasMetadata(const PropertyName &key, const PropertyName &metakey) const
{
    return this->_atoms().hasMetadata(key, metakey);
}

QStringList AtomMatchM::propertyKeys() const
{
    return this->_atoms().propertyKeys();
}

QStringList AtomMatchM::metadataKeys() const
{
    return this->_atoms().metadataKeys();
}

QStringList AtomMatchM::metadataKeys(const PropertyName &key) const
{
    return this->_atoms().metadataKeys(key);
}

AtomMatchM::operator SelectorM<Atom>() const
{
    return this->_atoms();
}
