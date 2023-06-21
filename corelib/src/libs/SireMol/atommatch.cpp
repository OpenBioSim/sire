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
    sds << match.matches << match.reference << static_cast<const Selector<Atom> &>(match);

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, AtomMatch &match)
{
    auto v = readHeader(ds, r_atommatch);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> match.matches >> match.reference >> static_cast<Selector<Atom> &>(match);
    }
    else
        throw SireStream::version_error(v, "1", r_atommatch, CODELOC);

    return ds;
}

AtomMatch::AtomMatch() : ConcreteProperty<AtomMatch, Selector<Atom>>()
{
}

AtomMatch::AtomMatch(const MoleculeView &molview)
    : ConcreteProperty<AtomMatch, Selector<Atom>>()
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
    : ConcreteProperty<AtomMatch, Selector<Atom>>()
{
    if (not m.isEmpty())
        this->operator=(AtomMatch(molview, QList<QList<qint64>>({m})));
}

AtomMatch::AtomMatch(const Selector<Atom> &molview, const QList<QList<qint64>> &m)
    : ConcreteProperty<AtomMatch, Selector<Atom>>()
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

        Selector<Atom>::operator=(reference(atoms));
    }
}

AtomMatch::AtomMatch(const AtomMatch &other)
    : ConcreteProperty<AtomMatch, Selector<Atom>>(other),
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
        Selector<Atom>::operator=(other);
    }

    return *this;
}

bool AtomMatch::operator==(const AtomMatch &other) const
{
    return matches == other.matches and Selector<Atom>::operator==(other);
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

////////
//////// Implementation of AtomMatchM
////////

static const RegisterMetaType<AtomMatchM> r_atommatchm;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const AtomMatchM &match)
{
    writeHeader(ds, r_atommatchm, 1);

    SharedDataStream sds(ds);
    sds << match.matches << static_cast<const SelectorM<Atom> &>(match);

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, AtomMatchM &match)
{
    auto v = readHeader(ds, r_atommatchm);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> match.matches >> static_cast<SelectorM<Atom> &>(match);
    }
    else
        throw SireStream::version_error(v, "1", r_atommatchm, CODELOC);

    return ds;
}

AtomMatchM::AtomMatchM()
    : ConcreteProperty<AtomMatchM, SelectorM<Atom>>()
{
}

AtomMatchM::AtomMatchM(const AtomMatch &match)
    : ConcreteProperty<AtomMatchM, SelectorM<Atom>>()
{
    if (not match.isEmpty())
    {
        this->operator=(AtomMatchM(QList<AtomMatch>({match})));
    }
}

AtomMatchM::AtomMatchM(const QList<AtomMatch> &m)
    : ConcreteProperty<AtomMatchM, SelectorM<Atom>>()
{
    QList<Selector<Atom>> atoms;

    for (const auto &match : m)
    {
        if (not match.isEmpty())
        {
            atoms.append(match);
            matches.append(match);
        }
    }

    SelectorM<Atom>::operator=(SelectorM<Atom>(atoms));
}

AtomMatchM::AtomMatchM(const SelectResult &mols)
    : ConcreteProperty<AtomMatchM, SelectorM<Atom>>()
{
    QList<AtomMatch> atoms;

    for (const auto &mol : mols)
    {
        atoms.append(AtomMatch(mol.read()));
    }

    this->operator=(AtomMatchM(atoms));
}

AtomMatchM::AtomMatchM(const AtomMatchM &other)
    : ConcreteProperty<AtomMatchM, SelectorM<Atom>>(other),
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
        SelectorM<Atom>::operator=(other);
    }

    return *this;
}

bool AtomMatchM::operator==(const AtomMatchM &other) const
{
    return matches == other.matches and
           SelectorM<Atom>::operator==(other);
}

bool AtomMatchM::operator!=(const AtomMatchM &other) const
{
    return not this->operator==(other);
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
            const auto mols = this->molecules();
            auto atoms = match.group(i);
            atoms.update(mols[atoms.data().number()].data());
            return atoms;
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

    const auto mols = this->molecules();

    for (auto match : matches)
    {
        match.update(mols[match.data().number()].data());
        ret += match.groups();
    }

    return ret;
}
