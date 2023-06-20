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

#include "SireBase/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

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

    QStringList parts;

    if (this->nGroups() <= 10)
    {
        for (int i = 0; i < this->nGroups(); ++i)
        {
            parts.append(QString("group %1:  num_atoms = %2").arg(i).arg(matches[i].count()));
        }
    }
    else
    {
        for (int i = 0; i < 5; ++i)
        {
            parts.append(QString("group %1:  num_atoms = %2").arg(i).arg(matches[i].count()));
        }

        parts.append("...");

        for (int i = this->size() - 5; i < this->size(); ++i)
        {
            parts.append(QString("group %1:  num_atoms = %2").arg(i).arg(matches[i].count()));
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
