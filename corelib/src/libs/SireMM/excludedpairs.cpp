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

#include "excludedpairs.h"

#include "SireMM/cljnbpairs.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireID;
using namespace SireError;
using namespace SireStream;

static const RegisterMetaType<ExcludedPairs> r_pairs;

QDataStream &operator<<(QDataStream &ds, const ExcludedPairs &pairs)
{
    writeHeader(ds, r_pairs, 1);

    SharedDataStream sds(ds);

    sds << pairs.minfo << pairs.excl_pairs
        << static_cast<const MolViewProperty &>(pairs);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, ExcludedPairs &pairs)
{
    VersionID v = readHeader(ds, r_pairs);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pairs.minfo >> pairs.excl_pairs >> static_cast<MolViewProperty &>(pairs);
    }
    else
        throw SireStream::version_error(v, "1", r_pairs, CODELOC);

    return ds;
}

ExcludedPairs::ExcludedPairs()
    : ConcreteProperty<ExcludedPairs, MolViewProperty>()
{
}

ExcludedPairs::ExcludedPairs(const MoleculeView &molecule, const PropertyMap &map)
    : ConcreteProperty<ExcludedPairs, MolViewProperty>()
{
    minfo = MoleculeInfo(molecule.data().info());

    // try to autogenerate the pairs from the CLJNBPairs and
    // connectivity properties
    const auto &cljscl = molecule.data().property(map["intrascale"]).asA<CLJNBPairs>();
    const auto &connectivity = molecule.data().property(map["connectivity"]).asA<Connectivity>();

    // if the connectivity is empty then we have to assume that all pairs are excluded
    if (connectivity.nConnections() == 0 and minfo.nAtoms() > 1)
    {
        const int nats = minfo.nAtoms();

        for (int i = 0; i < nats - 1; ++i)
        {
            for (int j = i + 1; j < nats; ++j)
            {
                this->excl_pairs.append(i);
                this->excl_pairs.append(j);
            }
        }

        return;
    }

    const auto bond_matrix = connectivity.getBondMatrix(4);

    // loop over all the scale factors - this is potentially quite slow
    const int ncg = cljscl.nGroups();

    for (int igrp = 0; igrp < ncg; ++igrp)
    {
        const int nats0 = minfo.nAtoms(CGIdx(igrp));

        for (int jgrp = igrp; jgrp < ncg; ++jgrp)
        {
            const auto &pairs = cljscl.get(CGIdx(igrp), CGIdx(jgrp));

            if (pairs.isEmpty())
            {
                // all of these pairs have the default value
                auto scl = pairs.defaultValue();

                if (scl == CLJScaleFactor(0))
                {
                    // check that all of the atoms are bonded
                    // (any that aren't must be excluded)
                    if (igrp == jgrp)
                    {
                        for (int i = 0; i < nats0 - 1; ++i)
                        {
                            qint64 iatm = minfo.atomIdx(CGAtomIdx(CGIdx(igrp), Index(i))).value();

                            const auto &row = bond_matrix.constData()[iatm];

                            for (int j = i + 1; j < nats0; ++j)
                            {
                                qint64 jatm = minfo.atomIdx(CGAtomIdx(CGIdx(igrp), Index(j))).value();

                                if (not row.constData()[jatm])
                                {
                                    // these atoms aren't bonded - must be excluded
                                    this->excl_pairs.append(iatm);
                                    this->excl_pairs.append(jatm);
                                }
                            }
                        }
                    }
                    else
                    {
                        const int nats1 = minfo.nAtoms(CGIdx(jgrp));

                        for (int i = 0; i < nats0; ++i)
                        {
                            qint64 iatm = minfo.atomIdx(CGAtomIdx(CGIdx(igrp), Index(i))).value();

                            const auto &row = bond_matrix.constData()[iatm];

                            for (int j = 0; j < nats1; ++j)
                            {
                                qint64 jatm = minfo.atomIdx(CGAtomIdx(CGIdx(jgrp), Index(j))).value();

                                if (not row.constData()[jatm])
                                {
                                    // these atoms aren't bonded - must be excluded
                                    this->excl_pairs.append(iatm);
                                    this->excl_pairs.append(jatm);
                                }
                            }
                        }
                    }
                }
                else if (scl != CLJScaleFactor(1))
                {
                    // check that all of the atoms are 1-4 bonded...
                    qDebug() << "WARNING: INVALID 1-4 SCALING FACTOR FOR NON-BONDED GROUPS!"
                             << igrp << jgrp << scl.toString();
                }
                // else scl == CLJScaleFactor(1), meaning they are all bonded
                // and none of them are excluded
            }
            else
            {
                // no default value - need to compare each value one by one
                if (igrp == jgrp)
                {
                    for (int i = 0; i < nats0 - 1; ++i)
                    {
                        qint64 iatm = minfo.atomIdx(CGAtomIdx(CGIdx(igrp), Index(i))).value();

                        const auto &row = bond_matrix.constData()[iatm];

                        for (int j = i + 1; j < nats0; ++j)
                        {
                            const auto &scl = pairs.get(i, j);

                            qint64 jatm = minfo.atomIdx(CGAtomIdx(CGIdx(igrp), Index(j))).value();

                            if (scl == CLJScaleFactor(0))
                            {
                                if (not row.constData()[jatm])
                                {
                                    // these atoms aren't bonded - must be excluded
                                    this->excl_pairs.append(iatm);
                                    this->excl_pairs.append(jatm);
                                }
                            }
                            else if (scl != CLJScaleFactor(1))
                            {
                                if (not row.constData()[jatm])
                                {
                                    qDebug() << "WARNING: INVALID 1-4 SCALING FACTOR FOR NON-BONDED ATOMS!"
                                             << minfo.name(AtomIdx(iatm)).value()
                                             << minfo.name(AtomIdx(jatm)).value() << scl.toString();
                                }
                            }
                        }
                    }
                }
                else
                {
                    const int nats1 = minfo.nAtoms(CGIdx(jgrp));

                    for (int i = 0; i < nats0; ++i)
                    {
                        qint64 iatm = minfo.atomIdx(CGAtomIdx(CGIdx(igrp), Index(i))).value();

                        const auto &row = bond_matrix.constData()[iatm];

                        for (int j = 0; j < nats1; ++j)
                        {
                            const auto &scl = pairs.get(i, j);

                            qint64 jatm = minfo.atomIdx(CGAtomIdx(CGIdx(jgrp), Index(j))).value();

                            if (scl == CLJScaleFactor(0))
                            {
                                if (not row.constData()[jatm])
                                {
                                    // these atoms aren't bonded - must be excluded
                                    this->excl_pairs.append(iatm);
                                    this->excl_pairs.append(jatm);
                                }
                            }
                            else if (scl != CLJScaleFactor(1))
                            {
                                if (not row.constData()[jatm])
                                {
                                    qDebug() << "WARNING: INVALID 1-4 SCALING FACTOR FOR NON-BONDED ATOMS!"
                                             << minfo.name(AtomIdx(iatm)).value()
                                             << minfo.name(AtomIdx(jatm)).value() << scl.toString();
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

ExcludedPairs::ExcludedPairs(const ExcludedPairs &other)
    : ConcreteProperty<ExcludedPairs, MolViewProperty>(other),
      minfo(other.minfo), excl_pairs(other.excl_pairs)
{
}

ExcludedPairs::~ExcludedPairs()
{
}

const char *ExcludedPairs::typeName()
{
    return QMetaType::typeName(qMetaTypeId<ExcludedPairs>());
}

const char *ExcludedPairs::what() const
{
    return ExcludedPairs::typeName();
}

ExcludedPairs &ExcludedPairs::operator=(const ExcludedPairs &other)
{
    if (this != &other)
    {
        minfo = other.minfo;
        excl_pairs = other.excl_pairs;
    }

    return *this;
}

bool ExcludedPairs::operator==(const ExcludedPairs &other) const
{
    return this->minfo == other.minfo and this->excl_pairs == other.excl_pairs;
}

bool ExcludedPairs::operator!=(const ExcludedPairs &other) const
{
    return not this->operator==(other);
}

ExcludedPairs *ExcludedPairs::clone() const
{
    return new ExcludedPairs(*this);
}

int ExcludedPairs::nExcludedPairs() const
{
    return this->excl_pairs.count() / 2;
}

int ExcludedPairs::count() const
{
    return this->nExcludedPairs();
}

std::tuple<AtomIdx, AtomIdx> ExcludedPairs::operator[](int i) const
{
    i = Index(i).map(this->count());

    return std::make_tuple(AtomIdx(this->excl_pairs[2 * i]),
                           AtomIdx(this->excl_pairs[(2 * i) + 1]));
}

QString ExcludedPairs::toString() const
{
    if (this->excl_pairs.isEmpty())
        return QObject::tr("ExcludedPairs::empty");
    else
        return QObject::tr("ExcludedPair( nExcludedPairs = %1 )")
            .arg(this->nExcludedPairs());
}

MoleculeInfo ExcludedPairs::info() const
{
    return this->minfo;
}

bool ExcludedPairs::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return this->minfo == MoleculeInfo(molinfo);
}

void ExcludedPairs::updateBondMatrix(QVector<QVector<bool>> &bond_matrix) const
{
    if (this->excl_pairs.isEmpty())
        return;

    const int nats = this->minfo.nAtoms();

    if (bond_matrix.count() != nats)
    {
        throw SireError::incompatible_error(QObject::tr(
                                                "Incompatible bond matrix? %1 versus %2.")
                                                .arg(bond_matrix.count())
                                                .arg(nats),
                                            CODELOC);
    }

    auto *row_data = bond_matrix.data();

    const auto *pair_data = this->excl_pairs.constData();

    for (int i = 0; i < this->excl_pairs.count(); i += 2)
    {
        auto atom0 = pair_data[i];
        auto atom1 = pair_data[i + 1];

        auto &row0 = row_data[atom0];
        auto &row1 = row_data[atom1];

        if (row0.count() != nats or row1.count() != nats)
        {
            throw SireError::program_bug(QObject::tr(
                                             "Non-square bond matrix? %1 %2 %3")
                                             .arg(nats)
                                             .arg(row0.count())
                                             .arg(row1.count()),
                                         CODELOC);
        }

        row0.data()[atom1] = true;
        row1.data()[atom0] = true;
    }
}

int ExcludedPairs::getIndex(qint64 atom0, qint64 atom1) const
{
    for (int i = 0; i < this->excl_pairs.count(); i += 2)
    {
        qint64 p0 = this->excl_pairs[i];
        qint64 p1 = this->excl_pairs[i + 1];

        if ((atom0 == p0 and atom1 == p1) or (atom0 == p1 and atom1 == p0))
            return i;
    }

    return -1;
}

bool ExcludedPairs::areExcluded(const AtomID &atom0, const AtomID &atom1) const
{
    qint64 idx0 = this->minfo.atomIdx(atom0).value();
    qint64 idx1 = this->minfo.atomIdx(atom1).value();

    return this->getIndex(idx0, idx1) != -1;
}

void ExcludedPairs::setExcluded(const AtomID &atom0, const AtomID &atom1,
                                bool are_excluded)
{
    qint64 idx0 = this->minfo.atomIdx(atom0).value();
    qint64 idx1 = this->minfo.atomIdx(atom1).value();

    int idx = this->getIndex(idx0, idx1);

    if (are_excluded)
    {
        if (idx == -1)
        {
            this->excl_pairs.append(idx0);
            this->excl_pairs.append(idx1);
        }
    }
    else if (idx != -1)
    {
        // remove this pair
        this->excl_pairs.removeAt(idx);
        this->excl_pairs.removeAt(idx);
    }
}
