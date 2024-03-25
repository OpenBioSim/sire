/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "cljnbpairs.h"

#include "SireMol/moleculeinfo.h"

#include "SireBase/parallel.h"

#include "SireStream/datastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireID;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<CoulombNBPairs> r_coulnbpairs;
static const RegisterMetaType<LJNBPairs> r_ljnbpairs;
static const RegisterMetaType<CLJNBPairs> r_cljnbpairs;

////////
//////// Fully instantiate the template class
////////

namespace SireMM
{
    template class AtomPairs<CoulombScaleFactor>;
    template class CGAtomPairs<CoulombScaleFactor>;

    template class AtomPairs<LJScaleFactor>;
    template class CGAtomPairs<LJScaleFactor>;

    template class AtomPairs<CLJScaleFactor>;
    template class CGAtomPairs<CLJScaleFactor>;
} // namespace SireMM

////////
//////// Implementation of CoulombScaleFactor
////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CoulombScaleFactor &sclfac)
{
    ds << sclfac.cscl;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CoulombScaleFactor &sclfac)
{
    ds >> sclfac.cscl;

    return ds;
}

/** Construct with the Coulomb scale factor equal to 'scl' */
CoulombScaleFactor::CoulombScaleFactor(double scl) : cscl(scl)
{
}

/** Copy constructor */
CoulombScaleFactor::CoulombScaleFactor(const CoulombScaleFactor &other) : cscl(other.cscl)
{
}

/** Destructor */
CoulombScaleFactor::~CoulombScaleFactor()
{
}

/** Copy assignment operator */
CoulombScaleFactor &CoulombScaleFactor::operator=(const CoulombScaleFactor &other)
{
    cscl = other.cscl;

    return *this;
}

/** Comparison operator */
bool CoulombScaleFactor::operator==(const CoulombScaleFactor &other) const
{
    return cscl == other.cscl;
}

/** Comparison operator */
bool CoulombScaleFactor::operator!=(const CoulombScaleFactor &other) const
{
    return cscl != other.cscl;
}

/** Return the Coulomb parameter scaling factor */
double CoulombScaleFactor::coulomb() const
{
    return cscl;
}

////////
//////// Implementation of LJScaleFactor
////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const LJScaleFactor &sclfac)
{
    ds << sclfac.ljscl;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, LJScaleFactor &sclfac)
{
    ds >> sclfac.ljscl;

    return ds;
}

/** Construct with the LJ scale factor equal to 'scl' */
LJScaleFactor::LJScaleFactor(double scl) : ljscl(scl)
{
}

/** Copy constructor */
LJScaleFactor::LJScaleFactor(const LJScaleFactor &other) : ljscl(other.ljscl)
{
}

/** Destructor */
LJScaleFactor::~LJScaleFactor()
{
}

/** Copy assignment operator */
LJScaleFactor &LJScaleFactor::operator=(const LJScaleFactor &other)
{
    ljscl = other.ljscl;

    return *this;
}

/** Comparison operator */
bool LJScaleFactor::operator==(const LJScaleFactor &other) const
{
    return ljscl == other.ljscl;
}

/** Comparison operator */
bool LJScaleFactor::operator!=(const LJScaleFactor &other) const
{
    return ljscl != other.ljscl;
}

/** Return the LJ parameter scaling factor */
double LJScaleFactor::lj() const
{
    return ljscl;
}

////////
//////// Implementation of CLJScaleFactor
////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CLJScaleFactor &sclfac)
{
    ds << static_cast<const CoulombScaleFactor &>(sclfac) << static_cast<const LJScaleFactor &>(sclfac);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CLJScaleFactor &sclfac)
{
    ds >> static_cast<CoulombScaleFactor &>(sclfac) >> static_cast<LJScaleFactor &>(sclfac);

    return ds;
}

/** Construct with both the Coulomb and LJ scale factors equal to 'scl' */
CLJScaleFactor::CLJScaleFactor(double scl) : CoulombScaleFactor(scl), LJScaleFactor(scl)
{
}

/** Construct with 'scale_coul' Coulomb scaling, and 'scale_lj'
    LJ scaling. */
CLJScaleFactor::CLJScaleFactor(double scale_coul, double scale_lj)
    : CoulombScaleFactor(scale_coul), LJScaleFactor(scale_lj)
{
}

/** Copy constructor */
CLJScaleFactor::CLJScaleFactor(const CLJScaleFactor &other) : CoulombScaleFactor(other), LJScaleFactor(other)
{
}

/** Destructor */
CLJScaleFactor::~CLJScaleFactor()
{
}

/** Copy assignment operator */
CLJScaleFactor &CLJScaleFactor::operator=(const CLJScaleFactor &other)
{
    CoulombScaleFactor::operator=(other);
    LJScaleFactor::operator=(other);

    return *this;
}

/** Comparison operator */
bool CLJScaleFactor::operator==(const CLJScaleFactor &other) const
{
    return CoulombScaleFactor::operator==(other) and LJScaleFactor::operator==(other);
}

/** Comparison operator */
bool CLJScaleFactor::operator!=(const CLJScaleFactor &other) const
{
    return CoulombScaleFactor::operator!=(other) or LJScaleFactor::operator!=(other);
}

QString CLJScaleFactor::toString() const
{
    return QObject::tr("CLJScaleFactor( coulomb() == %1, lj() == %2 )").arg(coulomb()).arg(lj());
}

////////
//////// Implementation of CoulombNBPairs
////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CoulombNBPairs &coulnbpairs)
{
    writeHeader(ds, r_coulnbpairs, 1) << static_cast<const AtomPairs<CoulombScaleFactor> &>(coulnbpairs);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CoulombNBPairs &coulnbpairs)
{
    VersionID v = readHeader(ds, r_coulnbpairs);

    if (v == 1)
    {
        ds >> static_cast<AtomPairs<CoulombScaleFactor> &>(coulnbpairs);
    }
    else
        throw version_error(v, "1", r_coulnbpairs, CODELOC);

    return ds;
}

/** Null constructor */
CoulombNBPairs::CoulombNBPairs()
    : ConcreteProperty<CoulombNBPairs, AtomPairs<CoulombScaleFactor>>(CoulombScaleFactor(1))
{
}

/** Construct, using 'default_scale' for all of the atom-atom
    interactions in the molecule 'molinfo' */
CoulombNBPairs::CoulombNBPairs(const MoleculeInfoData &molinfo, const CoulombScaleFactor &default_scale)
    : ConcreteProperty<CoulombNBPairs, AtomPairs<CoulombScaleFactor>>(molinfo, default_scale)
{
}

/** Construct for the molecule viewed in 'molview' */
CoulombNBPairs::CoulombNBPairs(const MoleculeView &molview, const CoulombScaleFactor &default_scale)
    : ConcreteProperty<CoulombNBPairs, AtomPairs<CoulombScaleFactor>>(molview, default_scale)
{
}

/** Construct from the coulomb scaling factors in 'cljpairs' */
CoulombNBPairs::CoulombNBPairs(const CLJNBPairs &cljpairs)
    : ConcreteProperty<CoulombNBPairs, AtomPairs<CoulombScaleFactor>>(
          static_cast<const AtomPairs<CLJScaleFactor> &>(cljpairs))
{
}

/** Copy constructor */
CoulombNBPairs::CoulombNBPairs(const CoulombNBPairs &other)
    : ConcreteProperty<CoulombNBPairs, AtomPairs<CoulombScaleFactor>>(other)
{
}

/** Destructor */
CoulombNBPairs::~CoulombNBPairs()
{
}

/** Copy assignment operator */
CoulombNBPairs &CoulombNBPairs::operator=(const CoulombNBPairs &other)
{
    AtomPairs<CoulombScaleFactor>::operator=(other);
    return *this;
}

/** Copy from a CLJNBPairs object */
CoulombNBPairs &CoulombNBPairs::operator=(const CLJNBPairs &cljpairs)
{
    return this->operator=(CoulombNBPairs(cljpairs));
}

/** Comparison operator */
bool CoulombNBPairs::operator==(const CoulombNBPairs &other) const
{
    return AtomPairs<CoulombScaleFactor>::operator==(other);
}

/** Comparison operator */
bool CoulombNBPairs::operator!=(const CoulombNBPairs &other) const
{
    return AtomPairs<CoulombScaleFactor>::operator!=(other);
}

////////
//////// Implementation of LJNBPairs
////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const LJNBPairs &ljnbpairs)
{
    writeHeader(ds, r_ljnbpairs, 1) << static_cast<const AtomPairs<LJScaleFactor> &>(ljnbpairs);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, LJNBPairs &ljnbpairs)
{
    VersionID v = readHeader(ds, r_ljnbpairs);

    if (v == 1)
    {
        ds >> static_cast<AtomPairs<LJScaleFactor> &>(ljnbpairs);
    }
    else
        throw version_error(v, "1", r_ljnbpairs, CODELOC);

    return ds;
}

/** Null constructor */
LJNBPairs::LJNBPairs() : ConcreteProperty<LJNBPairs, AtomPairs<LJScaleFactor>>(LJScaleFactor(1))
{
}

/** Construct, using 'default_scale' for all of the atom-atom
    interactions in the molecule 'molinfo' */
LJNBPairs::LJNBPairs(const MoleculeInfoData &molinfo, const LJScaleFactor &default_scale)
    : ConcreteProperty<LJNBPairs, AtomPairs<LJScaleFactor>>(molinfo, default_scale)
{
}

/** Construct for the molecule viewed in 'molview' */
LJNBPairs::LJNBPairs(const MoleculeView &molview, const LJScaleFactor &default_scale)
    : ConcreteProperty<LJNBPairs, AtomPairs<LJScaleFactor>>(molview, default_scale)
{
}

/** Construct from the LJ scaling factors in 'cljpairs' */
LJNBPairs::LJNBPairs(const CLJNBPairs &cljpairs)
    : ConcreteProperty<LJNBPairs, AtomPairs<LJScaleFactor>>(static_cast<const AtomPairs<CLJScaleFactor> &>(cljpairs))
{
}

/** Copy constructor */
LJNBPairs::LJNBPairs(const LJNBPairs &other) : ConcreteProperty<LJNBPairs, AtomPairs<LJScaleFactor>>(other)
{
}

/** Destructor */
LJNBPairs::~LJNBPairs()
{
}

/** Copy assignment operator */
LJNBPairs &LJNBPairs::operator=(const LJNBPairs &other)
{
    AtomPairs<LJScaleFactor>::operator=(other);
    return *this;
}

/** Copy from a LJNBPairs object */
LJNBPairs &LJNBPairs::operator=(const CLJNBPairs &cljpairs)
{
    return this->operator=(LJNBPairs(cljpairs));
}

/** Comparison operator */
bool LJNBPairs::operator==(const LJNBPairs &other) const
{
    return AtomPairs<LJScaleFactor>::operator==(other);
}

/** Comparison operator */
bool LJNBPairs::operator!=(const LJNBPairs &other) const
{
    return AtomPairs<LJScaleFactor>::operator!=(other);
}

////////
//////// Implementation of CLJNBPairs
////////

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CLJNBPairs &cljnbpairs)
{
    writeHeader(ds, r_cljnbpairs, 1) << static_cast<const AtomPairs<CLJScaleFactor> &>(cljnbpairs);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CLJNBPairs &cljnbpairs)
{
    VersionID v = readHeader(ds, r_cljnbpairs);

    if (v == 1)
    {
        ds >> static_cast<AtomPairs<CLJScaleFactor> &>(cljnbpairs);
    }
    else
        throw version_error(v, "1", r_cljnbpairs, CODELOC);

    return ds;
}

/** Null constructor */
CLJNBPairs::CLJNBPairs() : ConcreteProperty<CLJNBPairs, AtomPairs<CLJScaleFactor>>(CLJScaleFactor(1, 1))
{
}

/** Construct, using 'default_scale' for all of the atom-atom
    interactions in the molecule 'molinfo' */
CLJNBPairs::CLJNBPairs(const MoleculeInfoData &molinfo, const CLJScaleFactor &default_scale)
    : ConcreteProperty<CLJNBPairs, AtomPairs<CLJScaleFactor>>(molinfo, default_scale)
{
}

/** Construct, using 'default_scale' for all of the atom-atom
    interactions in the molecule 'molinfo' */
CLJNBPairs::CLJNBPairs(const MoleculeInfo &molinfo, const CLJScaleFactor &default_scale)
    : ConcreteProperty<CLJNBPairs, AtomPairs<CLJScaleFactor>>(molinfo, default_scale)
{
}

/** Construct for the molecule viewed in 'molview' */
CLJNBPairs::CLJNBPairs(const MoleculeView &molview, const CLJScaleFactor &default_scale)
    : ConcreteProperty<CLJNBPairs, AtomPairs<CLJScaleFactor>>(molview, default_scale)
{
}

/** Construct, automatically setting the bonded pairs to 0 (bond and angled atoms),
    non-bonded pairs to 1, and 1-4 pairs to 'scale14' */
CLJNBPairs::CLJNBPairs(const Connectivity &connectivity, const CLJScaleFactor &scale14)
    : ConcreteProperty<CLJNBPairs, AtomPairs<CLJScaleFactor>>(connectivity.info(), CLJScaleFactor(1, 1))
{
    const auto molinfo = connectivity.info();

    if (molinfo.nAtoms() == 0)
        return;

    // create a list of connected CutGroups
    QList<std::tuple<CGIdx, CGIdx>> connected_cgroups;

    const int ncg = molinfo.nCutGroups();

    for (CGIdx i(0); i < ncg; ++i)
    {
        for (CGIdx j(i); j < ncg; ++j)
        {
            if (i == j or connectivity.areConnected(i, j))
                connected_cgroups.append(std::make_tuple(i, j));
        }
    }

    // now we have the set of connected pairs, loop through all pairs of atoms
    // of connected cutgroups and work out how they are bonded
    QMutex mutex;

    auto get_pairs = [&](std::tuple<CGIdx, CGIdx> cgpair)
    {
        // loop over all pairs of atoms between the two cutgroups
        const CGIdx cg0 = std::get<0>(cgpair);
        const CGIdx cg1 = std::get<1>(cgpair);

        const int nats0 = molinfo.nAtoms(cg0);
        const int nats1 = molinfo.nAtoms(cg1);

        // default is that atoms are not bonded, so scale factor is 1,1
        CGPairs pairs01(CLJScaleFactor(1, 1));
        pairs01.reserve(nats0, nats1);

        CGPairs pairs10 = pairs01;

        if (cg0 != cg1)
        {
            pairs10 = CGPairs(CLJScaleFactor(1, 1));
            pairs10.reserve(nats1, nats0);
        }

        for (int i = 0; i < nats0; ++i)
        {
            const CGAtomIdx cgidx0(cg0, Index(i));
            const AtomIdx atom0 = molinfo.atomIdx(cgidx0);

            for (int j = 0; j < nats1; ++j)
            {
                const CGAtomIdx cgidx1(cg1, Index(j));
                const AtomIdx atom1 = molinfo.atomIdx(cgidx1);

                int connection_type = connectivity.connectionType(atom0, atom1);

                if (connection_type > 0 and connection_type < 4)
                {
                    // this is either the same pair, bonded pair or angled pair
                    pairs01.set(i, j, CLJScaleFactor(0, 0));

                    if (cg0 != cg1)
                        pairs10.set(j, i, CLJScaleFactor(0, 0));
                }
                else if (connection_type == 4)
                {
                    // this is a 1-4 pair
                    pairs01.set(i, j, scale14);

                    if (cg0 != cg1)
                        pairs10.set(j, i, scale14);
                }
            }
        }

        // now update the global map for all cg/cg pairs
        QMutexLocker lkr(&mutex);
        cgpairs.set(cg0.value(), cg1.value(), pairs01);

        if (cg0 != cg1)
            cgpairs.set(cg1.value(), cg0.value(), pairs10);
    };

    bool uses_parallel = true;

    if (uses_parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, connected_cgroups.count()), [&](const tbb::blocked_range<int> &r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                get_pairs(connected_cgroups.at(i));
            } });
    }
    else
    {
        for (int i = 0; i < connected_cgroups.count(); ++i)
        {
            get_pairs(connected_cgroups.at(i));
        }
    }
}

/** Copy constructor */
CLJNBPairs::CLJNBPairs(const CLJNBPairs &other) : ConcreteProperty<CLJNBPairs, AtomPairs<CLJScaleFactor>>(other)
{
}

/** Destructor */
CLJNBPairs::~CLJNBPairs()
{
}

/** Copy assignment operator */
CLJNBPairs &CLJNBPairs::operator=(const CLJNBPairs &other)
{
    AtomPairs<CLJScaleFactor>::operator=(other);
    return *this;
}

/** Comparison operator */
bool CLJNBPairs::operator==(const CLJNBPairs &other) const
{
    return AtomPairs<CLJScaleFactor>::operator==(other);
}

/** Comparison operator */
bool CLJNBPairs::operator!=(const CLJNBPairs &other) const
{
    return AtomPairs<CLJScaleFactor>::operator!=(other);
}

QString CLJNBPairs::toString() const
{
    if (nAtoms() == 0)
        return QObject::tr("CLJNBPairs::null");

    return QObject::tr("CLJNBPairs( nAtoms() == %1, nGroups() == %2 )").arg(nAtoms()).arg(nGroups());
}

/** Return all of the excluded atoms for the atoms in the specified
 *  CutGroup, returned in a hash indexed by the AtomIdx of those
 *  atoms. This is equivalent to calling excludedAtoms individually,
 *  but is far more efficient if trying to get all of the
 *  excluded atoms in the whole molecule
 */
QHash<AtomIdx, QVector<AtomIdx>> CLJNBPairs::excludedAtoms(CGIdx cgidx) const
{
    const auto molinfo = info();

    const int nats = molinfo.nAtoms(cgidx);

    QHash<AtomIdx, QVector<AtomIdx>> all_cgpairs;

    if (nats == 0)
        return all_cgpairs;

    all_cgpairs.reserve(nats);

    if (molinfo.nCutGroups() == 1)
    {
        if (nats == 1)
        {
            all_cgpairs.insert(AtomIdx(0), QVector<AtomIdx>());
            return all_cgpairs;
        }
        else if (nats == 2)
        {
            all_cgpairs.insert(AtomIdx(0), QVector<AtomIdx>({AtomIdx(1)}));
            all_cgpairs.insert(AtomIdx(1), QVector<AtomIdx>({AtomIdx(0)}));
            return all_cgpairs;
        }
        else if (nats <= 3)
        {
            all_cgpairs.insert(AtomIdx(0), QVector<AtomIdx>({AtomIdx(1), AtomIdx(2)}));
            all_cgpairs.insert(AtomIdx(1), QVector<AtomIdx>({AtomIdx(0), AtomIdx(2)}));
            all_cgpairs.insert(AtomIdx(2), QVector<AtomIdx>({AtomIdx(0), AtomIdx(1)}));
            return all_cgpairs;
        }

        // we only need to worry about ourselves
        const auto cgpairs = this->get(CGIdx(0), CGIdx(0));

        if (cgpairs.isEmpty())
        {
            // all of the pairs have the same value (surprising!)
            const auto cljscl = cgpairs.defaultValue();

            if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
            {
                // all of the pairs are excluded for all atoms
                QVector<AtomIdx> all_excluded(nats);
                auto all_excluded_data = all_excluded.data();

                for (int i = 0; i < nats; ++i)
                {
                    all_excluded_data[i] = AtomIdx(i);
                }

                for (int i = 0; i < nats; ++i)
                {
                    all_cgpairs.insert(AtomIdx(i), all_excluded);
                    all_cgpairs[AtomIdx(i)].removeAll(AtomIdx(i));
                }
            }
        }
        else
        {
            // the pairs have different values, so add these in
            for (int i = 0; i < nats; ++i)
            {
                QVector<AtomIdx> all_excluded;

                for (int j = 0; j < nats; ++j)
                {
                    if (i != j)
                    {
                        const auto cljscl = cgpairs.get(i, j);

                        if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
                        {
                            // this pair is excluded
                            all_excluded.append(AtomIdx(j));
                        }
                    }
                }

                all_cgpairs.insert(AtomIdx(i), all_excluded);
            }
        }
    }
    else
    {
        QVector<AtomIdx> atomidxs(nats);

        for (int i = 0; i < nats; ++i)
        {
            atomidxs[i] = molinfo.atomIdx(CGAtomIdx(CGIdx(cgidx), Index(i)));
            all_cgpairs.insert(atomidxs[i], QVector<AtomIdx>());
        }

        const auto atomidxs_data = atomidxs.constData();

        for (int jcg = 0; jcg < molinfo.nCutGroups(); ++jcg)
        {
            const auto cgpairs = this->get(cgidx, CGIdx(jcg));

            if (cgpairs.isEmpty())
            {
                const auto cljscl = cgpairs.defaultValue();

                if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
                {
                    // all of the pairs are excluded for all atoms
                    const int other_nats = molinfo.nAtoms(CGIdx(jcg));
                    QVector<AtomIdx> all_excluded(other_nats);
                    auto all_excluded_data = all_excluded.data();

                    for (int j = 0; j < other_nats; ++j)
                    {
                        all_excluded_data[j] = molinfo.atomIdx(CGAtomIdx(CGIdx(jcg), Index(j)));
                    }

                    for (int i = 0; i < nats; ++i)
                    {
                        const auto &atomidx = atomidxs_data[i];

                        all_cgpairs[atomidx] += all_excluded;

                        // don't included the atom with itself
                        int idx = all_excluded.indexOf(atomidx);

                        if (idx >= 0)
                        {
                            all_cgpairs[atomidx].removeAll(atomidx);
                        }
                    }
                }
            }
            else
            {
                // the pairs have different values, so add these in
                const int other_nats = molinfo.nAtoms(CGIdx(jcg));

                for (int i = 0; i < nats; ++i)
                {
                    QVector<AtomIdx> all_excluded;

                    for (int j = 0; j < other_nats; ++j)
                    {
                        const auto cljscl = cgpairs.get(i, j);

                        if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
                        {
                            // this pair is excluded
                            all_excluded.append(molinfo.atomIdx(CGAtomIdx(CGIdx(jcg), Index(j))));
                        }
                    }

                    if (not all_excluded.isEmpty())
                    {
                        // don't include the atom excluded from itself
                        int idx = all_excluded.indexOf(atomidxs_data[i]);

                        if (idx >= 0)
                            all_excluded.remove(idx);

                        all_cgpairs[atomidxs_data[i]] += all_excluded;
                    }
                }
            }
        }
    }

    return all_cgpairs;
}

/** Return the excluded atoms for the atom matching ID 'atomid'. This
    returns all of the atoms for which the interaction with atomid is
    equal to zero */
QVector<AtomIdx> CLJNBPairs::excludedAtoms(const AtomID &atomid) const
{
    const auto molinfo = info();

    auto const cgatomidx = molinfo.cgAtomIdx(atomid);

    QVector<AtomIdx> excluded_atoms;

    // loop through all of the CGAtomPairs
    for (int i = 0; i < molinfo.nCutGroups(); ++i)
    {
        const auto cgpairs = this->get(cgatomidx.cutGroup(), CGIdx(i));

        if (cgpairs.isEmpty())
        {
            // all of the pairs have the same value
            const auto cljscl = cgpairs.defaultValue();

            if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
            {
                // all of the pairs are excluded!
                for (int j = 0; j < molinfo.nAtoms(CGIdx(i)); ++j)
                {
                    excluded_atoms.append(molinfo.atomIdx(CGAtomIdx(CGIdx(i), Index(j))));
                }
            }
        }
        else
        {
            // the pairs have different values, so add these in
            for (int j = 0; j < molinfo.nAtoms(CGIdx(i)); ++j)
            {
                const auto cljscl = cgpairs.get(cgatomidx.atom().value(), j);

                if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
                {
                    // this pair is excluded
                    excluded_atoms.append(molinfo.atomIdx(CGAtomIdx(CGIdx(i), Index(j))));
                }
            }
        }
    }

    // don't include the atom excluded from itself
    int idx = excluded_atoms.indexOf(molinfo.atomIdx(atomid));

    if (idx >= 0)
        excluded_atoms.remove(idx);

    return excluded_atoms;
}

/** Return the number of excluded atoms for the atom matching ID 'atomid'.
    This returns the number of atoms that do not interact with this atom
    using a non-bonded potential */
int CLJNBPairs::nExcludedAtoms(const AtomID &atomid) const
{
    const auto molinfo = info();

    auto const cgatomidx = molinfo.cgAtomIdx(atomid);

    int nexcluded = 0;

    // loop through all of the CGAtomPairs
    for (int i = 0; i < molinfo.nCutGroups(); ++i)
    {
        const auto cgpairs = this->get(cgatomidx.cutGroup(), CGIdx(i));

        if (cgpairs.isEmpty())
        {
            // all of the pairs have the same value
            const auto cljscl = cgpairs.defaultValue();

            if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
            {
                nexcluded += molinfo.nAtoms(CGIdx(i));
            }
        }
        else
        {
            // the pairs have different values, so add these in
            for (int j = 0; j < molinfo.nAtoms(CGIdx(i)); ++j)
            {
                const auto cljscl = cgpairs.get(cgatomidx.atom().value(), j);

                if (cljscl.coulomb() == 0 and cljscl.lj() == 0)
                {
                    // this pair is excluded
                    nexcluded += 1;
                }
            }
        }
    }

    // don't count the atom excluded with itself
    return nexcluded - 1;
}

/** Return the total number of atoms that are excluded from the internal
    non-bonded calculation. These are atoms that do not interact with any
    other atoms (e.g. because their nbscl factors to all other atoms in
    the molecule are zero) */
int CLJNBPairs::nExcludedAtoms() const
{
    // loop through all atoms and find those that don't interact with any others
    int nats = info().nAtoms();
    int nexcl = 0;

    for (AtomIdx i(0); i < nats; ++i)
    {
        if (this->nExcludedAtoms(i) == nats - 1)
        {
            nexcl += 1;
        }
    }

    return nexcl;
}

/** Return the IDs of atoms that don't interact with any other atom in
    the intramolecular non-bonded calculation (their scale factors to all
    other atoms is zero) */
QVector<AtomIdx> CLJNBPairs::excludedAtoms() const
{
    // loop through all atoms and find those that don't interact with any others
    int nats = info().nAtoms();
    QVector<AtomIdx> excl;

    for (AtomIdx i(0); i < nats; ++i)
    {
        if (this->nExcludedAtoms(i) == nats - 1)
        {
            excl.append(i);
        }
    }

    return excl;
}

const char *CoulombScaleFactor::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CoulombScaleFactor>());
}

const char *LJScaleFactor::typeName()
{
    return QMetaType::typeName(qMetaTypeId<LJScaleFactor>());
}

const char *CLJScaleFactor::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CLJScaleFactor>());
}

const char *CLJNBPairs::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CLJNBPairs>());
}

const char *CoulombNBPairs::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CoulombNBPairs>());
}

const char *LJNBPairs::typeName()
{
    return QMetaType::typeName(qMetaTypeId<LJNBPairs>());
}

/** Merge this property with another property */
PropertyList CLJNBPairs::merge(const MolViewProperty &other,
                               const AtomIdxMapping &mapping,
                               const QString &ghost,
                               const SireBase::PropertyMap &map) const
{
    if (not other.isA<CLJNBPairs>())
    {
        throw SireError::incompatible_error(QObject::tr("Cannot merge %1 with %2 as they are different types.")
                                                .arg(this->what())
                                                .arg(other.what()),
                                            CODELOC);
    }

    if (not ghost.isEmpty())
    {
        Console::warning(QObject::tr("The ghost parameter '%1' for CLJNBPairs parameters is ignored").arg(ghost));
    }

    const CLJNBPairs &ref = *this;
    const CLJNBPairs &pert = other.asA<CLJNBPairs>();

    CLJNBPairs prop0 = ref;
    CLJNBPairs prop1 = ref;

    // we now go through all of the atoms that are mapped and set the
    // CLJ NB pair to the right value for each end state. We copy the
    // values from the alternate end state for ghost atoms, as we can
    // assume that the ghost atoms will have the same bonding
    // arrangement as in their end state
    for (auto it1 = mapping.begin(); it1 != mapping.end(); ++it1)
    {
        const auto &atom_a = *it1;

        for (auto it2 = it1 + 1; it2 != mapping.end(); ++it2)
        {
            const auto &atom_b = *it2;

            if (atom_a.isUnmappedIn0() or atom_b.isUnmappedIn0())
            {
                // this pair does not exist in the reference state
                if (atom_a.isUnmappedIn1() or atom_b.isUnmappedIn1())
                {
                    // this pair does not exist in the perturbed state either.
                    // This pair should not interact with each other
                    prop0.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(), CLJScaleFactor(0, 0));
                    prop1.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(), CLJScaleFactor(0, 0));
                }
                else
                {
                    // set both end states to the value in the perturbed state
                    const auto &scl = pert.get(atom_a.cgAtomIdx1(), atom_b.cgAtomIdx1());

                    prop0.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(), scl);
                    prop1.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(), scl);
                }
            }
            else if (atom_a.isUnmappedIn1() or atom_b.isUnmappedIn1())
            {
                if (atom_a.isUnmappedIn0() or atom_b.isUnmappedIn0())
                {
                    // this pair does not exist in the reference state
                    // This pair should not interact with each other
                    prop0.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(), CLJScaleFactor(0, 0));
                    prop1.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(), CLJScaleFactor(0, 0));
                }
                else
                {
                    // set both end states to the value in the reference state
                    const auto &scl = ref.get(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0());
                    // already set in the reference state
                    prop1.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(), scl);
                }
            }
            else
            {
                // we only need to update the pertubed state to equal the
                // value from the perturbed parameters
                prop1.set(atom_a.cgAtomIdx0(), atom_b.cgAtomIdx0(),
                          pert.get(atom_a.cgAtomIdx1(), atom_b.cgAtomIdx1()));
            }
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

    This will only copy the values of pairs of atoms that are
    successfully matched - all other pairs will have the default
    value of this AtomPairs object.

    \throw SireError::incompatible_error
*/
SireBase::PropertyPtr CLJNBPairs::_pvt_makeCompatibleWith(
    const MoleculeInfoData &other_info, const AtomMatcher &atommatcher) const
{
    // if the atommatcher doesn't change order and the new molecule info
    // has the same number of atoms in the same number of cutgroups, then
    // there is nothing that we need to do
    if (not atommatcher.changesOrder(this->info(), other_info))
    {
        bool same_arrangement = true;

        // ensure that the number of atoms and number of cutgroups are the same
        if (this->info().nAtoms() == other_info.nAtoms() and this->info().nCutGroups() == other_info.nCutGroups())
        {
            for (int i = 0; i < other_info.nCutGroups(); ++i)
            {
                if (this->info().nAtoms(CGIdx(i)) != other_info.nAtoms(CGIdx(i)))
                {
                    same_arrangement = false;
                    break;
                }
            }
        }

        if (same_arrangement)
        {
            // there is no change in the atom order - this AtomPairs object is still valid,
            // create a copy of the object and update the molinfo
            CLJNBPairs ret(*this);
            ret.molinfo = other_info;
            return ret;
        }
    }

    QHash<AtomIdx, AtomIdx> matched_atoms = atommatcher.match(this->info(), other_info);

    // check to see if the AtomIdx to AtomIdx map changes - if not, then
    // we can return a copy of this object with the new molinfo
    bool same_mapping = true;

    for (auto it = matched_atoms.begin(); it != matched_atoms.end(); ++it)
    {
        if (it.key() != it.value())
        {
            same_mapping = false;
            break;
        }
    }

    if (same_mapping)
    {
        CLJNBPairs ret(*this);
        ret.molinfo = other_info;
        return ret;
    }
    else
        return this->_pvt_makeCompatibleWith(other_info, matched_atoms);
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom map in
    'map' to match atoms from the current molecule to the atoms in the
    molecule whose layout is in 'molinfo'

    This will only copy the values of pairs of atoms that are
    successfully matched - all other pairs will have the default
    value of this AtomPairs object.

    \throw SireError::incompatible_error
*/
SireBase::PropertyPtr CLJNBPairs::_pvt_makeCompatibleWith(
    const MoleculeInfoData &other_info, const QHash<AtomIdx, AtomIdx> &map) const
{
    const auto &this_info = this->info();

    // create a map from CGAtomIdx to CGAtomIdx for both states
    // Only insert values where they have changed - use a null
    // value to indicate that the atom does not exist in the new map
    QHash<CGAtomIdx, CGAtomIdx> cg_map;
    cg_map.reserve(map.count());

    QSet<CGIdx> changed_cgroups, deleted_cgroups, mapped_cgroups;
    const int ncg = this_info.nCutGroups();
    mapped_cgroups.reserve(ncg);
    changed_cgroups.reserve(ncg);
    deleted_cgroups.reserve(ncg);

    for (CGIdx i(0); i < ncg; ++i)
    {
        deleted_cgroups.insert(i);
    }

    for (auto it = map.begin(); it != map.end(); ++it)
    {
        CGAtomIdx atom0 = this_info.cgAtomIdx(it.key());
        CGAtomIdx atom1;

        if (not it.value().isNull())
            atom1 = other_info.cgAtomIdx(it.value());

        if (not atom1.isNull())
        {
            deleted_cgroups.remove(atom1.cutGroup());
            mapped_cgroups.insert(atom0.cutGroup());
        }

        if (atom0 != atom1)
        {
            // this has changed
            cg_map.insert(atom0, atom1);
            changed_cgroups.insert(atom0.cutGroup());

            if (not atom1.isNull())
                changed_cgroups.insert(atom1.cutGroup());
        }
    }

    if (cg_map.isEmpty())
    {
        // nothing has changed - we don't need to do any work
        CLJNBPairs ret(*this);
        ret.molinfo = other_info;
        return ret;
    }

    // there are some changes - start by creating a completely
    // empty set of pairs, using a default value of 1,1
    CLJNBPairs ret(other_info, CLJScaleFactor(1, 1));

    // now go through all of the atom pairs, in CGIdx order, and
    // copy where we can from this object to the new object, and
    // if not possible, then copy individual values
    for (CGIdx i(0); i < ncg; ++i)
    {
        bool changed_i = changed_cgroups.contains(i);
        const int nats_i = this_info.nAtoms(i);

        if (not mapped_cgroups.contains(i))
        {
            // this CutGroup has been deleted
            continue;
        }

        for (CGIdx j(i); j < ncg; ++j)
        {
            if (not mapped_cgroups.contains(j))
            {
                // this CutGroup has been deleted
                continue;
            }

            bool changed_j = changed_cgroups.contains(j);

            const auto &cgpairs = this->get(i, j);

            if (not(changed_i or changed_j))
            {
                // nothing has changed, so copy in the original values (only if the CutGroup
                // pair hasn't been deleted)
                if (not(deleted_cgroups.contains(i) or deleted_cgroups.contains(j)))
                    ret.cgpairs.set(i, j, cgpairs);

                continue;
            }

            // there's change, so just copy the values for all atom pairs
            const int nats_j = this_info.nAtoms(j);

            auto new_cgpairs = CGPairs(CLJScaleFactor(1, 1));
            bool changed_atom_pair = false;

            for (int atom_i = 0; atom_i < nats_i; ++atom_i)
            {
                auto new_cgidx_i = cg_map.value(CGAtomIdx(i, Index(atom_i)), CGAtomIdx(i, Index(atom_i)));

                if (new_cgidx_i.isNull())
                    // this atom isn't mapped, so don't copy any values
                    continue;

                for (int atom_j = 0; atom_j < nats_j; ++atom_j)
                {
                    auto new_cgidx_j = cg_map.value(CGAtomIdx(j, Index(atom_j)), CGAtomIdx(j, Index(atom_j)));

                    if (new_cgidx_j.isNull())
                    {
                        // this atom isn't mapped, so don't copy any values
                        continue;
                    }

                    // get the current value at the current index
                    const auto &scl0 = cgpairs.get(atom_i, atom_j);

                    // set the new value at the new index
                    if (new_cgidx_i.cutGroup() == i and new_cgidx_j.cutGroup() == j)
                    {
                        // this is in the current CutGroup pair, so can set directly
                        new_cgpairs.set(new_cgidx_i.atom().value(), new_cgidx_j.atom().value(), scl0);
                        changed_atom_pair = true;
                    }
                    else
                    {
                        // this is in a completely different CutGroup pair!
                        ret.set(new_cgidx_i, new_cgidx_j, scl0);
                    }
                }
            }

            // save the cgpairs
            if (changed_atom_pair)
            {
                ret.cgpairs.set(i, j, new_cgpairs);
            }
        }
    }

    return ret;
}
