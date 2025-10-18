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

#include <QDataStream>
#include <QElapsedTimer>

#include <boost/assert.hpp>

#include "atommatcher.h"
#include "atomselection.h"
#include "connectivity.h"
#include "moleculedata.h"
#include "moleculeinfo.h"
#include "moleculeinfodata.h"
#include "moleculeview.h"
#include "atomidxmapping.h"

#include "angleid.h"
#include "bondid.h"
#include "dihedralid.h"
#include "improperid.h"

#include "SireMol/errors.h"

#include "SireBase/errors.h"
#include "SireBase/console.h"
#include "SireBase/parallel.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDebug>

using namespace SireStream;
using namespace SireMol;
using namespace SireBase;

//////
////// Implementation of detail::IDPair
//////

SIRE_ALWAYS_INLINE QDataStream &operator<<(QDataStream &ds, const SireMol::detail::IDPair &idpair)
{
    ds << idpair.atom0 << idpair.atom1;
    return ds;
}

SIRE_ALWAYS_INLINE QDataStream &operator>>(QDataStream &ds, SireMol::detail::IDPair &idpair)
{
    ds >> idpair.atom0 >> idpair.atom1;
    return ds;
}

SireMol::detail::IDPair::IDPair(quint32 atm0, quint32 atm1) : atom0(atm0), atom1(atm1)
{
    if (atm0 > atm1)
    {
        qSwap(atom0, atom1);
    }
}

SireMol::detail::IDPair::IDPair(const SireMol::detail::IDPair &other) : atom0(other.atom0), atom1(other.atom1)
{
}

SireMol::detail::IDPair::~IDPair()
{
}

SireMol::detail::IDPair &SireMol::detail::IDPair::operator=(const SireMol::detail::IDPair &other)
{
    atom0 = other.atom0;
    atom1 = other.atom1;
    return *this;
}

bool SireMol::detail::IDPair::operator==(const SireMol::detail::IDPair &other) const
{
    return atom0 == other.atom0 and atom1 == other.atom1;
}

bool SireMol::detail::IDPair::operator!=(const SireMol::detail::IDPair &other) const
{
    return atom0 != other.atom0 or atom1 != other.atom1;
}

bool SireMol::detail::IDPair::operator<(const SireMol::detail::IDPair &other) const
{
    return atom0 < other.atom0 or ((atom0 == other.atom0) and (atom1 < other.atom1));
}

bool SireMol::detail::IDPair::operator<=(const SireMol::detail::IDPair &other) const
{
    return atom0 < other.atom0 or ((atom0 == other.atom0) and (atom1 <= other.atom1));
}

bool SireMol::detail::IDPair::operator>(const SireMol::detail::IDPair &other) const
{
    return atom0 > other.atom0 or ((atom0 == other.atom0) and (atom1 > other.atom1));
}

bool SireMol::detail::IDPair::operator>=(const SireMol::detail::IDPair &other) const
{
    return atom0 > other.atom0 or ((atom0 == other.atom0) and (atom1 >= other.atom1));
}

/////////
///////// Implementation of ConnectivityBase
/////////

static const RegisterMetaType<ConnectivityBase> r_conbase(MAGIC_ONLY, "SireMol::ConnectivityBase");

/** Serialise ConnectivityBase */
QDataStream &operator<<(QDataStream &ds, const ConnectivityBase &conbase)
{
    writeHeader(ds, r_conbase, 4);

    SharedDataStream sds(ds);

    SharedDataPointer<MoleculeInfoData> d(conbase.minfo);

    sds << conbase.connected_atoms << conbase.connected_res << conbase.bond_props << conbase.ang_props
        << conbase.dih_props << conbase.imp_props << d << static_cast<const MolViewProperty &>(conbase);

    return ds;
}

/** Deserialise a MoleculeBonds */
QDataStream &operator>>(QDataStream &ds, ConnectivityBase &conbase)
{
    VersionID v = readHeader(ds, r_conbase);

    if (v == 4)
    {
        SharedDataStream sds(ds);

        SharedDataPointer<MoleculeInfoData> d;

        sds >> conbase.connected_atoms >> conbase.connected_res >> conbase.bond_props >> conbase.ang_props >>
            conbase.dih_props >> conbase.imp_props >> d >> static_cast<MolViewProperty &>(conbase);

        conbase.minfo = d;
    }
    else if (v == 3)
    {
        SharedDataStream sds(ds);

        SharedDataPointer<MoleculeInfoData> d;

        sds >> conbase.connected_atoms >> conbase.connected_res >> conbase.bond_props >> d >>
            static_cast<MolViewProperty &>(conbase);

        conbase.ang_props.clear();
        conbase.dih_props.clear();
        conbase.imp_props.clear();

        conbase.minfo = d;
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);

        SharedDataPointer<MoleculeInfoData> d;

        conbase.bond_props.clear();
        conbase.ang_props.clear();
        conbase.dih_props.clear();
        conbase.imp_props.clear();

        sds >> conbase.connected_atoms >> conbase.connected_res >> d >> static_cast<MolViewProperty &>(conbase);

        conbase.minfo = d;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        SharedDataPointer<MoleculeInfoData> d;

        conbase.bond_props.clear();
        conbase.ang_props.clear();
        conbase.dih_props.clear();
        conbase.imp_props.clear();

        sds >> conbase.connected_atoms >> conbase.connected_res >> d >> static_cast<Property &>(conbase);

        conbase.minfo = d;
    }
    else
        throw version_error(v, "1,2", r_conbase, CODELOC);

    return ds;
}

/** Null constructor */
ConnectivityBase::ConnectivityBase() : MolViewProperty()
{
}

/** Return the info object that describes the molecule for which this connectivity applies */
MoleculeInfo ConnectivityBase::info() const
{
    return minfo;
}

/** Construct the connectivity for molecule described by
    the passed info object */
ConnectivityBase::ConnectivityBase(const MoleculeInfo &info) : MolViewProperty(), minfo(info)
{
    if (minfo.nAtoms() > 0)
    {
        connected_atoms.resize(minfo.nAtoms());
        connected_atoms.squeeze();
    }

    if (minfo.nResidues() > 0)
    {
        connected_res.resize(minfo.nResidues());
        connected_res.squeeze();
    }
}

/** Construct the connectivity for molecule described by
    the passed info object */
ConnectivityBase::ConnectivityBase(const MoleculeData &moldata) : MolViewProperty(), minfo(moldata.info())
{
    if (info().nAtoms() > 0)
    {
        connected_atoms.resize(info().nAtoms());
        connected_atoms.squeeze();
    }

    if (info().nResidues() > 0)
    {
        connected_res.resize(info().nResidues());
        connected_res.squeeze();
    }
}

/** Copy constructor */
ConnectivityBase::ConnectivityBase(const ConnectivityBase &other)
    : MolViewProperty(other), connected_atoms(other.connected_atoms), connected_res(other.connected_res),
      bond_props(other.bond_props), ang_props(other.ang_props), dih_props(other.dih_props), imp_props(other.imp_props),
      minfo(other.minfo)
{
}

/** Destructor */
ConnectivityBase::~ConnectivityBase()
{
}

/** Copy assignment operator */
ConnectivityBase &ConnectivityBase::operator=(const ConnectivityBase &other)
{
    if (this != &other)
    {
        connected_atoms = other.connected_atoms;
        connected_res = other.connected_res;
        bond_props = other.bond_props;
        ang_props = other.ang_props;
        dih_props = other.dih_props;
        imp_props = other.imp_props;
        minfo = other.minfo;
    }

    return *this;
}

/** Comparison operator */
bool ConnectivityBase::operator==(const ConnectivityBase &other) const
{
    return minfo == other.minfo and connected_atoms == other.connected_atoms and bond_props == other.bond_props and
           ang_props == other.ang_props and dih_props == other.dih_props and imp_props == other.imp_props;
}

/** Comparison operator */
bool ConnectivityBase::operator!=(const ConnectivityBase &other) const
{
    return not operator==(other);
}

bool ConnectivityBase::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return minfo == MoleculeInfo(molinfo);
}

PropertyPtr ConnectivityBase::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                      const AtomMatcher &atommatcher) const
{
    try
    {
        if (not atommatcher.changesOrder(this->info(), molinfo))
        {
            // the order of the atoms remains the same - this means that the
            // AtomIdx indicies are still valid
            Connectivity ret;
            ret.connected_atoms = connected_atoms;
            ret.connected_atoms.resize(molinfo.nAtoms());
            ret.connected_res = connected_res;
            ret.connected_res.resize(molinfo.nResidues());
            ret.bond_props = bond_props;
            ret.minfo = MoleculeInfo(molinfo);
            return ret;
        }

        QHash<AtomIdx, AtomIdx> matched_atoms = atommatcher.match(this->info(), molinfo);

        // see if there is any actual change - if not, then we can just return
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
            Connectivity ret;
            ret.connected_atoms = connected_atoms;
            ret.connected_atoms.resize(molinfo.nAtoms());
            ret.connected_res = connected_res;
            ret.connected_res.resize(molinfo.nResidues());
            ret.bond_props = bond_props;
            ret.minfo = MoleculeInfo(molinfo);

            if (ret.minfo.nAtoms() < this->minfo.nAtoms())
            {
                const int nats = ret.minfo.nAtoms();
                const int nres = ret.minfo.nResidues();

                // atoms have been removed, so need to remove any references
                // to atoms that no longer exist
                for (int i = 0; i < ret.connected_atoms.count(); ++i)
                {
                    QSet<AtomIdx> &connected = ret.connected_atoms[i];

                    for (auto it = connected.begin(); it != connected.end();)
                    {
                        if (it->value() >= nats)
                        {
                            it = connected.erase(it);
                        }
                        else
                        {
                            ++it;
                        }
                    }
                }

                for (int i = 0; i < ret.connected_res.count(); ++i)
                {
                    QSet<ResIdx> &connected = ret.connected_res[i];

                    for (auto it = connected.begin(); it != connected.end();)
                    {
                        if (it->value() >= nres)
                        {
                            it = connected.erase(it);
                        }
                        else
                        {
                            ++it;
                        }
                    }
                }

                // remove any bond properties that refer to atoms that no longer exist
                for (auto it = ret.bond_props.begin(); it != ret.bond_props.end();)
                {
                    if (it.key().atom0 >= nats or it.key().atom1 >= nats)
                    {
                        it = ret.bond_props.erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }

                // clear the angle, dihedral and improper properties as
                // these are too complex to make compatible
                ret.ang_props.clear();
                ret.dih_props.clear();
                ret.imp_props.clear();
            }

            return ret;
        }
        else
            return this->_pvt_makeCompatibleWith(molinfo, matched_atoms);
    }
    catch (const SireError::exception &e)
    {
        throw;
        return Connectivity();
    }
}

PropertyPtr ConnectivityBase::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                      const QHash<AtomIdx, AtomIdx> &map) const
{
    try
    {
        ConnectivityEditor editor;
        editor.minfo = MoleculeInfo(molinfo);
        editor.connected_atoms = QVector<QSet<AtomIdx>>(molinfo.nAtoms());
        editor.connected_res = QVector<QSet<ResIdx>>(molinfo.nResidues());

        for (int i = 0; i < connected_atoms.count(); ++i)
        {
            AtomIdx old_idx(i);

            AtomIdx new_idx = map.value(old_idx, AtomIdx(-1));

            if (new_idx != -1)
            {
                foreach (AtomIdx old_bond, this->connectionsTo(old_idx))
                {
                    AtomIdx new_bond = map.value(old_bond, AtomIdx(-1));

                    if (new_bond != -1)
                    {
                        if (new_bond > new_idx)
                        {
                            editor.connect(new_idx, new_bond);
                            Properties props = bond_props.value(SireMol::detail::IDPair(old_idx, old_bond));

                            if (not props.isEmpty())
                            {
                                editor.bond_props[SireMol::detail::IDPair(new_idx, new_bond)] = props;
                            }
                        }
                    }
                }
            }
        }

        // clear the angle, dihedral and improper properties as
        // these are too complex to make compatible
        editor.ang_props.clear();
        editor.dih_props.clear();
        editor.imp_props.clear();

        return editor.commit();
    }
    catch (const SireError::exception &e)
    {
        throw;
        return Connectivity();
    }
}

static QString atomString(const MoleculeInfoData &molinfo, AtomIdx atom)
{
    if (molinfo.isWithinResidue(atom))
        return QString("%1:%2:%3")
            .arg(molinfo.name(atom))
            .arg(molinfo.name(molinfo.parentResidue(atom)))
            .arg(molinfo.number(molinfo.parentResidue(atom)));
    else
        return QString("%1:%2").arg(molinfo.name(atom)).arg(atom);
}

QString ConnectivityBase::toString() const
{
    QStringList lines;

    lines.append(QObject::tr("Connectivity: nConnections() == %1.").arg(this->nConnections()));

    if (nConnections() == 0)
        return lines.at(0);

    if (not connected_res.isEmpty())
    {
        lines.append(QObject::tr("Connected residues:"));

        for (int i = 0; i < std::min(5, connected_res.count()); ++i)
        {
            QStringList resnums;

            foreach (ResIdx j, connected_res.at(i))
            {
                resnums.append(QString("%1:%2").arg(minfo.name(j)).arg(minfo.number(j)));
            }

            if (not connected_res.at(i).isEmpty())
                lines.append(QObject::tr("  * Residue %1:%2 bonded to %3.")
                                 .arg(minfo.name(ResIdx(i)))
                                 .arg(minfo.number(ResIdx(i)))
                                 .arg(resnums.join(" ")));
        }

        if (connected_res.count() > 5)
            lines.append("   ...");
    }

    if (not connected_atoms.isEmpty())
    {
        lines.append(QObject::tr("Connected atoms:"));

        for (int i = 0; i < std::min(10, connected_atoms.count()); ++i)
        {
            QStringList atoms;

            foreach (AtomIdx j, connected_atoms.at(i))
            {
                atoms.append(::atomString(info(), j));
            }

            if (not connected_atoms.at(i).isEmpty())
                lines.append(QObject::tr("  * Atom %1 bonded to %2.")
                                 .arg(::atomString(info(), AtomIdx(i)))
                                 .arg(atoms.join(" ")));
        }

        if (connected_atoms.count() > 10)
            lines.append("   ...");
    }

    return lines.join("\n");
}

/** Return a PDB format CONECT record for this connectivity object. */
QString ConnectivityBase::toCONECT(int offset) const
{
    QStringList lines;

    if (nConnections() == 0)
        return lines.at(0);

    if (not connected_atoms.isEmpty())
    {
        // Create a set to store the bonds that have been seen. A pair is
        // used to identify the atom indices in the bond.
        QSet<QPair<int, int>> seen_bonds;

        for (int i = 0; i < connected_atoms.count(); ++i)
        {
            QStringList atoms;

            // Convert set to a list and sort the atom indices.
            auto list = connected_atoms[i].values();
            std::sort(list.begin(), list.end());

            // Loop over each atom and add to the record.
            foreach (AtomIdx j, list)
            {
                // Store outer loop index.
                int x = i;

                // Convert to int.
                int y = j.value();

                // Make sure the smaller index is the lead.
                if (x > y)
                {
                    int tmp = x;
                    x = y;
                    y = tmp;
                }

                // Create the bond index pair.
                QPair<int, int> bond_pair(x, y);

                // Don't duplicate records.
                if (not seen_bonds.contains(bond_pair))
                {
                    atoms.append(QString("%1").arg(j.value() + 1 + offset).rightJustified(4));
                    seen_bonds.insert(bond_pair);
                }
            }

            // Append the record.
            if (not connected_atoms.at(i).isEmpty() and not atoms.isEmpty())
                lines.append(QObject::tr("CONECT %1 %2")
                                 .arg(QString("%1").arg(AtomIdx(i + 1 + offset).value()).rightJustified(4))
                                 .arg(atoms.join(" ")));
        }
    }

    // Ensure that lines are 80 characters wide.
    for (auto &line : lines)
        line.resize(80, ' ');

    return lines.join("\n");
}

/** Return the total number of connections between atoms
    in this connectivity object */
int ConnectivityBase::nConnections() const
{
    int n = 0;

    const QSet<AtomIdx> *connected_atoms_array = connected_atoms.constData();
    int nats = connected_atoms.count();

    for (int i = 0; i < nats; ++i)
    {
        foreach (AtomIdx j, connected_atoms_array[i])
        {
            if (i < j)
            {
                // only count connections when the i atom is less
                // than j - this avoids double counting the bond
                //  i-j and j-i
                ++n;
            }
        }
    }

    return n;
}

/** Return the indicies of atoms connected to the atom at index 'atomidx'.
    This returns an empty set if there are no atoms connected to
    this atom

    \throw SireError::invalid_index
*/
const QSet<AtomIdx> &ConnectivityBase::connectionsTo(AtomIdx atomidx) const
{
    return connected_atoms.constData()[atomidx.map(connected_atoms.count())];
}

/** Return the indicies of atoms connected to the atom identified
    by 'resid' - this returns an empty set if there are no connections
    to this atom

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
const QSet<AtomIdx> &ConnectivityBase::connectionsTo(const AtomID &atomid) const
{
    return connected_atoms.constData()[info().atomIdx(atomid)];
}

/** Return the indicies of the residues connected to the residue at
    index 'residx'. This returns an empty set if there are no residues
    connected to this residue

    \throw SireError::invalid_index
*/
const QSet<ResIdx> &ConnectivityBase::connectionsTo(ResIdx residx) const
{
    return connected_res.constData()[residx.map(connected_res.count())];
}

/** Return the indicies of the residues connectd to the residue
    identified by 'resid'. This returns an empty set if there are
    no residues connected to this residue.

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
const QSet<ResIdx> &ConnectivityBase::connectionsTo(const ResID &resid) const
{
    return connected_res.constData()[info().resIdx(resid)];
}

/** Return the number of connections to the atom at index 'atomidx'

    \throw SireError::index_error
*/
int ConnectivityBase::nConnections(AtomIdx atomidx) const
{
    return this->connectionsTo(atomidx).count();
}

/** Return the number of connections to the atom with ID 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(const AtomID &atomid) const
{
    return this->connectionsTo(atomid).count();
}

/** Return the number of connections to the residue at index 'residx'

    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(ResIdx residx) const
{
    return this->connectionsTo(residx).count();
}

/** Return the number of connections to the residue identified
    by 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(const ResID &resid) const
{
    return this->connectionsTo(resid).count();
}

/** Return the number of atom connections between the residues at
    indicies 'res0' and 'res1'

    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(ResIdx res0, ResIdx res1) const
{
    if (not this->areConnected(res0, res1))
        return 0;

    int n = 0;

    QList<AtomIdx> atoms0 = info().getAtomsIn(res0);
    QList<AtomIdx> atoms1 = info().getAtomsIn(res1);

    foreach (AtomIdx atom0, atoms0)
    {
        const QSet<AtomIdx> &connected = this->connectionsTo(atom0);

        foreach (AtomIdx atom1, atoms1)
        {
            if (connected.contains(atom1))
                ++n;
        }
    }

    return n;
}

/** Return the number of atom connections between the residues
    identified by 'res0' and 'res1'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
int ConnectivityBase::nConnections(const ResID &res0, const ResID &res1) const
{
    return this->nConnections(info().resIdx(res0), info().resIdx(res1));
}

/** Return whether or not the atoms at indicies 'atom0' and 'atom1'
    are connected

    \throw SireError::invalid_index
*/
bool ConnectivityBase::areConnected(AtomIdx atom0, AtomIdx atom1) const
{
    return this->connectionsTo(atom0).contains(AtomIdx(atom1.map(connected_atoms.count())));
}

/** Return whether or not the atoms identified by 'atom0' and 'atom1'
    are connected

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ConnectivityBase::areConnected(const AtomID &atom0, const AtomID &atom1) const
{
    return this->connectionsTo(atom0).contains(info().atomIdx(atom1));
}

/** Return whether or not the two atoms are bonded together */
bool ConnectivityBase::areBonded(AtomIdx atom0, AtomIdx atom1) const
{
    return areConnected(atom0, atom1);
}

/** Return whether or not the two atoms are angled together */
bool ConnectivityBase::areAngled(AtomIdx atom0, AtomIdx atom2) const
{
    atom0 = minfo.atomIdx(atom0);
    atom2 = minfo.atomIdx(atom2);

    if (atom0 == atom2)
        return false;

    foreach (const AtomIdx &atom1, connected_atoms[atom0.value()])
    {
        if (connected_atoms[atom2.value()].contains(atom1))
            return true;
    }

    return false;
}

/** Return whether or not the two atoms are dihedraled together */
bool ConnectivityBase::areDihedraled(AtomIdx atom0, AtomIdx atom3) const
{
    atom0 = minfo.atomIdx(atom0);
    atom3 = minfo.atomIdx(atom3);

    if (atom0 == atom3)
        return false;

    foreach (const AtomIdx &atom1, connected_atoms[atom0.value()])
    {
        if (atom1.value() != atom3.value())
        {
            foreach (const AtomIdx &atom2, connected_atoms[atom3.value()])
            {
                if (atom2.value() != atom0.value())
                {
                    if (connected_atoms[atom1.value()].contains(atom2))
                        return true;
                }
            }
        }
    }

    return false;
}

/** Return whether or not the two atoms are bonded together */
bool ConnectivityBase::areBonded(const AtomID &atom0, const AtomID &atom1) const
{
    return areBonded(minfo.atomIdx(atom0), minfo.atomIdx(atom1));
}

/** Return whether or not the two atoms are angled together */
bool ConnectivityBase::areAngled(const AtomID &atom0, const AtomID &atom2) const
{
    return areAngled(minfo.atomIdx(atom0), minfo.atomIdx(atom2));
}

/** Return whether or not the two atoms are bonded together */
bool ConnectivityBase::areDihedraled(const AtomID &atom0, const AtomID &atom3) const
{
    return areDihedraled(minfo.atomIdx(atom0), minfo.atomIdx(atom3));
}

/** Return the connection type of the passed two atoms. This returns;

    0 = the atoms are "non-bonded", meaning no bond, angle or dihedral relationship
    1 = it is the same atom
    2 = the two atoms are directly bonded
    3 = the two atoms are part of an angle
    4 = the two atoms are part of a dihedral

    The closest connection is returned, e.g. atom pairs in rings that are simultanouesly
    bonded, angled and dihedraled are reported as bonded.
*/
int ConnectivityBase::connectionType(AtomIdx atom0, AtomIdx atom1) const
{
    atom0 = atom0.map(minfo.nAtoms());
    atom1 = atom1.map(minfo.nAtoms());

    if (atom0 == atom1)
        return 1;

    else if (areBonded(atom0, atom1))
        return 2;

    else if (areAngled(atom0, atom1))
        return 3;

    else if (areDihedraled(atom0, atom1))
        return 4;

    else
        return 0;
}

/** Return the connection type of the passed two atoms. This returns;

    0 = the atoms are "non-bonded", meaning no bond, angle or dihedral relationship
    1 = it is the same atom
    2 = the two atoms are directly bonded
    3 = the two atoms are part of an angle
    4 = the two atoms are part of a dihedral

    The closest connection is returned, e.g. atom pairs in rings that are simultanouesly
    bonded, angled and dihedraled are reported as bonded.
*/
int ConnectivityBase::connectionType(const AtomID &atom0, const AtomID &atom1) const
{
    return this->connectionType(minfo.atomIdx(atom0), minfo.atomIdx(atom1));
}

/** Return whether or not the residues at indicies 'res0' and 'res1'
    are connected

    \throw SireError::invalid_index
*/
bool ConnectivityBase::areConnected(ResIdx res0, ResIdx res1) const
{
    return this->connectionsTo(res0).contains(ResIdx(res1.map(connected_res.count())));
}

/** Return whether the residues identified by 'res0' and 'res1' are connected */
bool ConnectivityBase::areConnected(const ResID &res0, const ResID &res1) const
{
    return this->connectionsTo(res0).contains(info().resIdx(res1));
}

/** Return whether or not the CutGroups at indicies 'cg0' and 'cg1' are
    connected */
bool ConnectivityBase::areConnected(CGIdx cg0, CGIdx cg1) const
{
    cg0 = cg0.map(minfo.nCutGroups());
    cg1 = cg1.map(minfo.nCutGroups());

    if (cg0 == cg1)
    {
        return true;
    }
    else if (minfo.isResidueCutting())
    {
        // can work from the residues
        return this->areConnected(ResIdx(cg0), ResIdx(cg1));
    }
    else
    {
        // need to look at all pairs of atoms in both CutGroups
        for (const auto &atom0 : minfo.getAtomsIn(cg0))
        {
            for (const auto &atom1 : minfo.getAtomsIn(cg1))
            {
                if (this->areConnected(atom0, atom1))
                    return true;
            }
        }

        return false;
    }
}

/** Return whether or not the CutGroups at indicies 'cg0' and 'cg1' are
    connected */
bool ConnectivityBase::areConnected(const CGID &cg0, const CGID &cg1) const
{
    return this->areConnected(minfo.cgIdx(cg0), minfo.cgIdx(cg1));
}

/** Non-checking version of Connectivity::connectedTo(AtomIdx) */
const QSet<AtomIdx> &ConnectivityBase::_pvt_connectedTo(AtomIdx atom) const
{
    return connected_atoms.constData()[atom];
}

/** Internal recursive function used to find all paths between two atoms */
QList<QList<AtomIdx>> ConnectivityBase::_pvt_findPaths(AtomIdx cursor, const AtomIdx end_atom, QSet<AtomIdx> &done,
                                                       int max_length) const
{
    // create the list containing all paths from the cursor atom to the end atom
    QList<QList<AtomIdx>> all_paths;

    if (not done.contains(cursor))
    {
        // we have not traced through this atom before...
        done.insert(cursor);

        // loop through all atoms bonded to the cursor
        foreach (const AtomIdx &bonded_to_cursor, this->_pvt_connectedTo(cursor))
        {
            if (bonded_to_cursor == end_atom)
            {
                // we have found a path to the end atom. Return a single list containing
                // cursor and end_atom, so that the functions that call this can then add their
                // atoms to create all of the paths
                QList<AtomIdx> path;
                path.append(cursor);
                path.append(end_atom);
                all_paths.append(path);
            }
            else
            {
                if (max_length >= 0)
                {
                    if (done.count() + 1 >= max_length)
                    {
                        continue;
                    }
                }

                QSet<AtomIdx> new_done = done;

                QList<QList<AtomIdx>> paths = this->_pvt_findPaths(bonded_to_cursor, end_atom, new_done, max_length);

                if (not paths.isEmpty())
                {
                    for (QList<QList<AtomIdx>>::iterator it = paths.begin(); it != paths.end(); ++it)
                    {
                        (*it).prepend(cursor);
                        all_paths.append(*it);
                    }
                }
            }
        }
    }

    return all_paths;
}

/** Return all possible bonded paths between two atoms. This returns an empty
    list if there are no bonded paths between the two atoms */
QList<QList<AtomIdx>> ConnectivityBase::findPaths(AtomIdx atom0, AtomIdx atom1) const
{
    atom0 = atom0.map(minfo.nAtoms());
    atom1 = atom1.map(minfo.nAtoms());

    if (atom0 == atom1)
        return QList<QList<AtomIdx>>();

    QSet<AtomIdx> done;
    done.reserve(minfo.nAtoms());

    return this->_pvt_findPaths(atom0, atom1, done);
}

/** Return all possible bonded paths between two atoms where the path has
    a maximum length. This returns an empty list if there are no bonded
    paths between the two atoms */
QList<QList<AtomIdx>> ConnectivityBase::findPaths(AtomIdx atom0, AtomIdx atom1, int max_length) const
{
    atom0 = atom0.map(minfo.nAtoms());
    atom1 = atom1.map(minfo.nAtoms());

    if (atom0 == atom1)
        return QList<QList<AtomIdx>>();

    QSet<AtomIdx> done;

    return this->_pvt_findPaths(atom0, atom1, done, max_length);
}

/** Find the shortest bonded path between two atoms. This returns an empty
    list if there is no bonded path between these two atoms */
QList<AtomIdx> ConnectivityBase::findPath(AtomIdx atom0, AtomIdx atom1) const
{
    QList<QList<AtomIdx>> paths = findPaths(atom0, atom1);

    // Sort the paths so that the result is reproducible.
    std::stable_sort(paths.begin(), paths.end());

    QList<AtomIdx> shortest;

    foreach (const QList<AtomIdx> &path, paths)
    {
        if (shortest.isEmpty())
            shortest = path;

        else if (shortest.count() > path.count())
            shortest = path;
    }

    return shortest;
}

/** Find the shortest bonded path between two atoms where the path has
    a maximum length. This returns an empty list if there is no bonded
    path between these two atoms */
QList<AtomIdx> ConnectivityBase::findPath(AtomIdx atom0, AtomIdx atom1, int max_length) const
{
    QList<QList<AtomIdx>> paths = findPaths(atom0, atom1, max_length);

    // Sort the paths so that the result is reproducible.
    std::stable_sort(paths.begin(), paths.end());

    QList<AtomIdx> shortest;

    foreach (const QList<AtomIdx> &path, paths)
    {
        if (shortest.isEmpty())
            shortest = path;

        else if (shortest.count() > path.count())
            shortest = path;
    }

    return shortest;
}

/** Return all possible bonded paths between two atoms. This returns an empty
    list if there are no bonded paths between the two atoms */
QList<AtomIdx> ConnectivityBase::findPath(const AtomID &atom0, const AtomID &atom1) const
{
    return this->findPath(minfo.atomIdx(atom0), minfo.atomIdx(atom1));
}

/** Return all possible bonded paths between two atoms. This returns an empty
    list if there are no bonded paths between the two atoms where the path has
    a maximum length.*/
QList<AtomIdx> ConnectivityBase::findPath(const AtomID &atom0, const AtomID &atom1, int max_length) const
{
    return this->findPath(minfo.atomIdx(atom0), minfo.atomIdx(atom1), max_length);
}

/** Find the shortest bonded path between two atoms. This returns an empty
    list if there is no bonded path between these two atoms */
QList<QList<AtomIdx>> ConnectivityBase::findPaths(const AtomID &atom0, const AtomID &atom1) const
{
    return this->findPaths(minfo.atomIdx(atom0), minfo.atomIdx(atom1));
}

/** Find the shortest bonded path between two atoms where the path has
    a maximum length. This returns an empty list if there is no bonded
    path between these two atoms */
QList<QList<AtomIdx>> ConnectivityBase::findPaths(const AtomID &atom0, const AtomID &atom1, int max_length) const
{
    return this->findPaths(minfo.atomIdx(atom0), minfo.atomIdx(atom1), max_length);
}

// We only consider rings consisting of less than 12 atoms.
const int MAX_RING_SIZE = 12;

/** This function returns whether or not the atom is in a ring. */
bool ConnectivityBase::inRing(AtomIdx atom) const
{
    // Loop over all atoms connected to this atom.
    for (const auto &atm : this->connectionsTo(atom))
    {
        // Find all of the paths between the two atoms.
        QList<QList<AtomIdx>> paths = findPaths(atom, atm, MAX_RING_SIZE);

        // The atom is part of a ring.
        if (paths.count() > 1)
            return true;
    }

    // If we get this far, then the atom isn't part of a ring.
    return false;
}

/** This function returns whether or not the two passed atoms are connected */
bool ConnectivityBase::inRing(AtomIdx atom0, AtomIdx atom1) const
{
    QList<QList<AtomIdx>> paths = findPaths(atom0, atom1, MAX_RING_SIZE);

    // if there is more than one path between the atoms then they must
    // be part of a ring
    return (paths.count() > 1);
}

/** This function returns whether or not the three passed atoms are connected
    via a ring */
bool ConnectivityBase::inRing(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2) const
{
    QList<QList<AtomIdx>> paths = findPaths(atom0, atom2, MAX_RING_SIZE);

    atom1 = atom1.map(minfo.nAtoms());

    if (paths.count() > 1)
    {
        // the three are part of the same ring if any of the paths contains
        // atom1
        foreach (const QList<AtomIdx> &path, paths)
        {
            if (path.contains(atom1))
                return true;
        }
    }

    return false;
}

/** This function returns whether or not the four passed atoms are connected
    via a same ring */
bool ConnectivityBase::inRing(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3) const
{
    QList<QList<AtomIdx>> paths = findPaths(atom0, atom3, MAX_RING_SIZE);

    atom1 = atom1.map(minfo.nAtoms());
    atom2 = atom2.map(minfo.nAtoms());

    if (paths.count() > 1)
    {
        // the four are part of the same ring if atom1 and atom2 are contained
        // in any of the paths
        bool have_atom1 = false;
        bool have_atom2 = false;

        foreach (const QList<AtomIdx> &path, paths)
        {
            if (path.contains(atom1))
            {
                have_atom1 = true;
                if (have_atom2)
                    return this->inRing(atom1, atom2);
            }

            if (path.contains(atom2))
            {
                have_atom2 = true;
                if (have_atom1)
                    return this->inRing(atom1, atom2);
            }
        }
    }

    return false;
}

/** This function returns whether or not the atom is in a ring */
bool ConnectivityBase::inRing(const AtomID &atom) const
{
    return this->inRing(info().atomIdx(atom));
}

/** This function returns whether or not the two passed atoms are connected
    via a ring */
bool ConnectivityBase::inRing(const AtomID &atom0, const AtomID &atom1) const
{
    return this->inRing(info().atomIdx(atom0), info().atomIdx(atom1));
}

/** This function returns whether or not the three passed atoms are connected
    via a ring */
bool ConnectivityBase::inRing(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2) const
{
    return this->inRing(info().atomIdx(atom0), info().atomIdx(atom1), info().atomIdx(atom2));
}

/** This function returns whether or not the two passed atoms are connected
    via a ring */
bool ConnectivityBase::inRing(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2, const AtomID &atom3) const
{
    return this->inRing(info().atomIdx(atom0), info().atomIdx(atom1), info().atomIdx(atom2), info().atomIdx(atom3));
}

/** This function returns whether or not the two atoms in the passed bond
    are both part of the same ring */
bool ConnectivityBase::inRing(const BondID &bond) const
{
    return this->inRing(bond.atom0(), bond.atom1());
}

/** This function returns whether or not the three atoms in the passed angle
    are all part of the same ring */
bool ConnectivityBase::inRing(const AngleID &angle) const
{
    return this->inRing(angle.atom0(), angle.atom1(), angle.atom2());
}

/** This function returns whether or not the four atoms in the passed dihedral
    are all part of the same ring */
bool ConnectivityBase::inRing(const DihedralID &dihedral) const
{
    return this->inRing(dihedral.atom0(), dihedral.atom1(), dihedral.atom2(), dihedral.atom3());
}

/** This is a recursive function that traces all atoms that can trace their bonding
    to 'strt' but that don't go through any of the atoms in 'root',
    and to add those atoms to 'group'.

    \throw SireMol::ring_error
*/
void ConnectivityBase::traceRoute(AtomIdx start, QSet<AtomIdx> &root, QSet<AtomIdx> &group) const
{
    // add this atom to the group
    group.insert(start);
    root.insert(start);

    // now see if any of its bonded atoms need to be added
    const QSet<AtomIdx> &bonded_atoms = this->_pvt_connectedTo(start);

    // loop over every bond that involves the 'start' atom
    for (QSet<AtomIdx>::const_iterator it = bonded_atoms.constBegin(); it != bonded_atoms.constEnd(); ++it)
    {
        // if this is a root atom then ignore it, as we don't
        // want to build a path through this atom
        if (root.contains(*it))
        {
            continue;
        }
        // has this atom already been selected?
        else if (group.contains(*it))
        {
            // yes, this atom is already included!
            continue;
        }
        else
        {
            // now we can trace the atoms from the 'other' atom...
            this->traceRoute(*it, root, group);
        }
    }

    // ok - we have added all of the atoms that are connected to this atom. We
    // have finished with this atom, so we can return.
    return;
}

/** This is a recursive function that traces all atoms that can trace their bonding
    to 'strt' but that don't go through any of the atoms in 'root',
    and to add those atoms to 'group'. This traces only atoms that are contained
    in 'selected_atoms'

    \throw SireMol::ring_error
*/
void ConnectivityBase::traceRoute(const AtomSelection &selected_atoms, AtomIdx start, QSet<AtomIdx> &root,
                                  QSet<AtomIdx> &group) const
{
    if (not selected_atoms.selected(start))
        // the atom is not one of the selected atoms
        return;

    // add this atom to the group
    group.insert(start);
    root.insert(start);

    // now see if any of its bonded atoms need to be added
    const QSet<AtomIdx> &bonded_atoms = this->_pvt_connectedTo(start);

    // loop over every bond that involves the 'start' atom
    for (QSet<AtomIdx>::const_iterator it = bonded_atoms.constBegin(); it != bonded_atoms.constEnd(); ++it)
    {
        // if this is a root atom then ignore it, as we don't
        // want to move backwards!
        if (root.contains(*it))
            continue;

        // has this atom or residue already been selected?
        else if (group.contains(*it))
            // yes, this atom or residue is already included!
            continue;

        else
            // now we can trace the atoms from the 'other' atom...
            this->traceRoute(selected_atoms, *it, root, group);
    }

    // ok - we have added all of the atoms that are connected to this atom. We
    // have finished with this atom, so we can return.
    return;
}

/** Return the two AtomSelection objects corresponding to the atoms
    selected in 'group0' and 'group1' */
tuple<AtomSelection, AtomSelection> ConnectivityBase::selectGroups(const QSet<AtomIdx> &group0,
                                                                   const QSet<AtomIdx> &group1) const
{
    AtomSelection grp0(minfo);
    AtomSelection grp1(minfo);

    tuple<AtomSelection, AtomSelection> groups(grp0.selectOnly(group0), grp1.selectOnly(group1));

    return groups;
}

/** Split this molecule into two parts about the atoms
    atom0 and atom1. For example;

    For example;

          C1--C2--C3--C4--C5--C6--C7--C8

    Splitting C3 and C4 would result in two groups, {C1,C2,C3} and {C4,C5,C6,C7,C8}

          C1\
          C2-C4--C5
          C3/

    Splitting C4 and C5 would result in two groups, {C1,C2,C3,C4} and {C5}

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1) const
{
    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    // map the atoms
    int nats = minfo.nAtoms();
    atom0 = atom0.map(nats);
    atom1 = atom1.map(nats);

    if (atom0 == atom1)
        throw SireMol::ring_error(QObject::tr("You cannot split a molecule into two parts using the same atom! (%1).")
                                      .arg(::atomString(info(), atom0)),
                                  CODELOC);

    // make sure that there is sufficient space for the
    // selections - this prevents mallocs while tracing
    // the bonds
    group0.reserve(nats);
    group1.reserve(nats);
    root0.reserve(nats);
    root1.reserve(nats);

    // add the two atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom1);
    root0.insert(atom0);
    root0.insert(atom1);
    root1.insert(atom1);
    root1.insert(atom0);

    // add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if (bonded_atom != atom1)
        {
            this->traceRoute(bonded_atom, root0, group0);
        }
    }

    // remove atom1 from group0, in case it was found as part of a ring
    group0.remove(atom1);

    // now add the atoms bonded to atom1 to group1
    bool has_rings = false;

    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
    {
        if (bonded_atom != atom0)
        {
            if (group0.contains(bonded_atom))
            {
                has_rings = true;
            }

            this->traceRoute(bonded_atom, root1, group1);
        }
    }

    group1.remove(atom0);

    // if there is any overlap in the two sets then that means that
    // the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (bonded_atom != atom1 and group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                // qDebug() << "DISCONNECTING(0)" << minfo.name(bonded_atom) << minfo.name(atom0);
            }
        }

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
        {
            if (bonded_atom != atom0 and group0.contains(bonded_atom))
            {
                editor.disconnect(atom1, bonded_atom);
                // qDebug() << "DISCONNECTING(1)" << minfo.name(bonded_atom) << minfo.name(atom1);
            }
        }

        // release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();

        // split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;

        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( minfo.name(atom) );
        }

        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( minfo.name(atom) );
        }

        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/

        return this->selectGroups(group0, group1);
    }
}

/** Split the molecule into two parts about the bond between atom0 and atom1.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1) const
{
    return this->split(minfo.atomIdx(atom0), minfo.atomIdx(atom1));
}

/** Split the molecule into two parts about the bond 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const BondID &bond) const
{
    return this->split(bond.atom0(), bond.atom1());
}

/** Split the selected atoms of this molecule into two parts about the atoms
    atom0 and atom1. For example;

    For example;

          C1--C2--C3--C4--C5--C6--C7--C8

    Splitting C3 and C4 would result in two groups, {C1,C2,C3} and {C4,C5,C6,C7,C8}

          C1\
          C2-C4--C5
          C3/

    Splitting C4 and C5 would result in two groups, {C1,C2,C3,C4} and {C5}

    However splitting C1 and C5 would add a bond between C1 and C5. This would mean
    than C1-C4-C5 would form a ring, so an exception would be thrown.

    Note that both atom0 and atom1 *must* be selected as part of 'selected_atoms'
    or else a missing_atom exception will be thrown.

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1,
                                                            const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(minfo);

    if (selected_atoms.selectedAll())
        return this->split(atom0, atom1);

    selected_atoms.assertSelected(atom0);
    selected_atoms.assertSelected(atom1);

    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    // make sure that there is sufficient space for the
    // selections - this prevents mallocs while tracing
    // the bonds
    group0.reserve(selected_atoms.nSelected());
    group1.reserve(selected_atoms.nSelected());
    root0.reserve(selected_atoms.nSelected());
    root1.reserve(selected_atoms.nSelected());

    // map the atoms
    atom0 = atom0.map(minfo.nAtoms());
    atom1 = atom1.map(minfo.nAtoms());

    if (atom0 == atom1)
        throw SireMol::ring_error(QObject::tr("You cannot split a molecule into two parts using the same atom! (%1).")
                                      .arg(::atomString(info(), atom0)),
                                  CODELOC);

    // add the two atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom1);
    root0.insert(atom0);
    root0.insert(atom1);
    root1.insert(atom1);
    root1.insert(atom0);

    // add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if ((bonded_atom != atom1) and selected_atoms.selected(bonded_atom))
        {
            this->traceRoute(selected_atoms, bonded_atom, root0, group0);
        }
    }

    // now add the atoms bonded to atom1 to group1
    bool has_rings = false;

    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
    {
        if ((bonded_atom != atom0) and selected_atoms.selected(bonded_atom))
        {
            if (group0.contains(bonded_atom))
                has_rings = true;

            this->traceRoute(selected_atoms, bonded_atom, root1, group1);
        }
    }

    // if there is any overlap in the two sets then that means that
    // the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (bonded_atom != atom1 and group1.contains(bonded_atom))
                editor.disconnect(atom0, bonded_atom);
        }

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom1))
        {
            if (bonded_atom != atom0 and group0.contains(bonded_atom))
                editor.disconnect(atom1, bonded_atom);
        }

        // release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();

        // split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1);
    }
    else
        return this->selectGroups(group0, group1);
}

/** Split the selected atoms of this molecule about the atoms
    'atom0' and 'atom1'

    \throw SireMol::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/

tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                                                            const AtomSelection &selected_atoms) const
{
    return this->split(minfo.atomIdx(atom0), minfo.atomIdx(atom1), selected_atoms);
}

/** Split the selected atoms of this molecule into two parts
    about the bond 'bond'

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const BondID &bond,
                                                            const AtomSelection &selected_atoms) const
{
    return this->split(bond.atom0(), bond.atom1(), selected_atoms);
}

/** Split this molecule into three parts about the atoms
    'atom0', 'atom1' and 'atom2'.

    For example;

      C1   C3--C8
        \ /
        C2
       /  \
     C4   C5--C6-C7

    Splitting C5,C2,C3 would return {C5,C6,C7} in one group, and {C3,C8}
    in the other group.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2) const
{
    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    // map the atoms
    int nats = minfo.nAtoms();
    atom0 = atom0.map(nats);
    atom1 = atom1.map(nats);
    atom2 = atom2.map(nats);

    if (atom0 == atom1 or atom0 == atom2 or atom1 == atom2)
        throw SireMol::ring_error(
            QObject::tr("You cannot split a molecule into two parts using the same atoms! "
                        "(%1, %2, %3).")
                .arg(::atomString(info(), atom0), ::atomString(info(), atom1), ::atomString(info(), atom2)),
            CODELOC);

    // make sure that there is sufficient space for the
    // selections - this prevents mallocs while tracing
    // the bonds
    group0.reserve(nats);
    group1.reserve(nats);
    root0.reserve(nats);
    root1.reserve(nats);

    // add the end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom2);
    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);

    // add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if (bonded_atom != atom1)
        {
            this->traceRoute(bonded_atom, root0, group0);
        }
    }

    // now add the atoms bonded to atom1 to group1
    bool has_rings = false;

    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
    {
        if (bonded_atom != atom1)
        {
            if (group0.contains(bonded_atom))
                has_rings = true;

            this->traceRoute(bonded_atom, root1, group1);
        }
    }

    // if there is any overlap in the two sets then that means that
    // the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                // qDebug() << "DISCONNECTING(0)" << minfo.name(bonded_atom) << minfo.name(atom0);
            }
        }

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom2, bonded_atom);
                // qDebug() << "DISCONNECTING(1)" << minfo.name(bonded_atom) << minfo.name(atom2);
            }
        }

        // release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();

        // split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;

        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( minfo.name(atom) );
        }

        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( minfo.name(atom) );
        }

        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/

        return this->selectGroups(group0, group1);
    }
}

/** Split the molecule into two parts based on the three supplied atoms

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                                                            const AtomID &atom2) const
{
    return this->split(minfo.atomIdx(atom0), minfo.atomIdx(atom1), minfo.atomIdx(atom2));
}

/** Split the molecule into two parts based on the supplied angle

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AngleID &angle) const
{
    return this->split(angle.atom0(), angle.atom1(), angle.atom2());
}

/** Split the selected atoms of this molecule into three parts about the atoms
    'atom0', 'atom1' and 'atom2'.

    Note that all three atoms must be contained in the selection or else
    a missing_atom exception will be thrown

    An exception will be thrown if it is not possible to split the molecule
    unambiguously in two, as the angle is part of a ring.

    For example;

      C1   C3--C8
        \ /
        C2
       /  \
     C4   C5--C6-C7

    Splitting C5,C2,C3 would return {C5,C6,C7} in one group, and {C3,C8} in the other group.

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2,
                                                            const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(minfo);

    if (selected_atoms.selectedAll())
        return this->split(atom0, atom1, atom2);

    selected_atoms.assertSelected(atom0);
    selected_atoms.assertSelected(atom1);
    selected_atoms.assertSelected(atom2);

    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    // make sure that there is sufficient space for the
    // selections - this prevents mallocs while tracing
    // the bonds
    group0.reserve(selected_atoms.nSelected());
    group1.reserve(selected_atoms.nSelected());
    root0.reserve(selected_atoms.nSelected());
    root1.reserve(selected_atoms.nSelected());

    // map the atoms
    atom0 = atom0.map(minfo.nAtoms());
    atom1 = atom1.map(minfo.nAtoms());
    atom2 = atom2.map(minfo.nAtoms());

    if (atom0 == atom1 or atom0 == atom2 or atom1 == atom2)
        throw SireMol::ring_error(
            QObject::tr("You cannot split a molecule into two parts using the same atoms! "
                        "(%1, %2, %3).")
                .arg(::atomString(info(), atom0), ::atomString(info(), atom1), ::atomString(info(), atom2)),
            CODELOC);

    // add the two end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom2);

    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);

    // add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if ((bonded_atom != atom1) and selected_atoms.selected(bonded_atom))
        {
            this->traceRoute(selected_atoms, bonded_atom, root0, group0);
        }
    }

    // now add the atoms bonded to atom1 to group1
    bool has_rings = false;
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
    {
        if ((bonded_atom != atom1) and selected_atoms.selected(bonded_atom))
        {
            if (group0.contains(bonded_atom))
                has_rings = true;

            this->traceRoute(selected_atoms, bonded_atom, root1, group1);
        }
    }

    // if there is any overlap in the two sets then that means that
    // the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                // qDebug() << "DISCONNECTING(0)" << minfo.name(bonded_atom) << minfo.name(atom0);
            }
        }

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom2))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom2, bonded_atom);
                // qDebug() << "DISCONNECTING(1)" << minfo.name(bonded_atom) << minfo.name(atom2);
            }
        }

        // release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();

        // split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;

        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( minfo.name(atom) );
        }

        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( minfo.name(atom) );
        }

        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/

        return this->selectGroups(group0, group1);
    }
}

/** Split the selected atoms of the molecule into two groups around the
    three supplied atoms

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                                                            const AtomID &atom2,
                                                            const AtomSelection &selected_atoms) const
{
    return this->split(minfo.atomIdx(atom0), minfo.atomIdx(atom1), minfo.atomIdx(atom2), selected_atoms);
}

/** Split the selected atoms 'selected_atoms' of this molecule
    into two parts based on the angle identified in
    'angle'. This splits the molecule about atom0() and atom2()
    of the angle, ignoring atom atom1().

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AngleID &angle,
                                                            const AtomSelection &selected_atoms) const
{
    return this->split(angle.atom0(), angle.atom1(), angle.atom2(), selected_atoms);
}

/** Split this molecule into two parts based on the passed atoms.
    This splits the molecule between atom0 and atom3, ignoring
    atom1 and atom2.

    C1   C4--C5--C6
      \ /
      C2    C8--C9
     /  \  /
    C3   C7
           \
            C10--C11

    Splitting C4,C2,C7,C10 will return {C4,C5,C6} and {C10,C11}.
    If this molecule had been split by just Bond(C2,C7) using the above
    function, then the first returned group would
    be {C1,C2,C3,C4,C5,C6}, while the second group would be {C7,C8,C9,C10,C11}.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2,
                                                            AtomIdx atom3) const
{
    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    // map the atoms
    int nats = minfo.nAtoms();
    atom0 = atom0.map(nats);
    atom1 = atom1.map(nats);
    atom2 = atom2.map(nats);
    atom3 = atom3.map(nats);

    if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or atom1 == atom2 or atom1 == atom3 or atom2 == atom3)
        throw SireMol::ring_error(QObject::tr("You cannot split a molecule into two parts using the same atoms! "
                                              "(%1, %2, %3, %4).")
                                      .arg(::atomString(info(), atom0), ::atomString(info(), atom1),
                                           ::atomString(info(), atom2), ::atomString(info(), atom3)),
                                  CODELOC);

    // make sure that there is sufficient space for the
    // selections - this prevents mallocs while tracing
    // the bonds
    group0.reserve(nats);
    group1.reserve(nats);
    root0.reserve(nats);
    root1.reserve(nats);

    // add the end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom3);

    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root0.insert(atom3);

    root1.insert(atom3);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);

    // add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if (bonded_atom != atom1)
        {
            this->traceRoute(bonded_atom, root0, group0);
        }
    }

    // now add the atoms bonded to atom1 to group1
    bool has_rings = false;
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
    {
        if (bonded_atom != atom2)
        {
            if (group0.contains(bonded_atom))
                has_rings = true;

            this->traceRoute(bonded_atom, root1, group1);
        }
    }

    // if there is any overlap in the two sets then that means that
    // the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                // qDebug() << "DISCONNECTING(0)" << minfo.name(bonded_atom) << minfo.name(atom0);
            }
        }

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom3, bonded_atom);
                // qDebug() << "DISCONNECTING(1)" << minfo.name(bonded_atom) << minfo.name(atom3);
            }
        }

        // release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();

        // split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2, atom3);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;

        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( minfo.name(atom) );
        }

        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( minfo.name(atom) );
        }

        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/

        return this->selectGroups(group0, group1);
    }
}

/** Split this molecule into two parts based on the passed atoms.
    This splits the molecule between atom0 and atom3, ignoring
    atom1 and atom2.

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                                                            const AtomID &atom2, const AtomID &atom3) const
{
    return this->split(minfo.atomIdx(atom0), minfo.atomIdx(atom1), minfo.atomIdx(atom2), minfo.atomIdx(atom3));
}

/** Split this molecule into two parts based on the dihedral identified in
    'dihedral'. This splits the molecule about atom0() and atom3()
    of the dihedral, ignoring atoms atom1() and atom2().

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const DihedralID &dihedral) const
{
    return this->split(dihedral.atom0(), dihedral.atom1(), dihedral.atom2(), dihedral.atom3());
}

/** Split the selected atoms of this molecule into two parts
    based on the passed atoms.

    This splits the molecule between atom0 and atom3, ignoring
    atom1 and atom2.

    All four atoms must be selected in 'selected_atoms' or else
    a missing_atom exception will be thrown

    C1   C4--C5--C6
      \ /
      C2    C8--C9
     /  \  /
    C3   C7
           \
            C10--C11

    Splitting C4,C2,C7,C10 will return {C4,C5,C6} and {C10,C11}.
    If this molecule had been split by just Bond(C2,C7) using the above
    function, then the first returned group would
    be {C1,C2,C3,C4,C5,C6}, while the second group would be {C7,C8,C9,C10,C11}.

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3,
                                                            const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(minfo);

    if (selected_atoms.selectedAll())
        return this->split(atom0, atom1, atom2, atom3);

    selected_atoms.assertSelected(atom0);
    selected_atoms.assertSelected(atom1);
    selected_atoms.assertSelected(atom2);
    selected_atoms.assertSelected(atom3);

    QSet<AtomIdx> group0, group1;
    QSet<AtomIdx> root0, root1;

    // make sure that there is sufficient space for the
    // selections - this prevents mallocs while tracing
    // the bonds
    group0.reserve(selected_atoms.nSelected());
    group1.reserve(selected_atoms.nSelected());
    root0.reserve(selected_atoms.nSelected());
    root1.reserve(selected_atoms.nSelected());

    // map the atoms
    atom0 = atom0.map(minfo.nAtoms());
    atom1 = atom1.map(minfo.nAtoms());
    atom2 = atom2.map(minfo.nAtoms());
    atom3 = atom3.map(minfo.nAtoms());

    if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or atom1 == atom2 or atom1 == atom3 or atom2 == atom3)
        throw SireMol::ring_error(QObject::tr("You cannot split a molecule into two parts using the same atoms! "
                                              "(%1, %2, %3, %4).")
                                      .arg(::atomString(info(), atom0), ::atomString(info(), atom1),
                                           ::atomString(info(), atom2), ::atomString(info(), atom3)),
                                  CODELOC);

    // add the two end atoms to their respective groups
    group0.insert(atom0);
    group1.insert(atom3);

    root0.insert(atom0);
    root0.insert(atom1);
    root0.insert(atom2);
    root0.insert(atom3);

    root1.insert(atom3);
    root1.insert(atom2);
    root1.insert(atom1);
    root1.insert(atom0);

    // add the atoms bonded to atom0 to group0
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
    {
        if ((bonded_atom != atom1) and selected_atoms.selected(bonded_atom))
        {
            this->traceRoute(selected_atoms, bonded_atom, root0, group0);
        }
    }

    // now add the atoms bonded to atom1 to group1
    bool has_rings = true;
    foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
    {
        if ((bonded_atom != atom2) and selected_atoms.selected(bonded_atom))
        {
            if (group0.contains(bonded_atom))
                has_rings = true;

            this->traceRoute(selected_atoms, bonded_atom, root1, group1);
        }
    }

    // if there is any overlap in the two sets then that means that
    // the two atoms are part of a ring
    if (has_rings)
    {
        ConnectivityEditor editor = Connectivity(*this).edit();

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom0))
        {
            if (group1.contains(bonded_atom))
            {
                editor.disconnect(atom0, bonded_atom);
                // qDebug() << "DISCONNECTING(0)" << minfo.name(bonded_atom) << minfo.name(atom0);
            }
        }

        foreach (const AtomIdx &bonded_atom, this->_pvt_connectedTo(atom3))
        {
            if (group0.contains(bonded_atom))
            {
                editor.disconnect(atom3, bonded_atom);
                // qDebug() << "DISCONNECTING(1)" << minfo.name(bonded_atom) << minfo.name(atom3);
            }
        }

        // release memory to make sure that we don't recursively fill up the stack
        group0 = group1 = root0 = root1 = QSet<AtomIdx>();

        // split the molecule again, with the ring bonds now broken
        return editor.commit().split(atom0, atom1, atom2, atom3);
    }
    else
    {
        /*QSet<AtomName> names0;
        QSet<AtomName> names1;

        foreach (const AtomIdx &atom, group0)
        {
            names0.insert( minfo.name(atom) );
        }

        foreach (const AtomIdx &atom, group1)
        {
            names1.insert( minfo.name(atom) );
        }

        qDebug() << "group0" << Sire::toString(names0);
        qDebug() << "group1" << Sire::toString(names1);*/

        return this->selectGroups(group0, group1);
    }
}

/** Split the selected atoms 'selected_atoms' of this molecule
    into two parts based on the passed atoms. This splits
    the molecule between atom0 and atom3, ignoring atom1 and
    atom2.

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const AtomID &atom0, const AtomID &atom1,
                                                            const AtomID &atom2, const AtomID &atom3,
                                                            const AtomSelection &selected_atoms) const
{
    return this->split(minfo.atomIdx(atom0), minfo.atomIdx(atom1), minfo.atomIdx(atom2), minfo.atomIdx(atom3),
                       selected_atoms);
}

/** Split the selected atoms 'selected_atoms' of this molecule
    into two parts based on the dihedral identified in
    'dihedral'. This splits the molecule about atom0() and atom3()
    of the dihedral, ignoring atoms atom1() and atom2().

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const DihedralID &dihedral,
                                                            const AtomSelection &selected_atoms) const
{
    return this->split(dihedral.atom0(), dihedral.atom1(), dihedral.atom2(), dihedral.atom3(), selected_atoms);
}

/** Split this molecule into two parts based on the improper angle
    identified by 'improper'. This splits the molecule about
    bond between atom0() and atom1() of the improper

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const ImproperID &improper) const
{
    return this->split(improper.atom0(), improper.atom1());
}

/** Split the selected atoms in 'selected_atoms' in this molecule
    into two parts based on the improper angle
    identified by 'improper'. This splits the molecule about
    bond between atom0() and atom1() of the improper

    \throw SireError::incompatible_error
    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMol::ring_error
*/
tuple<AtomSelection, AtomSelection> ConnectivityBase::split(const ImproperID &improper,
                                                            const AtomSelection &selected_atoms) const
{
    return this->split(improper.atom0(), improper.atom1(), selected_atoms);
}

QList<SireMol::detail::IDPair> _getBonds(const QVector<QSet<AtomIdx>> &connections)
{
    const QSet<AtomIdx> *connections_array = connections.constData();
    const int nats = connections.count();

    QList<SireMol::detail::IDPair> bonds;

    bonds.reserve(nats * 4);

    if (nats > 100)
    {
        tbb::spin_mutex mutex;

        tbb::parallel_for(tbb::blocked_range<int>(0, nats), [&](const tbb::blocked_range<int> &r)
                          {
            QList<SireMol::detail::IDPair> my_bonds;
            my_bonds.reserve(4 * (r.end() - r.begin()));

            for (quint32 i = r.begin(); i < r.end(); ++i)
            {
                for (const AtomIdx &j_idx : connections_array[i])
                {
                    const quint32 j = j_idx.value();

                    if (i < j)
                    {
                        // only count connections when the i atom is less
                        // than j - this avoids double counting the bond
                        //  i-j and j-i
                        my_bonds.append(SireMol::detail::IDPair(i, j));
                    }
                }
            }

            tbb::spin_mutex::scoped_lock lock(mutex);
            bonds += my_bonds; });
    }
    else
    {
        for (quint32 i = 0; i < nats; ++i)
        {
            for (const AtomIdx &j_idx : connections_array[i])
            {
                const quint32 j = j_idx.value();

                if (i < j)
                {
                    // only count connections when the i atom is less
                    // than j - this avoids double counting the bond
                    //  i-j and j-i
                    bonds.append(SireMol::detail::IDPair(i, j));
                }
            }
        }
    }

    return bonds;
}

/** Return the list of bonds present in this connectivity */
QList<BondID> ConnectivityBase::getBonds() const
{
    auto bonds = _getBonds(this->connected_atoms);

    std::sort(bonds.begin(), bonds.end());

    QList<BondID> ret;
    ret.reserve(bonds.count());

    for (const auto &bond : bonds)
    {
        ret.append(BondID(AtomIdx(bond.atom0), AtomIdx(bond.atom1)));
    }

    return ret;
}

/** Return the list of bonds in the connectivity containing atom */
QList<BondID> ConnectivityBase::getBonds(const AtomID &atom) const
{
    const auto atoms = this->minfo.map(atom);

    QList<SireMol::detail::IDPair> bonds;

    for (const auto &bond : _getBonds(this->connected_atoms))
    {
        for (const auto &a : atoms)
        {
            quint32 i = a.value();

            if (bond.atom0 == i or bond.atom1 == i)
            {
                bonds.append(bond);
                break;
            }
        }
    }

    std::sort(bonds.begin(), bonds.end());

    QList<BondID> ret;
    ret.reserve(bonds.count());

    for (const auto &bond : bonds)
    {
        ret.append(BondID(AtomIdx(bond.atom0), AtomIdx(bond.atom1)));
    }

    return ret;
}

/** Return all of the connections that involve the passed atoms - if exclusive is true,
 *  then return only connections where both atoms are present in the list.
 */
QList<BondID> ConnectivityBase::getBonds(const QList<AtomIdx> &atoms, bool exclusive) const
{
    QList<BondID> ret;

    const QSet<AtomIdx> atoms_set(atoms.begin(), atoms.end());

    if (exclusive)
    {
        for (const auto &bond : this->getBonds())
        {
            const auto atom0 = this->minfo.atomIdx(bond.atom0());
            const auto atom1 = this->minfo.atomIdx(bond.atom1());

            if (atoms_set.contains(atom0) and atoms_set.contains(atom1))
            {
                ret.append(bond);
            }
        }
    }
    else
    {
        for (const auto &bond : this->getBonds())
        {
            const auto atom0 = this->minfo.atomIdx(bond.atom0());
            const auto atom1 = this->minfo.atomIdx(bond.atom1());

            if (atoms_set.contains(atom0) or atoms_set.contains(atom1))
            {
                ret.append(bond);
            }
        }
    }

    return ret;
}

namespace SireMol
{
    namespace detail
    {

        class IDTriple
        {
        public:
            IDTriple(quint32 a0 = 0, quint32 a1 = 0, quint32 a2 = 0) : atom0(a0), atom1(a1), atom2(a2)
            {
            }

            IDTriple(const IDTriple &other) : atom0(other.atom0), atom1(other.atom1), atom2(other.atom2)
            {
            }

            ~IDTriple()
            {
            }

            IDTriple &operator=(const IDTriple &other)
            {
                atom0 = other.atom0;
                atom1 = other.atom1;
                atom2 = other.atom2;
                return *this;
            }

            bool operator==(const IDTriple &other) const
            {
                return atom0 == other.atom0 and atom1 == other.atom1 and atom2 == other.atom2;
            }

            bool operator!=(const IDTriple &other) const
            {
                return not operator==(other);
            }

            bool operator<(const IDTriple &other) const
            {
                return atom0 < other.atom0 or ((atom0 == other.atom0) and
                                               ((atom1 < other.atom1) or ((atom1 == other.atom1) and (atom2 < other.atom2))));
            }

            quint32 atom0;
            quint32 atom1;
            quint32 atom2;
        };

        class IDQuad
        {
        public:
            IDQuad(quint32 a0 = 0, quint32 a1 = 0, quint32 a2 = 0, quint32 a3 = 0) : atom0(a0), atom1(a1), atom2(a2), atom3(a3)
            {
            }

            IDQuad(const IDQuad &other) : atom0(other.atom0), atom1(other.atom1), atom2(other.atom2), atom3(other.atom3)
            {
            }

            ~IDQuad()
            {
            }

            IDQuad &operator=(const IDQuad &other)
            {
                atom0 = other.atom0;
                atom1 = other.atom1;
                atom2 = other.atom2;
                atom3 = other.atom3;
                return *this;
            }

            bool operator==(const IDQuad &other) const
            {
                return atom0 == other.atom0 and atom1 == other.atom1 and atom2 == other.atom2 and atom3 == other.atom3;
            }

            bool operator!=(const IDQuad &other) const
            {
                return not operator==(other);
            }

            bool operator<(const IDQuad &other) const
            {
                return atom0 < other.atom0 or ((atom0 == other.atom0) and
                                               ((atom1 < other.atom1) or
                                                ((atom1 == other.atom1) and
                                                 ((atom2 < other.atom2) or ((atom2 == other.atom2 and atom3 < other.atom3))))));
            }

            quint32 atom0;
            quint32 atom1;
            quint32 atom2;
            quint32 atom3;
        };

    } // namespace detail
} // namespace SireMol

QList<SireMol::detail::IDTriple> _getAngles(const QVector<QSet<AtomIdx>> &connections)
{
    const int nats = connections.count();

    QList<SireMol::detail::IDTriple> angles;
    angles.reserve(3 * nats);

    const QSet<AtomIdx> *connections_array = connections.constData();

    if (nats > 100)
    {
        tbb::spin_mutex mutex;

        tbb::parallel_for(tbb::blocked_range<int>(0, nats), [&](const tbb::blocked_range<int> &r)
                          {
            QList<SireMol::detail::IDTriple> my_angs;
            my_angs.reserve(3 * (r.end() - r.begin()));

            for (quint32 i = r.begin(); i < r.end(); ++i)
            {
                for (const AtomIdx &j_idx : connections_array[i])
                {
                    const quint32 j = j_idx.value();

                    for (const auto &k_idx : connections_array[j])
                    {
                        quint32 k = k_idx.value();

                        if (k < i and k != j)
                        {
                            // only add angles in numeric order
                            my_angs.append(SireMol::detail::IDTriple(i, j, k));
                        }
                    }
                }
            }

            tbb::spin_mutex::scoped_lock lock(mutex);
            angles += my_angs; });
    }
    else
    {
        for (quint32 i = 0; i < nats; ++i)
        {
            for (const auto &j_idx : connections_array[i])
            {
                quint32 j = j_idx.value();

                for (const auto &k_idx : connections_array[j])
                {
                    quint32 k = k_idx.value();

                    if (k < i and k != j)
                    {
                        // only add angles in numeric order
                        angles.append(SireMol::detail::IDTriple(i, j, k));
                    }
                }
            }
        }
    }

    return angles;
}

/** Return a list of angles defined by the connectivity*/
QList<AngleID> ConnectivityBase::getAngles() const
{
    auto angles = _getAngles(this->connected_atoms);

    std::sort(angles.begin(), angles.end());

    QList<AngleID> ret;
    ret.reserve(angles.count());

    for (const auto &angle : angles)
    {
        ret.append(AngleID(AtomIdx(angle.atom0), AtomIdx(angle.atom1), AtomIdx(angle.atom2)));
    }

    return ret;
}

/** Return a list of angles defined by the connectivity that involve atom0 and atom1*/
QList<AngleID> ConnectivityBase::getAngles(const AtomID &atom0, const AtomID &atom1) const
{
    const auto atoms0 = this->minfo.map(atom0);
    const auto atoms1 = this->minfo.map(atom1);

    QList<SireMol::detail::IDTriple> angles;

    for (const auto &angle : _getAngles(this->connected_atoms))
    {
        bool found = false;

        for (const auto &a : atoms0)
        {
            quint32 i = a.value();

            if (angle.atom0 == i or angle.atom1 == i or angle.atom2 == i)
            {
                for (const auto &b : atoms1)
                {
                    quint32 j = b.value();

                    if ((angle.atom0 == i and angle.atom1 == j) or (angle.atom1 == i and angle.atom2 == j) or
                        (angle.atom0 == j and angle.atom1 == i) or (angle.atom1 == j and angle.atom2 == i))
                    {
                        found = true;
                        break;
                    }
                }
            }

            if (found)
                break;
        }

        if (found)
            angles.append(angle);
    }

    std::sort(angles.begin(), angles.end());

    QList<AngleID> ret;
    ret.reserve(angles.count());

    for (const auto &angle : angles)
    {
        ret.append(AngleID(AtomIdx(angle.atom0), AtomIdx(angle.atom1), AtomIdx(angle.atom2)));
    }

    return ret;
}

/** Return a list of angles defined by the connectivity that involve atom0*/
QList<AngleID> ConnectivityBase::getAngles(const AtomID &atom0) const
{
    const auto atoms = this->minfo.map(atom0);

    QList<SireMol::detail::IDTriple> angles;

    for (const auto &angle : _getAngles(this->connected_atoms))
    {
        for (const auto &a : atoms)
        {
            quint32 i = a.value();

            if (angle.atom0 == i or angle.atom1 == i or angle.atom2 == i)
            {
                angles.append(angle);
                break;
            }
        }
    }

    std::sort(angles.begin(), angles.end());

    QList<AngleID> ret;
    ret.reserve(angles.count());

    for (const auto &angle : angles)
    {
        ret.append(AngleID(AtomIdx(angle.atom0), AtomIdx(angle.atom1), AtomIdx(angle.atom2)));
    }

    return ret;
}

QList<SireMol::detail::IDQuad> _getDihedrals(const QVector<QSet<AtomIdx>> &connections)
{
    const int nats = connections.count();

    QList<SireMol::detail::IDQuad> dihedrals;
    dihedrals.reserve(4 * nats);

    const QSet<AtomIdx> *connections_array = connections.constData();

    if (nats > 100)
    {
        tbb::spin_mutex mutex;

        tbb::parallel_for(tbb::blocked_range<int>(0, nats), [&](const tbb::blocked_range<int> &r)
                          {
            QList<SireMol::detail::IDQuad> my_dihs;
            my_dihs.reserve(3 * (r.end() - r.begin()));

            for (quint32 i = r.begin(); i < r.end(); ++i)
            {
                for (const auto &j_idx : connections_array[i])
                {
                    quint32 j = j_idx.value();

                    for (const auto &k_idx : connections_array[j])
                    {
                        quint32 k = k_idx.value();

                        for (const auto &l_idx : connections_array[k])
                        {
                            quint32 l = l_idx.value();

                            if (l < i and l != k and l != j and k != j and k != i and j != i)
                            {
                                // only add dihedrals in numeric order
                                my_dihs.append(SireMol::detail::IDQuad(i, j, k, l));
                            }
                        }
                    }
                }
            }

            tbb::spin_mutex::scoped_lock lock(mutex);
            dihedrals += my_dihs; });
    }
    else
    {
        for (quint32 i = 0; i < nats; ++i)
        {
            for (const auto &j_idx : connections_array[i])
            {
                quint32 j = j_idx.value();

                for (const auto &k_idx : connections_array[j])
                {
                    quint32 k = k_idx.value();

                    for (const auto &l_idx : connections_array[k])
                    {
                        quint32 l = l_idx.value();

                        if (l < i and l != k and l != j and k != j and k != i and j != i)
                        {
                            // only add dihedrals in numeric order
                            dihedrals.append(SireMol::detail::IDQuad(i, j, k, l));
                        }
                    }
                }
            }
        }
    }

    return dihedrals;
}

/** Return a list of dihedrals defined by the connectivity*/
QList<DihedralID> ConnectivityBase::getDihedrals() const
{
    auto dihedrals = _getDihedrals(this->connected_atoms);

    QList<DihedralID> ret;
    ret.reserve(dihedrals.count());

    for (const auto &dihedral : dihedrals)
    {
        ret.append(DihedralID(AtomIdx(dihedral.atom0), AtomIdx(dihedral.atom1), AtomIdx(dihedral.atom2),
                              AtomIdx(dihedral.atom3)));
    }

    return ret;
}

/** Return a list of dihedrals defined by the connectivity that involve atom0, atom1 and atom2*/
QList<DihedralID> ConnectivityBase::getDihedrals(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2) const
{
    const auto atoms0 = this->minfo.map(atom0);
    const auto atoms1 = this->minfo.map(atom1);
    const auto atoms2 = this->minfo.map(atom2);

    QList<SireMol::detail::IDQuad> dihedrals;

    for (const auto &dih : _getDihedrals(this->connected_atoms))
    {
        bool found = false;

        for (const auto &atom0 : atoms0)
        {
            quint32 i = atom0.value();

            if (dih.atom0 == i or dih.atom1 == i or dih.atom2 == i or dih.atom3 == i)
            {
                for (const auto &atom1 : atoms1)
                {
                    quint32 j = atom1.value();

                    if ((dih.atom0 == i and dih.atom1 == j) or (dih.atom1 == i and dih.atom2 == j) or
                        (dih.atom2 == i and dih.atom3 == j) or (dih.atom0 == j and dih.atom1 == i) or
                        (dih.atom1 == j and dih.atom2 == i) or (dih.atom2 == j and dih.atom3 == i))
                    {
                        for (const auto &atom2 : atoms2)
                        {
                            quint32 k = atom2.value();

                            if ((dih.atom0 == i and dih.atom1 == j and dih.atom2 == k) or
                                (dih.atom1 == i and dih.atom2 == j and dih.atom3 == k) or
                                (dih.atom0 == k and dih.atom2 == j and dih.atom3 == i) or
                                (dih.atom1 == k and dih.atom2 == j and dih.atom3 == i))
                            {
                                found = true;
                                break;
                            }
                        }
                    }

                    if (found)
                        break;
                }
            }

            if (found)
                break;
        }

        if (found)
            dihedrals.append(dih);
    }

    QList<DihedralID> ret;
    ret.reserve(dihedrals.count());

    for (const auto &dihedral : dihedrals)
    {
        ret.append(DihedralID(AtomIdx(dihedral.atom0), AtomIdx(dihedral.atom1), AtomIdx(dihedral.atom2),
                              AtomIdx(dihedral.atom3)));
    }

    return ret;
}

/** Return a list of dihedrals defined by the connectivity that involve atom0 and atom1*/
QList<DihedralID> ConnectivityBase::getDihedrals(const AtomID &atom0, const AtomID &atom1) const
{
    const auto atoms0 = this->minfo.map(atom0);
    const auto atoms1 = this->minfo.map(atom1);

    QList<SireMol::detail::IDQuad> dihedrals;

    for (const auto &dih : _getDihedrals(this->connected_atoms))
    {
        bool found = false;

        for (const auto &atom0 : atoms0)
        {
            quint32 i = atom0.value();

            if (dih.atom0 == i or dih.atom1 == i or dih.atom2 == i or dih.atom3 == i)
            {
                for (const auto &atom1 : atoms1)
                {
                    quint32 j = atom1.value();

                    if ((dih.atom0 == i and dih.atom1 == j) or (dih.atom1 == i and dih.atom2 == j) or
                        (dih.atom2 == i and dih.atom3 == j) or (dih.atom0 == j and dih.atom1 == i) or
                        (dih.atom1 == j and dih.atom2 == i) or (dih.atom2 == j and dih.atom3 == i))
                    {
                        found = true;
                        break;
                    }
                }
            }

            if (found)
                break;
        }

        if (found)
            dihedrals.append(dih);
    }

    QList<DihedralID> ret;
    ret.reserve(dihedrals.count());

    for (const auto &dihedral : dihedrals)
    {
        ret.append(DihedralID(AtomIdx(dihedral.atom0), AtomIdx(dihedral.atom1), AtomIdx(dihedral.atom2),
                              AtomIdx(dihedral.atom3)));
    }

    return ret;
}

/** Return a list of dihedrals defined by the connectivity that involve atom0*/
QList<DihedralID> ConnectivityBase::getDihedrals(const AtomID &atom0) const
{
    const auto atoms0 = this->minfo.map(atom0);

    QList<SireMol::detail::IDQuad> dihedrals;

    for (const auto &dih : _getDihedrals(this->connected_atoms))
    {
        for (const auto &atom0 : atoms0)
        {
            quint32 i = atom0.value();

            if (dih.atom0 == i or dih.atom1 == i or dih.atom2 == i or dih.atom3 == i)
            {
                dihedrals.append(dih);
                break;
            }
        }
    }

    QList<DihedralID> ret;
    ret.reserve(dihedrals.count());

    for (const auto &dihedral : dihedrals)
    {
        ret.append(DihedralID(AtomIdx(dihedral.atom0), AtomIdx(dihedral.atom1), AtomIdx(dihedral.atom2),
                              AtomIdx(dihedral.atom3)));
    }

    return ret;
}

/** Return a matrix (organised by AtomIdx) that says which atoms are bonded between
    order 'start' and order 'end' (e.g. if order is two, it returns true for each atom pair that
    are bonded together, if order is three, then true for each atom pair that are
    bonded or angled together, if order is four, then true for each atom pair
    that are bonded, angled or dihedraled) */
QVector<QVector<bool>> ConnectivityBase::getBondMatrix(int start, int end) const
{
    if (start < 0)
        start = 0;

    if (end < 0)
        end = 0;

    if (start > end)
        qSwap(start, end);

    QVector<QVector<bool>> ret;

    const int nats = minfo.nAtoms();

    if (this->nConnections() == 0)
    {
        // there are no internal connections - we have to assume that
        // everything is bonded to everything else. This is because
        // we are likely to be an empty connectivity for, e.g.
        // a solvent molecule
        ret = QVector<QVector<bool>>(nats, QVector<bool>(nats, true));
        return ret;
    }

    if (nats == 0)
        return ret;

    else
    {
        ret = QVector<QVector<bool>>(nats);
        ret.squeeze();

        QVector<bool> row;

        if (start == 0)
            row = QVector<bool>(nats, true);
        else
            row = QVector<bool>(nats, false);

        row.squeeze();

        for (int i = 0; i < nats; ++i)
        {
            ret.data()[i] = row;
        }
    }

    if (start <= 0)
        return ret;

    for (int order = start; order <= end; ++order)
    {
        if (order == 1)
        {
            for (int i = 0; i < nats; ++i)
            {
                ret[i][i] = true;
            }
        }

        if (order == 2)
        {
            for (int i = 0; i < nats; ++i)
            {
                QVector<bool> &row = ret.data()[i];

                for (QSet<AtomIdx>::const_iterator it = connected_atoms[i].constBegin();
                     it != connected_atoms[i].constEnd(); ++it)
                {
                    row[it->value()] = true;
                }
            }
        }

        if (order == 3)
        {
            for (int atm0 = 0; atm0 < nats; ++atm0)
            {
                QVector<bool> &row = ret.data()[atm0];

                for (QSet<AtomIdx>::const_iterator it = connected_atoms[atm0].constBegin();
                     it != connected_atoms[atm0].constEnd(); ++it)
                {
                    const int atm1 = it->value();

                    for (QSet<AtomIdx>::const_iterator it2 = connected_atoms[atm1].constBegin();
                         it2 != connected_atoms[atm1].constEnd(); ++it2)
                    {
                        const int atm2 = it2->value();

                        if (atm2 != atm0)
                            row[atm2] = true;
                    }
                }
            }
        }

        if (order == 4)
        {
            for (int atm0 = 0; atm0 < nats; ++atm0)
            {
                QVector<bool> &row = ret.data()[atm0];

                for (QSet<AtomIdx>::const_iterator it = connected_atoms[atm0].constBegin();
                     it != connected_atoms[atm0].constEnd(); ++it)
                {
                    const int atm1 = it->value();

                    for (QSet<AtomIdx>::const_iterator it2 = connected_atoms[atm1].constBegin();
                         it2 != connected_atoms[atm1].constEnd(); ++it2)
                    {
                        const int atm2 = it2->value();

                        if (atm2 != atm0)
                        {
                            for (QSet<AtomIdx>::const_iterator it3 = connected_atoms[atm2].constBegin();
                                 it3 != connected_atoms[atm2].constEnd(); ++it3)
                            {
                                const int atm3 = it3->value();

                                if (atm3 != atm0 and atm3 != atm1 and atm3 != atm2)
                                {
                                    row[atm3] = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (order > 4)
        {
            qDebug() << "Cannot build a bond matrix for values greater than 4";
            break;
        }
    }

    return ret;
}

/** Return a matrix (organised by AtomIdx) that says which atoms are bonded up to
    order 'order' (e.g. if order is two, it returns true for each atom pair that
    are bonded together, if order is three, then true for each atom pair that are
    bonded or angled together, if order is four, then true for each atom pair
    that are bonded, angled or dihedraled) */
QVector<QVector<bool>> ConnectivityBase::getBondMatrix(int order) const
{
    if (order <= 0)
        return getBondMatrix(0, 0);
    else
        return getBondMatrix(1, order);
}

/** Return all of the property keys for all of the bonds */
QStringList ConnectivityBase::propertyKeys() const
{
    QSet<QString> keys;

    for (const auto &key : this->bond_props.keys())
    {
        const auto &props = this->bond_props[key];

        for (const auto &k : props.propertyKeys())
        {
            keys.insert(k);
        }
    }

    for (const auto &key : this->ang_props.keys())
    {
        const auto &props = this->ang_props[key];

        for (const auto &k : props.propertyKeys())
        {
            keys.insert(k);
        }
    }

    for (const auto &key : this->dih_props.keys())
    {
        const auto &props = this->dih_props[key];

        for (const auto &k : props.propertyKeys())
        {
            keys.insert(k);
        }
    }

    for (const auto &key : this->imp_props.keys())
    {
        const auto &props = this->imp_props[key];

        for (const auto &k : props.propertyKeys())
        {
            keys.insert(k);
        }
    }

    QStringList ret = keys.values();
    ret.sort();

    return ret;
}

/** Return the properties of the passed bond */
Properties ConnectivityBase::properties(const BondID &bond) const
{
    auto id = SireMol::detail::IDPair(this->minfo.atomIdx(bond.atom0()), this->minfo.atomIdx(bond.atom1()));

    return this->bond_props.value(id);
}

/** Return whether the specified bond has a property at key 'key' */
bool ConnectivityBase::hasProperty(const BondID &bond, const PropertyName &key) const
{
    return this->properties(bond).hasProperty(key);
}

/** Return the type of the property for the specified bond at key 'key' */
const char *ConnectivityBase::propertyType(const BondID &bond, const PropertyName &key) const
{
    return this->properties(bond).propertyType(key);
}

/** Return the property keys for the specified bond */
QStringList ConnectivityBase::propertyKeys(const BondID &bond) const
{
    return this->properties(bond).propertyKeys();
}

/** Return the specified property of the specified bond */
const Property &ConnectivityBase::property(const BondID &bond, const PropertyName &key) const
{
    return this->properties(bond).property(key);
}

/** Return the specified property of the specified bond, or
    'default_value' if such a property is not defined
 */
const Property &ConnectivityBase::property(const BondID &bond, const PropertyName &key,
                                           const Property &default_value) const
{
    return this->properties(bond).property(key, default_value);
}

/** Assert that the specified bond has the specified property */
void ConnectivityBase::assertHasProperty(const BondID &bond, const PropertyName &key) const
{
    if (not this->hasProperty(bond, key))
        throw SireBase::missing_property(QObject::tr("Bond %1 "
                                                     "does not have a valid property at key \"%2\".")
                                             .arg(bond.toString())
                                             .arg(key.toString()),
                                         CODELOC);
}

AngleID _to_canonical(const AngleID &angle, const MoleculeInfoData &info)
{
    auto a0 = info.atomIdx(angle.atom0());
    auto a1 = info.atomIdx(angle.atom1());
    auto a2 = info.atomIdx(angle.atom2());

    if (a0 < a2)
        qSwap(a0, a2);

    return AngleID(a0, a1, a2);
}

/** Return the properties of the passed angle */
Properties ConnectivityBase::properties(const AngleID &angle) const
{
    return this->ang_props.value(_to_canonical(angle, this->info()));
}

/** Return whether the specified angle has a property at key 'key' */
bool ConnectivityBase::hasProperty(const AngleID &angle, const PropertyName &key) const
{
    return this->properties(angle).hasProperty(key);
}

/** Return the type of the property for the specified angle at key 'key' */
const char *ConnectivityBase::propertyType(const AngleID &angle, const PropertyName &key) const
{
    return this->properties(angle).propertyType(key);
}

/** Return the property keys for the specified angle */
QStringList ConnectivityBase::propertyKeys(const AngleID &angle) const
{
    return this->properties(angle).propertyKeys();
}

/** Return the specified property of the specified angle */
const Property &ConnectivityBase::property(const AngleID &angle, const PropertyName &key) const
{
    return this->properties(angle).property(key);
}

/** Return the specified property of the specified angle, or
    'default_value' if such a property is not defined
 */
const Property &ConnectivityBase::property(const AngleID &angle, const PropertyName &key,
                                           const Property &default_value) const
{
    return this->properties(angle).property(key, default_value);
}

/** Assert that the specified angle has the specified property */
void ConnectivityBase::assertHasProperty(const AngleID &angle, const PropertyName &key) const
{
    if (not this->hasProperty(angle, key))
        throw SireBase::missing_property(QObject::tr("Angle %1 "
                                                     "does not have a valid property at key \"%2\".")
                                             .arg(angle.toString())
                                             .arg(key.toString()),
                                         CODELOC);
}

DihedralID _to_canonical(const DihedralID &dihedral, const MoleculeInfoData &info)
{
    auto a0 = info.atomIdx(dihedral.atom0());
    auto a1 = info.atomIdx(dihedral.atom1());
    auto a2 = info.atomIdx(dihedral.atom2());
    auto a3 = info.atomIdx(dihedral.atom3());

    if (a0 < a3)
    {
        qSwap(a0, a3);
        qSwap(a1, a2);
    }

    return DihedralID(a0, a1, a2, a3);
}

/** Return the properties of the passed dihedral */
Properties ConnectivityBase::properties(const DihedralID &dihedral) const
{
    return this->dih_props.value(_to_canonical(dihedral, this->info()));
}

/** Return whether the specified dihedral has a property at key 'key' */
bool ConnectivityBase::hasProperty(const DihedralID &dihedral, const PropertyName &key) const
{
    return this->properties(dihedral).hasProperty(key);
}

/** Return the type of the property for the specified dihedral at key 'key' */
const char *ConnectivityBase::propertyType(const DihedralID &dihedral, const PropertyName &key) const
{
    return this->properties(dihedral).propertyType(key);
}

/** Return the property keys for the specified dihedral */
QStringList ConnectivityBase::propertyKeys(const DihedralID &dihedral) const
{
    return this->properties(dihedral).propertyKeys();
}

/** Return the specified property of the specified dihedral */
const Property &ConnectivityBase::property(const DihedralID &dihedral, const PropertyName &key) const
{
    return this->properties(dihedral).property(key);
}

/** Return the specified property of the specified dihedral, or
    'default_value' if such a property is not defined
 */
const Property &ConnectivityBase::property(const DihedralID &dihedral, const PropertyName &key,
                                           const Property &default_value) const
{
    return this->properties(dihedral).property(key, default_value);
}

/** Assert that the specified angle has the specified property */
void ConnectivityBase::assertHasProperty(const DihedralID &dihedral, const PropertyName &key) const
{
    if (not this->hasProperty(dihedral, key))
        throw SireBase::missing_property(QObject::tr("Dihedral %1 "
                                                     "does not have a valid property at key \"%2\".")
                                             .arg(dihedral.toString())
                                             .arg(key.toString()),
                                         CODELOC);
}

ImproperID _to_canonical(const ImproperID &improper, const MoleculeInfoData &info)
{
    auto a0 = info.atomIdx(improper.atom0());
    auto a1 = info.atomIdx(improper.atom1());
    auto a2 = info.atomIdx(improper.atom2());
    auto a3 = info.atomIdx(improper.atom3());

    return ImproperID(a0, a1, a2, a3);
}

/** Return the properties of the passed improper */
Properties ConnectivityBase::properties(const ImproperID &improper) const
{
    return this->imp_props.value(_to_canonical(improper, this->info()));
}

/** Return whether the specified improper has a property at key 'key' */
bool ConnectivityBase::hasProperty(const ImproperID &improper, const PropertyName &key) const
{
    return this->properties(improper).hasProperty(key);
}

/** Return the type of the property for the specified improper at key 'key' */
const char *ConnectivityBase::propertyType(const ImproperID &improper, const PropertyName &key) const
{
    return this->properties(improper).propertyType(key);
}

/** Return the property keys for the specified improper */
QStringList ConnectivityBase::propertyKeys(const ImproperID &improper) const
{
    return this->properties(improper).propertyKeys();
}

/** Return the specified property of the specified improper */
const Property &ConnectivityBase::property(const ImproperID &improper, const PropertyName &key) const
{
    return this->properties(improper).property(key);
}

/** Return the specified property of the specified improper, or
    'default_value' if such a property is not defined
 */
const Property &ConnectivityBase::property(const ImproperID &improper, const PropertyName &key,
                                           const Property &default_value) const
{
    return this->properties(improper).property(key, default_value);
}

/** Assert that the specified angle has the specified property */
void ConnectivityBase::assertHasProperty(const ImproperID &improper, const PropertyName &key) const
{
    if (not this->hasProperty(improper, key))
        throw SireBase::missing_property(QObject::tr("Improper %1 "
                                                     "does not have a valid property at key \"%2\".")
                                             .arg(improper.toString())
                                             .arg(key.toString()),
                                         CODELOC);
}

/** Return whether or not this atom is in a ring */
static bool is_on_ring(const AtomIdx &atom, const Connectivity &conn)
{
    if (conn.inRing(atom))
    {
        // loop over all atoms connected to this atom
        for (const auto &neighbour : conn.connectionsTo(atom))
        {
            // if the neighbour is in the ring, then this atom is in the ring
            if (not(conn.inRing(neighbour, atom)))
                return true;
        }
    }

    return false;
}

/** Return whether or not there is a change in ring between the passed
 *  two atoms in the passed two connectivities */
static bool is_ring_size_changed(const Connectivity &conn0,
                                 const Connectivity &conn1,
                                 const AtomIdxMappingEntry &atom0,
                                 const AtomIdxMappingEntry &atom1,
                                 int max_ring_size = 12)
{
    // Have a ring changed size? If so, then the minimum path size between
    // two atoms will have changed.

    // Work out the paths connecting the atoms in the two end states.
    auto paths0 = conn0.findPaths(atom0.atomIdx0(), atom1.atomIdx0(), max_ring_size);
    auto paths1 = conn1.findPaths(atom0.atomIdx1(), atom1.atomIdx1(), max_ring_size);

    // Initialise the ring size in each end state.
    auto ring0 = -1;
    auto ring1 = -1;

    // Determine the minimum path in the lambda = 0 state.
    if (paths0.count() > 1)
    {
        QVector<int> path_lengths0;

        for (const auto &path : paths0)
        {
            path_lengths0.append(path.count());
        }

        ring0 = *(std::min_element(path_lengths0.begin(), path_lengths0.end()));
    }

    if (ring0 == -1)
        return false;

    // Determine the minimum path in the lambda = 1 state.
    if (paths1.count() > 1)
    {
        QVector<int> path_lengths1;

        for (const auto &path : paths1)
        {
            path_lengths1.append(path.count());
        }

        ring1 = *(std::min_element(path_lengths1.begin(), path_lengths1.end()));
    }

    // Return whether the ring has changed size.
    return ring0 != ring1;
}

/** Return whether any ring that both the atoms are on is broken/formed
 *  during the merge */
static bool is_ring_broken(const Connectivity &conn0,
                           const Connectivity &conn1,
                           const AtomIdxMappingEntry &atom0,
                           const AtomIdxMappingEntry &atom1)
{
    // Have we opened/closed a ring? This means that both atoms are part of a
    // ring in one end state (either in it, or on it), whereas at least one
    // isn't in the other state.

    if (not(atom0.isMappedInBoth() and atom1.isMappedInBoth()))
        // either atom isn't in both end states, so cannot be part
        // of a ring that is in both end states...
        return false;

    const auto idx0 = atom0.atomIdx0();
    const auto idy0 = atom1.atomIdx0();

    const auto idx1 = atom0.atomIdx1();
    const auto idy1 = atom1.atomIdx1();

    // Whether each atom is in a ring in both end states.
    const auto in_ring_idx0 = conn0.inRing(idx0);
    const auto in_ring_idy0 = conn0.inRing(idy0);
    const auto in_ring_idx1 = conn1.inRing(idx1);
    const auto in_ring_idy1 = conn1.inRing(idy1);

    // Whether each atom is on a ring in both end states.
    const auto on_ring_idx0 = is_on_ring(idx0, conn0);
    const auto on_ring_idy0 = is_on_ring(idy0, conn0);
    const auto on_ring_idx1 = is_on_ring(idx1, conn1);
    const auto on_ring_idy1 = is_on_ring(idy1, conn1);

    // Both atoms are in a ring in one end state and at least one isn't in the other.
    if ((in_ring_idx0 & in_ring_idy0) ^ (in_ring_idx1 & in_ring_idy1))
        return true;

    // Both atoms are on a ring in one end state and at least one isn't in the other.
    if ((on_ring_idx0 & on_ring_idy0 & (conn0.connectionType(idx0, idy0) == 4)) ^ (on_ring_idx1 & on_ring_idy1 & (conn1.connectionType(idx1, idy1) == 4)))
    {
        // Make sure that the change isn't a result of ring growth, i.e. one of
        // the atoms isn't in a ring in one end state, while its "on" ring status
        // has changed between states.
        if (not((in_ring_idx0 | in_ring_idx1) & (on_ring_idx0 ^ on_ring_idx1) or (in_ring_idy0 | in_ring_idy1) & (on_ring_idy0 ^ on_ring_idy1)))
        {
            return true;
        }
    }

    // Both atoms are in or on a ring in one state and at least one isn't in the other.
    if (((in_ring_idx0 | on_ring_idx0) & (in_ring_idy0 | on_ring_idy0) & (conn0.connectionType(idx0, idy0) == 3)) ^
        ((in_ring_idx1 | on_ring_idx1) & (in_ring_idy1 | on_ring_idy1) & (conn1.connectionType(idx1, idy1) == 3)))
    {
        auto iscn0 = conn0.connectionsTo(idx0);
        iscn0.intersect(conn0.connectionsTo(idy0));

        if (iscn0.count() != 1)
            return true;

        auto common_idx = *(iscn0.constBegin());
        iscn0.remove(common_idx);

        const auto in_ring_bond0 = conn0.inRing(idx0, common_idx) or conn0.inRing(idy0, common_idx);

        auto iscn1 = conn1.connectionsTo(idx1);
        iscn1.intersect(conn1.connectionsTo(idy1));

        if (iscn1.count() != 1)
            return true;

        common_idx = *(iscn1.constBegin());
        iscn1.remove(common_idx);

        const auto in_ring_bond1 = conn1.inRing(idx1, common_idx) or conn1.inRing(idy1, common_idx);

        if (in_ring_bond0 ^ in_ring_bond1)
            return true;
    }

    // If we get this far, then a ring wasn't broken.
    return false;
}

/** Merge this property with another property */
PropertyList ConnectivityBase::merge(const MolViewProperty &other,
                                     const AtomIdxMapping &mapping,
                                     const QString &ghost,
                                     const SireBase::PropertyMap &map) const
{
    if (not other.isA<ConnectivityBase>())
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

    const ConnectivityBase &ref = *this;
    const ConnectivityBase &pert = other.asA<ConnectivityBase>();

    auto prop0 = Connectivity(ref).edit();
    auto prop1 = Connectivity(ref).edit();

    // the prop0 properties are already correct

    // the prop1 properties are made by finding all of the atoms that
    // are involved in bonds in 'pert' and removing any bonds involving
    // only those atoms from 'prop1', and then adding back the matching
    // bonds from 'pert'. Use 'true' to only remove bonds where both
    // atoms are in the mapping
    prop1.disconnect(mapping.mappedIn1(), true);

    // get the mapping from the perturbed to reference states, including
    // atoms that don't exist in the reference state. In all cases,
    // the values are the indexes in the merged molecule
    auto map1to0 = mapping.map1to0(true);

    // now find all of the bonds in 'pert' where both atoms in the
    // bond are in map1to0.keys() - i.e. exist and are mapped from
    // the perturbed state
    const auto pert_bonds = pert.getBonds(map1to0.keys(), true);

    // connect those bonds together
    for (const auto &pert_bond : pert_bonds)
    {
        const auto atom0 = map1to0.value(info().atomIdx(pert_bond.atom0()));
        const auto atom1 = map1to0.value(info().atomIdx(pert_bond.atom1()));

        prop1.connect(atom0, atom1);

        if (mapping.isUnmappedIn0(atom0) or mapping.isUnmappedIn0(atom1))
        {
            // the prop0 properties are nearly correct - we just need to add
            // in a connection from 'pert' that involve the atoms that are not mapped
            // in the reference state - this way, those added atoms are
            // connected to the reference atoms
            prop0.connect(atom0, atom1);
        }
    }

    // now add in the connections to the perturbed state from the reference
    // state for any atoms that aren't mapped to the perturbed state.
    // This way, the removed atoms are connected to the perturbed atoms
    auto map0to1 = mapping.map0to1(true);

    const auto ref_bonds = ref.getBonds(map0to1.keys(), true);

    for (const auto &ref_bond : ref_bonds)
    {
        const auto atom0 = info().atomIdx(ref_bond.atom0());
        const auto atom1 = info().atomIdx(ref_bond.atom1());

        if (mapping.isUnmappedIn1(atom0) or mapping.isUnmappedIn1(atom1))
        {
            prop1.connect(atom0, atom1);
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

    SireBase::PropertyList ret;

    ret.append(prop0.commit());
    ret.append(prop1.commit());

    if (allow_ring_breaking and allow_ring_size_change)
    {
        // nothing more to do
        return ret;
    }

    // we need to check that the merge doesn't break or change rings - do this
    // by looping over pairs of atoms mapped in both states and checking,
    // if they are in a ring, if that ring has changed size or been broken or
    // formed by the merge
    for (auto it0 = mapping.constBegin(); it0 != mapping.constEnd(); ++it0)
    {
        const auto &atom0 = *it0;

        if (not atom0.isMappedInBoth())
            // this cannot be part of a ring, as it is not in both states
            continue;

        for (auto it1 = it0 + 1; it1 != mapping.constEnd(); ++it1)
        {
            const auto &atom1 = *it1;

            if (not atom1.isMappedInBoth())
                // this cannot be part of a ring, as it is not in both states
                continue;

            if (not allow_ring_size_change)
            {
                if (is_ring_size_changed(ref, pert, atom0, atom1))
                {
                    throw SireError::incompatible_error(QObject::tr("The merge has changed the size of a ring. To allow this "
                                                                    "perturbation, set the 'allow_ring_size_change' option "
                                                                    "to 'True'. Be aware that this perturbation may not work "
                                                                    "and a transition through an intermediate state may be "
                                                                    "preferable."),
                                                        CODELOC);
                }
            }

            if (not allow_ring_breaking)
            {
                if (is_ring_broken(ref, pert, atom0, atom1))
                {
                    throw SireError::incompatible_error(QObject::tr("The merge has changed the molecular connectivity "
                                                                    "but a ring didn't open/close or change size. "
                                                                    "If you want to proceed with this mapping pass "
                                                                    "'force=True'. You are warned that the resulting "
                                                                    "perturbation will likely be unstable."),
                                                        CODELOC);
                }
            }
        }
    }

    return ret;
}

/////////
///////// Implementation of Connectivity
/////////

static const RegisterMetaType<Connectivity> r_connectivity;

QDataStream &operator<<(QDataStream &ds, const Connectivity &conn)
{
    writeHeader(ds, r_connectivity, 1);

    ds << static_cast<const ConnectivityBase &>(conn);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Connectivity &conn)
{
    VersionID v = readHeader(ds, r_connectivity);

    if (v == 1)
    {
        ds >> static_cast<ConnectivityBase &>(conn);
    }
    else
        throw version_error(v, "1", r_connectivity, CODELOC);

    return ds;
}

/** Reduce the memory usage of this object to a minimum */
void Connectivity::squeeze()
{
    connected_atoms.squeeze();

    QSet<AtomIdx> *connected_atoms_array = connected_atoms.data();
    int nats = connected_atoms.count();

    for (int i = 0; i < nats; ++i)
    {
        connected_atoms_array[i].squeeze();
    }

    connected_res.squeeze();

    QSet<ResIdx> *connected_res_array = connected_res.data();
    int nres = connected_res.count();

    for (int i = 0; i < nres; ++i)
    {
        connected_res_array[i].squeeze();
    }
}

/** Null constructor */
Connectivity::Connectivity() : ConcreteProperty<Connectivity, ConnectivityBase>()
{
}

/** Construct the connectivity for the molecule whose data
    is in 'moldata' */
Connectivity::Connectivity(const MoleculeData &moldata) : ConcreteProperty<Connectivity, ConnectivityBase>(moldata)
{
}

/** Construct the connectivity for the passed molecule info */
Connectivity::Connectivity(const MoleculeInfo &molinfo) : ConcreteProperty<Connectivity, ConnectivityBase>(molinfo)
{
}

/** Construct the connectivity for the molecule viewed in the
    passed view. This automatically uses the bond hunting
    function to add all of the bonds for the atoms in this view */
Connectivity::Connectivity(const MoleculeView &molview, const BondHunter &bondhunter, const PropertyMap &map)
    : ConcreteProperty<Connectivity, ConnectivityBase>()
{
    this->operator=(bondhunter(molview, map));
}

/** Construct the connectivity from the passed editor */
Connectivity::Connectivity(const ConnectivityEditor &editor) : ConcreteProperty<Connectivity, ConnectivityBase>(editor)
{
    this->squeeze();
}

/** Private constructor allowing a ConnectivityBase to become a Connectivity */
Connectivity::Connectivity(const ConnectivityBase &base) : ConcreteProperty<Connectivity, ConnectivityBase>(base)
{
}

/** Copy constructor */
Connectivity::Connectivity(const Connectivity &other) : ConcreteProperty<Connectivity, ConnectivityBase>(other)
{
}

/** Destructor */
Connectivity::~Connectivity()
{
}

/** Copy assignment from another Connectivity object */
Connectivity &Connectivity::operator=(const Connectivity &other)
{
    ConnectivityBase::operator=(other);
    return *this;
}

/** Copy assignment from a ConnectivityEditor */
Connectivity &Connectivity::operator=(const ConnectivityEditor &editor)
{
    ConnectivityBase::operator=(editor);
    this->squeeze();
    return *this;
}

/** Comparison operator */
bool Connectivity::operator==(const Connectivity &other) const
{
    return ConnectivityBase::operator==(other);
}

/** Comparison operator */
bool Connectivity::operator!=(const Connectivity &other) const
{
    return ConnectivityBase::operator!=(other);
}

/** Return an editor that can edit a copy of this connectivity */
ConnectivityEditor Connectivity::edit() const
{
    return ConnectivityEditor(*this);
}

const char *Connectivity::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Connectivity>());
}

/////////
///////// Implementation of ConnectivityEditor
/////////

static const RegisterMetaType<ConnectivityEditor> r_conneditor;

QDataStream &operator<<(QDataStream &ds, const ConnectivityEditor &conn)
{
    writeHeader(ds, r_conneditor, 1);

    ds << static_cast<const ConnectivityBase &>(conn);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, ConnectivityEditor &conn)
{
    VersionID v = readHeader(ds, r_conneditor);

    if (v == 1)
    {
        ds >> static_cast<ConnectivityBase &>(conn);
    }
    else
        throw version_error(v, "1", r_conneditor, CODELOC);

    return ds;
}

/** Null constructor */
ConnectivityEditor::ConnectivityEditor() : ConcreteProperty<ConnectivityEditor, ConnectivityBase>()
{
}

/** Construct an editor to edit a copy of the passed
    Connectivity object */
ConnectivityEditor::ConnectivityEditor(const Connectivity &connectivity)
    : ConcreteProperty<ConnectivityEditor, ConnectivityBase>(connectivity)
{
}

/** Copy constructor */
ConnectivityEditor::ConnectivityEditor(const ConnectivityEditor &other)
    : ConcreteProperty<ConnectivityEditor, ConnectivityBase>(other)
{
}

/** Destructor */
ConnectivityEditor::~ConnectivityEditor()
{
}

/** Copy assignment operator */
ConnectivityEditor &ConnectivityEditor::operator=(const ConnectivityBase &other)
{
    ConnectivityBase::operator=(other);

    return *this;
}

/** Comparison operator */
bool ConnectivityEditor::operator==(const ConnectivityEditor &other) const
{
    return ConnectivityBase::operator==(other);
}

/** Comparison operator */
bool ConnectivityEditor::operator!=(const ConnectivityEditor &other) const
{
    return ConnectivityBase::operator!=(other);
}

/** Record the connection between the atoms at indicies 'atom0'
    and 'atom1'

    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::connect(AtomIdx atom0, AtomIdx atom1)
{
    AtomIdx atomidx0 = AtomIdx(atom0.map(connected_atoms.count()));
    AtomIdx atomidx1 = AtomIdx(atom1.map(connected_atoms.count()));

    if (atomidx0 == atomidx1)
        return *this;

    QSet<AtomIdx> *connected_atoms_array = connected_atoms.data();

    connected_atoms_array[atomidx0].insert(atomidx1);
    connected_atoms_array[atomidx1].insert(atomidx0);

    if (info().isWithinResidue(atomidx0) and info().isWithinResidue(atomidx1))
    {
        QSet<ResIdx> *connected_res_array = connected_res.data();

        ResIdx residx0 = info().parentResidue(atomidx0);
        ResIdx residx1 = info().parentResidue(atomidx1);

        if (residx0 != residx1)
        {
            connected_res_array[residx0].insert(residx1);
            connected_res_array[residx1].insert(residx0);
        }
    }

    return *this;
}

/** Record a connection between the atom identified by 'atom0' and
    the atom identified by 'atom1'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::connect(const AtomID &atom0, const AtomID &atom1)
{
    return this->connect(info().atomIdx(atom0), info().atomIdx(atom1));
}

/** Create a connection for the passed bond */
ConnectivityEditor &ConnectivityEditor::connect(const BondID &bond)
{
    return this->connect(bond.atom0(), bond.atom1());
}

/** Create a connection for the passed bonds */
ConnectivityEditor &ConnectivityEditor::connect(const QList<BondID> &bonds)
{
    for (const auto &bond : bonds)
    {
        this->connect(bond);
    }

    return *this;
}

/** Remove the connection between the atoms at indicies 'atom0'
    and 'atom1' - this does nothing if there isn't already a connection

    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::disconnect(AtomIdx atom0, AtomIdx atom1)
{
    if (this->areConnected(atom0, atom1))
    {
        AtomIdx atomidx0 = AtomIdx(atom0.map(connected_atoms.count()));
        AtomIdx atomidx1 = AtomIdx(atom1.map(connected_atoms.count()));

        connected_atoms[atomidx0].remove(atomidx1);
        connected_atoms[atomidx1].remove(atomidx0);

        // remove any properties associated with this bond
        this->bond_props.remove(SireMol::detail::IDPair(atomidx0, atomidx1));

        // now check to see if the residues are still connected
        ResIdx residx0 = info().parentResidue(atomidx0);
        ResIdx residx1 = info().parentResidue(atomidx1);

        if (this->nConnections(residx0, residx1) == 0)
        {
            connected_res[residx0].remove(residx1);
            connected_res[residx1].remove(residx0);
        }
    }

    return *this;
}

/** Disconnect the atoms that are identified by 'atom0' and 'atom1' -
    this does nothing if there isn't a connection between these atoms

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::disconnect(const AtomID &atom0, const AtomID &atom1)
{
    return this->disconnect(info().atomIdx(atom0), info().atomIdx(atom1));
}

/** Disconnect the atoms in the passed bond - this does nothing if the
 *  atoms aren't connected */
ConnectivityEditor &ConnectivityEditor::disconnect(const BondID &bond)
{
    return this->disconnect(bond.atom0(), bond.atom1());
}

/** Disconnect the atoms in the passed bonds - this does nothing for any
 *  of the atoms that aren't connected */
ConnectivityEditor &ConnectivityEditor::disconnect(const QList<BondID> &bonds)
{
    for (const auto &bond : bonds)
    {
        this->disconnect(bond);
    }

    return *this;
}

/** Disconnect any and all bonds involving the passed atoms. If exclusive is true,
 *  then this only removes connection where both atoms are in 'atoms', otherwise
 *  it removes connections which have one of more atoms in 'atoms'
 */
ConnectivityEditor &ConnectivityEditor::disconnect(const QList<AtomIdx> &atoms, bool exclusive)
{
    const QSet<AtomIdx> atoms_set(atoms.begin(), atoms.end());

    if (exclusive)
    {
        const auto bonds = this->getBonds();

        for (const auto &bond : bonds)
        {
            const auto atom0 = info().atomIdx(bond.atom0());
            const auto atom1 = info().atomIdx(bond.atom1());

            if (atoms_set.contains(atom0) and atoms_set.contains(atom1))
            {
                this->disconnect(atom0, atom1);
            }
        }
    }
    else
    {
        const auto bonds = this->getBonds();

        for (const auto &bond : bonds)
        {
            const auto atom0 = info().atomIdx(bond.atom0());
            const auto atom1 = info().atomIdx(bond.atom1());

            if (atoms_set.contains(atom0) or atoms_set.contains(atom1))
            {
                this->disconnect(atom0, atom1);
            }
        }
    }

    return *this;
}

/** Remove all of the connections to the atom at index 'atomidx'

    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::disconnectAll(AtomIdx atomidx)
{
    QSet<AtomIdx> connected = this->connectionsTo(atomidx);

    foreach (AtomIdx atom1, connected)
    {
        this->disconnect(atomidx, atom1);
    }

    return *this;
}

/** Remove all bonds from this molecule */
ConnectivityEditor &ConnectivityEditor::disconnectAll()
{
    bond_props.clear();
    connected_atoms.clear();
    connected_res.clear();

    if (info().nAtoms() > 0)
    {
        connected_atoms.resize(info().nAtoms());
        connected_atoms.squeeze();
    }

    if (info().nResidues() > 0)
    {
        connected_res.resize(info().nResidues());
        connected_res.squeeze();
    }

    return *this;
}

/** Remove all of the connections to the atom identified by 'atomid'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::disconnectAll(const AtomID &atomid)
{
    return this->disconnectAll(info().atomIdx(atomid));
}

/** Remove all of the connections that involve any of the atoms
    in the residue at index 'residx'

    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::disconnectAll(ResIdx residx)
{
    foreach (AtomIdx atomidx, info().getAtomsIn(residx))
    {
        this->disconnectAll(atomidx);
    }

    return *this;
}

/** Remove all of the connections that involve any of the atoms
    in the residue identified by 'resid'

    \throw SireMol::missing_residue
    \throw SireMol::duplicate_residue
    \throw SireError::invalid_index
*/
ConnectivityEditor &ConnectivityEditor::disconnectAll(const ResID &resid)
{
    return this->disconnectAll(info().resIdx(resid));
}

SireMol::detail::IDPair _get_id(const BondID &bond, const MoleculeInfo &minfo)
{
    return SireMol::detail::IDPair(minfo.atomIdx(bond.atom0()), minfo.atomIdx(bond.atom1()));
}

/** Set the property for the specified bond, at the specified key, to 'value' */
ConnectivityEditor &ConnectivityEditor::setProperty(const BondID &bond, const QString &key, const Property &value)
{
    auto atom0 = this->minfo.atomIdx(bond.atom0());
    auto atom1 = this->minfo.atomIdx(bond.atom1());

    if (not this->areConnected(atom0, atom1))
    {
        throw SireMol::missing_bond(QObject::tr("You cannot set the property %1 as the atoms in %2 "
                                                "are not connected.")
                                        .arg(key)
                                        .arg(bond.toString()),
                                    CODELOC);
    }

    auto id = SireMol::detail::IDPair(atom0, atom1);

    if (not this->bond_props.contains(id))
    {
        this->bond_props.insert(id, Properties());
    }

    this->bond_props[id].setProperty(key, value);

    return *this;
}

/** Set the property for the specified angle, at the specified key, to 'value' */
ConnectivityEditor &ConnectivityEditor::setProperty(const AngleID &angle, const QString &key, const Property &value)
{
    const auto id = _to_canonical(angle, this->minfo);

    if (not((this->areConnected(id.atom0(), id.atom1())) and (this->areConnected(id.atom1(), id.atom2()))))
    {
        throw SireMol::missing_angle(QObject::tr("You cannot set the property %1 as the atoms in %2 "
                                                 "are not connected.")
                                         .arg(key)
                                         .arg(angle.toString()),
                                     CODELOC);
    }

    if (not this->ang_props.contains(id))
    {
        this->ang_props.insert(id, Properties());
    }

    this->ang_props[id].setProperty(key, value);

    return *this;
}

/** Set the property for the specified dihedral, at the specified key, to 'value' */
ConnectivityEditor &ConnectivityEditor::setProperty(const DihedralID &dihedral, const QString &key,
                                                    const Property &value)
{
    const auto id = _to_canonical(dihedral, this->minfo);

    if (not((this->areConnected(id.atom0(), id.atom1())) and (this->areConnected(id.atom1(), id.atom2())) and
            (this->areConnected(id.atom2(), id.atom3()))))
    {
        throw SireMol::missing_angle(QObject::tr("You cannot set the property %1 as the atoms in %2 "
                                                 "are not connected.")
                                         .arg(key)
                                         .arg(dihedral.toString()),
                                     CODELOC);
    }

    if (not this->dih_props.contains(id))
    {
        this->dih_props.insert(id, Properties());
    }

    this->dih_props[id].setProperty(key, value);

    return *this;
}

/** Set the property for the specified improper, at the specified key, to 'value' */
ConnectivityEditor &ConnectivityEditor::setProperty(const ImproperID &improper, const QString &key,
                                                    const Property &value)
{
    const auto id = _to_canonical(improper, this->minfo);

    if (not this->imp_props.contains(id))
    {
        this->imp_props.insert(id, Properties());
    }

    this->imp_props[id].setProperty(key, value);

    return *this;
}

/** Remove the specified property from all bonds */
ConnectivityEditor &ConnectivityEditor::removeProperty(const QString &key)
{
    for (const auto &id : this->bond_props.keys())
    {
        this->bond_props[id].removeProperty(key);
    }

    for (const auto &id : this->ang_props.keys())
    {
        this->ang_props[id].removeProperty(key);
    }

    for (const auto &id : this->dih_props.keys())
    {
        this->dih_props[id].removeProperty(key);
    }

    for (const auto &id : this->imp_props.keys())
    {
        this->imp_props[id].removeProperty(key);
    }

    return *this;
}

/** Remove the specified property from the specified bond */
ConnectivityEditor &ConnectivityEditor::removeProperty(const BondID &bond, const QString &key)
{
    auto id = _get_id(bond, this->minfo);

    if (this->bond_props.contains(id))
    {
        this->bond_props[id].removeProperty(key);
    }

    return *this;
}

/** Remove the specified property from the specified angle */
ConnectivityEditor &ConnectivityEditor::removeProperty(const AngleID &angle, const QString &key)
{
    auto id = _to_canonical(angle, this->minfo);

    if (this->ang_props.contains(id))
    {
        this->ang_props[id].removeProperty(key);
    }

    return *this;
}

/** Remove the specified property from the specified dihedral */
ConnectivityEditor &ConnectivityEditor::removeProperty(const DihedralID &dihedral, const QString &key)
{
    auto id = _to_canonical(dihedral, this->minfo);

    if (this->dih_props.contains(id))
    {
        this->dih_props[id].removeProperty(key);
    }

    return *this;
}

/** Remove the specified property from the specified improper */
ConnectivityEditor &ConnectivityEditor::removeProperty(const ImproperID &imp, const QString &key)
{
    auto id = _to_canonical(imp, this->minfo);

    if (this->imp_props.contains(id))
    {
        this->imp_props[id].removeProperty(key);
    }

    return *this;
}

/** Take the specified property from the specified bond - this removes
    and returns the property if it exists. If it doesn't, then
    a NullProperty is returned
*/
PropertyPtr ConnectivityEditor::takeProperty(const BondID &bond, const QString &key)
{
    auto id = _get_id(bond, this->minfo);

    if (this->bond_props.contains(id))
    {
        PropertyPtr value = this->bond_props[id].property(key);
        this->bond_props[id].removeProperty(key);

        return value;
    }
    else
    {
        return NullProperty();
    }
}

/** Take the specified property from the specified angle - this removes
    and returns the property if it exists. If it doesn't, then
    a NullProperty is returned
*/
PropertyPtr ConnectivityEditor::takeProperty(const AngleID &angle, const QString &key)
{
    auto id = _to_canonical(angle, this->minfo);

    if (this->ang_props.contains(id))
    {
        PropertyPtr value = this->ang_props[id].property(key);
        this->ang_props[id].removeProperty(key);

        return value;
    }
    else
    {
        return NullProperty();
    }
}

/** Take the specified property from the specified dihedral - this removes
    and returns the property if it exists. If it doesn't, then
    a NullProperty is returned
*/
PropertyPtr ConnectivityEditor::takeProperty(const DihedralID &dihedral, const QString &key)
{
    auto id = _to_canonical(dihedral, this->minfo);

    if (this->dih_props.contains(id))
    {
        PropertyPtr value = this->dih_props[id].property(key);
        this->dih_props[id].removeProperty(key);

        return value;
    }
    else
    {
        return NullProperty();
    }
}

/** Take the specified property from the specified improper - this removes
    and returns the property if it exists. If it doesn't, then
    a NullProperty is returned
*/
PropertyPtr ConnectivityEditor::takeProperty(const ImproperID &improper, const QString &key)
{
    auto id = _to_canonical(improper, this->minfo);

    if (this->imp_props.contains(id))
    {
        PropertyPtr value = this->imp_props[id].property(key);
        this->imp_props[id].removeProperty(key);

        return value;
    }
    else
    {
        return NullProperty();
    }
}

/** Return the editied connectivity */
Connectivity ConnectivityEditor::commit() const
{
    return Connectivity(*this);
}

const char *ConnectivityEditor::typeName()
{
    return QMetaType::typeName(qMetaTypeId<ConnectivityEditor>());
}
