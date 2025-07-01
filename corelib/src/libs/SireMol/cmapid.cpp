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

#include "cmapid.h"

#include "moleculeinfodata.h"

#include "SireID/index.h"

#include "SireBase/property.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireID;
using namespace SireStream;

static const RegisterMetaType<CMAPID> r_cmapid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const CMAPID &cmapid)
{
    writeHeader(ds, r_cmapid, 1);

    SharedDataStream sds(ds);

    sds << cmapid.atm0 << cmapid.atm1 << cmapid.atm2 << cmapid.atm3 << cmapid.atm4;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, CMAPID &cmapid)
{
    VersionID v = readHeader(ds, r_cmapid);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> cmapid.atm0 >> cmapid.atm1 >> cmapid.atm2 >> cmapid.atm3 >> cmapid.atm4;
    }
    else
        throw version_error(v, "1", r_cmapid, CODELOC);

    return ds;
}

/** Null constructor */
CMAPID::CMAPID() : ID()
{
}

/** Construct a cmap between the two specified atoms. The order
    is important, as this cmap may be between two different
    molecules */
CMAPID::CMAPID(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2,
               const AtomID &atom3, const AtomID &atom4)
    : ID(), atm0(atom0), atm1(atom1), atm2(atom2), atm3(atom3), atm4(atom4)
{
}

/** Copy constructor */
CMAPID::CMAPID(const CMAPID &other)
    : ID(other), atm0(other.atm0), atm1(other.atm1), atm2(other.atm2),
      atm3(other.atm3), atm4(other.atm4)
{
}

/** Destructor */
CMAPID::~CMAPID()
{
}

/** Copy assignment operator */
CMAPID &CMAPID::operator=(const CMAPID &other)
{
    atm0 = other.atm0;
    atm1 = other.atm1;
    atm2 = other.atm2;
    atm3 = other.atm3;
    atm4 = other.atm4;

    return *this;
}

/** Comparison operator - the order is important */
bool CMAPID::operator==(const CMAPID &other) const
{
    return atm0 == other.atm0 and atm1 == other.atm1 and
           atm2 == other.atm2 and atm3 == other.atm3 and
           atm4 == other.atm4;
}

/** Comparison operator - the order is important */
bool CMAPID::operator!=(const CMAPID &other) const
{
    return atm0 != other.atm0 or atm1 != other.atm1 or
           atm2 != other.atm2 or atm3 != other.atm3 or
           atm4 != other.atm4;
}

const AtomID &CMAPID::operator[](int i) const
{
    i = Index(i).map(5);

    switch (i)
    {
    case 0:
        return atm0.base();
    case 1:
        return atm1.base();
    case 2:
        return atm2.base();
    case 3:
        return atm3.base();
    default:
        return atm4.base();
    }
}

/** Return the mirror of this CMAPID - i.e. if this is
    CMAPID(atom0, atom1, atom2, atom3, atom4), this returns
    CMAPID(atom4, atom3, atom2, atom1, atom0).

    This is useful if you know that CMAPID(atom0,atom1,atom2,atom3,atom4) equals
    CMAPID(atom4,atom3,atom2,atom1,atom0), e.g. you can now write;

    if (not (cmaps.contains(cmap) or cmaps.contains(cmap.mirror())) )
    {
        cmaps.insert(cmap);
    }

    or

    if (cmap == other_cmap or cmap == other_cmap.mirror())
    {
        //this is the same cmap
    }
*/
CMAPID CMAPID::mirror() const
{
    return CMAPID(atm4, atm3, atm2, atm1, atm0);
}

/** Return a hash for this ID */
uint CMAPID::hash() const
{
    return ((atm0.hash() * atm1.hash()) << 16) | ((atm2.hash() * atm3.hash() * atm4.hash()) & 0x0000FFFF);
}

/** Return a string representation of this ID */
QString CMAPID::toString() const
{
    return QString("CMAP( %1, %2, %3, %4, %5 )")
        .arg(atm0.toString(), atm1.toString(), atm2.toString(), atm3.toString(), atm4.toString());
}

/** Return whether this is a null ID */
bool CMAPID::isNull() const
{
    return atm0.isNull() and atm1.isNull() and atm2.isNull() and atm3.isNull() and atm4.isNull();
}

/** Comparison operator with another ID */
bool CMAPID::operator==(const SireID::ID &other) const
{
    const CMAPID *other_cmap = dynamic_cast<const CMAPID *>(&other);

    return other_cmap and this->operator==(*other_cmap);
}

/** Return the indicies of the five atoms in this cmap - this returns
    them in the order
    tuple(cmap.atom0(),cmap.atom1(),cmap.atom2(),cmap.atom3(),cmap.atom4())

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx, AtomIdx, AtomIdx, AtomIdx, AtomIdx> CMAPID::map(const MoleculeInfoData &molinfo) const
{
    return tuple<AtomIdx, AtomIdx, AtomIdx, AtomIdx, AtomIdx>(molinfo.atomIdx(atm0), molinfo.atomIdx(atm1),
                                                              molinfo.atomIdx(atm2), molinfo.atomIdx(atm3),
                                                              molinfo.atomIdx(atm4));
}

/** Return the indicies of the five atoms of this cmap, between the
    two molecules whose data is in 'mol0info' (containing cmap.atom0()),
    'mol1info' (containing cmap.atom1()), 'mol2info' (containing
    cmap.atom2()), 'mol3info' (containing cmap.atom3()) and
    'mol4info' (containing cmap.atom4())

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
tuple<AtomIdx, AtomIdx, AtomIdx, AtomIdx, AtomIdx> CMAPID::map(const MoleculeInfoData &mol0info,
                                                               const MoleculeInfoData &mol1info,
                                                               const MoleculeInfoData &mol2info,
                                                               const MoleculeInfoData &mol3info,
                                                               const MoleculeInfoData &mol4info) const
{
    return tuple<AtomIdx, AtomIdx, AtomIdx, AtomIdx, AtomIdx>(mol0info.atomIdx(atm0), mol1info.atomIdx(atm1),
                                                              mol2info.atomIdx(atm2), mol3info.atomIdx(atm3),
                                                              mol4info.atomIdx(atm4));
}

/** Return the ID of the first atom of the cmap */
const AtomID &CMAPID::atom0() const
{
    return atm0.base();
}

/** Return the ID of the second atom of the cmap */
const AtomID &CMAPID::atom1() const
{
    return atm1.base();
}

/** Return the ID of the third atom of the cmap */
const AtomID &CMAPID::atom2() const
{
    return atm2.base();
}

/** Return the ID of the fourth atom of the cmap */
const AtomID &CMAPID::atom3() const
{
    return atm3.base();
}

/** Return the ID of the fifth atom of the cmap */
const AtomID &CMAPID::atom4() const
{
    return atm4.base();
}

const char *CMAPID::typeName()
{
    return QMetaType::typeName(qMetaTypeId<CMAPID>());
}

CMAPID *CMAPID::clone() const
{
    return new CMAPID(*this);
}
