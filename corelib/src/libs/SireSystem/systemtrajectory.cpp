/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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

#include "systemtrajectory.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

////////
//////// Implementation of MolSystemTrajectory
////////

static const RegisterMetaType<MolSystemTrajectory> r_moltraj;

SIRESYSTEM_EXPORT QDataStream &operator<<(QDataStream &ds, const MolSystemTrajectory &traj)
{
    writeHeader(ds, r_moltraj, 1);

    // we don't stream the trajectory as it would be too big
    ds << static_cast<const TrajectoryData &>(traj);

    return ds;
}

SIRESYSTEM_EXPORT QDataStream &operator>>(QDataStream &ds, MolSystemTrajectory &traj)
{
    auto v = readHeader(ds, r_moltraj);

    if (v == 1)
    {
        // we don't stream the trajectory as it would be too big
        traj.clear();
        ds >> static_cast<TrajectoryData &>(traj);
    }
    else
        throw version_error(v, "1", r_moltraj, CODELOC);

    return ds;
}

MolSystemTrajectory::MolSystemTrajectory()
    : TrajectoryData(), start_atom(0), natoms(0)
{
}

MolSystemTrajectory::MolSystemTrajectory(const SystemTrajectory &trajectory,
                                         SireMol::MolNum molnum)
    : TrajectoryData(trajectory), start_atom(0), natoms(0)
{
    d = trajectory.d;

    auto it = trajectory.mol_atoms.constFind(molnum);

    if (it == trajectory.mol_atoms.constEnd())
    {
        throw SireMol::missing_molecule(QObject::tr(
                                            "There is no molecule with number %1 in the system")
                                            .arg(molnum.value()),
                                        CODELOC);
    }

    start_atom = it.value().first;
    natoms = it.value().second;
}

MolSystemTrajectory::MolSystemTrajectory(const MolSystemTrajectory &other)
    : TrajectoryData(other), d(other.d), start_atom(other.start_atom), natoms(other.natoms)
{
}

MolSystemTrajectory::~MolSystemTrajectory()
{
}

MolSystemTrajectory &MolSystemTrajectory::operator=(const MolSystemTrajectory &other)
{
    if (this != &other)
    {
        TrajectoryData::operator=(other);
        d = other.d;
        start_atom = other.start_atom;
        natoms = other.natoms;
    }

    return *this;
}

bool MolSystemTrajectory::operator==(const MolSystemTrajectory &other) const
{
    return TrajectoryData::operator==(other) &&
           d.get() == other.d.get() &&
           start_atom == other.start_atom &&
           natoms == other.natoms;
}

bool MolSystemTrajectory::operator!=(const MolSystemTrajectory &other) const
{
    return not this->operator==(other);
}

const char *MolSystemTrajectory::typeName()
{
    return QMetaType::typeName(qMetaTypeId<MolSystemTrajectory>());
}

const char *MolSystemTrajectory::what() const
{
    return MolSystemTrajectory::typeName();
}

MolSystemTrajectory *MolSystemTrajectory::clone() const
{
    return new MolSystemTrajectory(*this);
}

void MolSystemTrajectory::clear()
{
    d.reset();
    start_atom = 0;
    natoms = 0;
}

int MolSystemTrajectory::nFrames() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->nFrames();
}

int MolSystemTrajectory::nAtoms() const
{
    return natoms;
}

QStringList MolSystemTrajectory::filenames() const
{
    return QStringList();
}

Frame MolSystemTrajectory::getFrame(int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    return d->getFrame(i, start_atom, natoms);
}

Frame MolSystemTrajectory::getFrame(int i, const SireBase::LazyEvaluator &evaluator) const
{
    i = SireID::Index(i).map(this->nFrames());

    return d->getFrame(i, start_atom, natoms, evaluator);
}

bool MolSystemTrajectory::isEditable() const
{
    return false;
}

bool MolSystemTrajectory::_equals(const TrajectoryData &other) const
{
    const MolSystemTrajectory *p = dynamic_cast<const MolSystemTrajectory *>(&other);

    if (p)
        return this->operator==(*p);
    else
        return false;
}

////////
//////// Implementation of SystemTrajectory
////////
