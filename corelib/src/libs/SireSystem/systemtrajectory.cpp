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

#include "SireBase/lazyevaluator.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

namespace SireSystem
{
    class SystemFrames
    {
    public:
        SystemFrames();
        ~SystemFrames();

        int nFrames() const;
        int nAtoms() const;

        Frame getFrame(int i);
        Frame getFrame(int i, const LazyEvaluator &evaluator);

        void saveFrame(const Molecules &mols, const PropertyMap &map);
    };
}

SystemFrames::SystemFrames()
{
}

SystemFrames::~SystemFrames()
{
}

int SystemFrames::nAtoms() const
{
    return 0;
}

int SystemFrames::nFrames() const
{
    return 0;
}

Frame SystemFrames::getFrame(int i)
{
    return Frame();
}

Frame SystemFrames::getFrame(int i, const LazyEvaluator &evaluator)
{
    return Frame();
}

void SystemFrames::saveFrame(const Molecules &mols, const PropertyMap &map)
{
}

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

    return d->getFrame(i).subset(start_atom, natoms);
}

Frame MolSystemTrajectory::getFrame(int i, const SireBase::LazyEvaluator &evaluator) const
{
    i = SireID::Index(i).map(this->nFrames());

    auto key = QString("%1-%2").arg(qintptr(d.get())).arg(i);

    auto frame = evaluator.evaluate(key, [&]()
                                    { return d->getFrame(i); });

    return frame.read().asA<Frame>().subset(start_atom, natoms);
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

static const RegisterMetaType<SystemTrajectory> r_traj;

SIRESYSTEM_EXPORT QDataStream &operator<<(QDataStream &ds, const SystemTrajectory &traj)
{
    writeHeader(ds, r_traj, 1);

    // we don't stream the trajectory as it would be too big
    ds << static_cast<const TrajectoryData &>(traj);

    return ds;
}

SIRESYSTEM_EXPORT QDataStream &operator>>(QDataStream &ds, SystemTrajectory &traj)
{
    auto v = readHeader(ds, r_traj);

    if (v == 1)
    {
        // we don't stream the trajectory as it would be too big
        traj.clear();
        ds >> static_cast<TrajectoryData &>(traj);
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

SystemTrajectory::SystemTrajectory() : TrajectoryData()
{
}

SystemTrajectory::SystemTrajectory(const Molecules &mols, const PropertyMap &map)
    : TrajectoryData()
{
}

SystemTrajectory::SystemTrajectory(const SystemTrajectory &other)
    : TrajectoryData(other), d(other.d), mol_atoms(other.mol_atoms)
{
}

SystemTrajectory::~SystemTrajectory()
{
}

SystemTrajectory &SystemTrajectory::operator=(const SystemTrajectory &other)
{
    if (this != &other)
    {
        TrajectoryData::operator=(other);
        d = other.d;
        mol_atoms = other.mol_atoms;
    }

    return *this;
}

bool SystemTrajectory::operator==(const SystemTrajectory &other) const
{
    return TrajectoryData::operator==(other) &&
           d.get() == other.d.get() &&
           mol_atoms == other.mol_atoms;
}

bool SystemTrajectory::operator!=(const SystemTrajectory &other) const
{
    return not this->operator==(other);
}

const char *SystemTrajectory::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SystemTrajectory>());
}

const char *SystemTrajectory::what() const
{
    return SystemTrajectory::typeName();
}

SystemTrajectory *SystemTrajectory::clone() const
{
    return new SystemTrajectory(*this);
}

void SystemTrajectory::clear()
{
    d.reset();
    mol_atoms.clear();
}

bool SystemTrajectory::isCompatibleWith(const Molecules &mols,
                                        const PropertyMap &map) const
{
    if (d.get() == 0)
        return false;

    // make sure that all of the molecules exist in the hash
    // and the number of atoms match
    for (const auto &mol : mols)
    {
        const auto &moldata = mol.data();

        auto it = mol_atoms.constFind(moldata.number());

        if (it == mol_atoms.constEnd())
            return false;

        if (it.value().second != moldata.info().nAtoms())
            return false;
    }

    return true;
}

void SystemTrajectory::saveFrame(const Molecules &mols,
                                 const PropertyMap &map)
{
    if (d.get() == 0)
    {
        // create the data
        int natoms = 0;

        for (const auto &mol : mols)
        {
            const auto &moldata = mol.data();

            mol_atoms.insert(moldata.number(), qMakePair(natoms, moldata.info().nAtoms()));

            natoms += moldata.info().nAtoms();
        }

        d.reset(new SystemFrames());
    }

    // save the frame
    d->saveFrame(mols, map);
}

TrajectoryDataPtr SystemTrajectory::getTrajectory(MolNum molnum) const
{
    return TrajectoryDataPtr(new MolSystemTrajectory(*this, molnum));
}

int SystemTrajectory::nFrames() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->nFrames();
}

int SystemTrajectory::nAtoms() const
{
    if (d.get() == 0)
    {
        int natoms = 0;

        for (const auto &mol : mol_atoms)
            natoms += mol.second;

        return natoms;
    }
    else
    {
        return d->nAtoms();
    }
}

QStringList SystemTrajectory::filenames() const
{
    return QStringList();
}

Frame SystemTrajectory::getFrame(int i) const
{
    i = SireID::Index(i).map(this->nFrames());
    return d->getFrame(i);
}

Frame SystemTrajectory::getFrame(int i, const LazyEvaluator &evaluator) const
{
    i = SireID::Index(i).map(this->nFrames());
    return d->getFrame(i, evaluator);
}

bool SystemTrajectory::isEditable() const
{
    return false;
}

bool SystemTrajectory::_equals(const TrajectoryData &other) const
{
    const SystemTrajectory *p = dynamic_cast<const SystemTrajectory *>(&other);

    if (p)
        return this->operator==(*p);
    else
        return false;
}
