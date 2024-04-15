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
#include "SireBase/parallel.h"
#include "SireBase/pagecache.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <tbb/spin_mutex.h>

#include <QDir>

using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

namespace SireSystem
{
    std::weak_ptr<PageCache> shared_cache;

    /** Return a pointer to the PageCache that is shared by
     *  all SystemFrames objects
     */
    std::shared_ptr<PageCache> getSharedCache()
    {
        auto cache = shared_cache.lock();

        static QMutex mutex;

        if (not cache)
        {
            QMutexLocker locker(&mutex);
            cache = shared_cache.lock();

            if (not cache)
            {
                QString cache_dir = QDir::current().absoluteFilePath("trajectory_cache_XXXXXX");
                // use 32 MB pages
                cache = std::make_shared<PageCache>(cache_dir, 32 * 1024 * 1024);
                shared_cache = cache;
            }
        }

        return cache;
    }

    class SystemFrames
    {
    public:
        enum FrameType
        {
            EMPTY = 0x0000,
            COORDINATES = 0x0001,
            VELOCITIES = 0x0010,
            FORCES = 0x0100
        };

        SystemFrames();
        SystemFrames(const QList<MolNum> &molnums,
                     const Molecules &mols, const PropertyMap &map);
        ~SystemFrames();

        int nFrames() const;

        int nAtoms() const;
        int nAtoms(MolNum molnum) const;

        bool saveCoordinates() const;
        bool saveVelocities() const;
        bool saveForces() const;

        Frame getFrame(int i) const;
        Frame getFrame(int i, const LazyEvaluator &evaluator) const;

        Frame getFrame(MolNum molnum, int i) const;
        Frame getFrame(MolNum molnum, int i,
                       const LazyEvaluator &evaluator) const;

        void saveFrame(const Molecules &mols,
                       const Space &space,
                       SireUnits::Dimension::Time time,
                       const Properties &props,
                       const PropertyMap &map);

        bool isCompatibleWith(const QList<MolNum> &molnums,
                              const Molecules &mols,
                              const PropertyMap &map) const;

    private:
        tbb::spin_mutex &getMutex() const;
        Frame _lkr_getFrame(int i) const;

        /** The order the molecules should appear in the trajectory */
        QList<MolNum> molnums;

        /** The start index and number of atoms for each molecule
         *  in the system. This is the same for all frames
         */
        QHash<SireMol::MolNum, QPair<int, int>> mol_atoms;

        /** The data for the trajectory */
        QVector<PageCache::Handle> frames;

        /** Pointer to the shared pagecache */
        std::shared_ptr<PageCache> cache;

        /** Mutex to serialize access to the data of this class */
        tbb::spin_mutex mutex;

        /** The current (unpacked) frame */
        Frame current_frame;

        /** The index of the current frame */
        int current_frame_index;

        /** The total number of atoms */
        int natoms;

        /** What type of frame data - coordinates, velocities, forces */
        int frame_type;
    };
}

tbb::spin_mutex &SystemFrames::getMutex() const
{
    return *(const_cast<tbb::spin_mutex *>(&mutex));
}

SystemFrames::SystemFrames()
    : current_frame_index(-1), natoms(0), frame_type(EMPTY)
{
    cache = getSharedCache();
}

SystemFrames::SystemFrames(const QList<MolNum> &nums,
                           const Molecules &mols, const PropertyMap &map)
    : molnums(nums), current_frame_index(-1), natoms(0), frame_type(EMPTY)
{
    cache = getSharedCache();

    tbb::spin_mutex::scoped_lock lock(getMutex());

    for (const auto &molnum : molnums)
    {
        auto it = mols.constFind(molnum);

        if (it == mols.constEnd())
            throw SireMol::missing_molecule(QObject::tr(
                                                "There is no molecule with number %1 in the system")
                                                .arg(molnum.value()),
                                            CODELOC);

        const auto &moldata = it.value().data();

        mol_atoms.insert(molnum, qMakePair(natoms, moldata.info().nAtoms()));

        natoms += moldata.info().nAtoms();
    }

    // should work out if coords, vels and/or forces should be saved
    // based on the properties in the map
    frame_type = COORDINATES;

    bool save_coordinates = true;
    bool save_velocities = false;
    bool save_forces = false;

    if (map.specified("save_coordinates"))
    {
        save_coordinates = map["save_coordinates"].value().asABoolean();
    }

    if (map.specified("save_velocities"))
    {
        save_velocities = map["save_velocities"].value().asABoolean();
    }

    if (map.specified("save_forces"))
    {
        save_forces = map["save_forces"].value().asABoolean();
    }

    if (save_coordinates)
        frame_type |= COORDINATES;

    if (save_velocities)
        frame_type |= VELOCITIES;

    if (save_forces)
        frame_type |= FORCES;
}

SystemFrames::~SystemFrames()
{
}

bool SystemFrames::saveCoordinates() const
{
    return frame_type & COORDINATES;
}

bool SystemFrames::saveVelocities() const
{
    return frame_type & VELOCITIES;
}

bool SystemFrames::saveForces() const
{
    return frame_type & FORCES;
}

int SystemFrames::nAtoms() const
{
    return natoms;
}

int SystemFrames::nAtoms(MolNum molnum) const
{
    tbb::spin_mutex::scoped_lock lock(getMutex());

    auto it = mol_atoms.constFind(molnum);

    if (it == mol_atoms.constEnd())
    {
        throw SireMol::missing_molecule(QObject::tr(
                                            "There is no molecule with number %1 in the system")
                                            .arg(molnum.value()),
                                        CODELOC);
    }

    return it.value().second;
}

int SystemFrames::nFrames() const
{
    tbb::spin_mutex::scoped_lock lock(getMutex());
    return frames.count();
}

/** It is only compatible if we have the same molecules with
 *  the same number of atoms. This is because we cannot cope
 *  with molecules being added or removed from the system,
 *  or with the number of atoms changing. These events will
 *  trigger the creation of a new SystemFrames higher in the
 *  stack
 */
bool SystemFrames::isCompatibleWith(const QList<MolNum> &nums,
                                    const Molecules &mols,
                                    const PropertyMap &map) const
{
    tbb::spin_mutex::scoped_lock lock(getMutex());

    if (molnums != nums or mols.nMolecules() != mol_atoms.size())
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

Frame SystemFrames::_lkr_getFrame(int i) const
{
    try
    {
        i = Index(i).map(frames.count());
    }
    catch (...)
    {
        throw SireError::invalid_index(
            QObject::tr("Invalid frame index %1. Number of frames is %2.")
                .arg(i)
                .arg(frames.count()),
            CODELOC);
    }

    if (i == current_frame_index)
    {
        return current_frame;
    }

    auto data = frames.at(i).fetch();

    QDataStream ds(data);

    Frame frame;

    // auto start_time = std::chrono::high_resolution_clock::now();
    ds >> frame;
    // auto end_time = std::chrono::high_resolution_clock::now();

    // qDebug() << "Loading frame" << i << "with" << frames.at(i).size() << "bytes" << frame.numBytes();
    // qDebug() << "Deserialization time" << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "microseconds";

    const_cast<SystemFrames *>(this)->current_frame = frame;
    const_cast<SystemFrames *>(this)->current_frame_index = i;

    return current_frame;
}

Frame SystemFrames::getFrame(int i) const
{
    tbb::spin_mutex::scoped_lock lock(getMutex());
    return this->_lkr_getFrame(i);
}

Frame SystemFrames::getFrame(int i, const LazyEvaluator &evaluator) const
{
    auto key = QString("%1-%2").arg(qintptr(this)).arg(i);

    auto frame = evaluator.evaluate(key, [&]()
                                    { return this->getFrame(i); });

    return frame.read().asA<Frame>();
}

Frame SystemFrames::getFrame(MolNum molnum, int i) const
{
    tbb::spin_mutex::scoped_lock lock(getMutex());

    auto it = mol_atoms.constFind(molnum);

    if (it == mol_atoms.constEnd())
    {
        throw SireMol::missing_molecule(QObject::tr(
                                            "There is no molecule with number %1 in the system")
                                            .arg(molnum.value()),
                                        CODELOC);
    }

    return this->_lkr_getFrame(i).subset(it.value().first, it.value().second);
}

Frame SystemFrames::getFrame(MolNum molnum, int i,
                             const LazyEvaluator &evaluator) const
{
    tbb::spin_mutex::scoped_lock lock(getMutex());
    auto it = mol_atoms.constFind(molnum);

    if (it == mol_atoms.constEnd())
    {
        throw SireMol::missing_molecule(QObject::tr(
                                            "There is no molecule with number %1 in the system")
                                            .arg(molnum.value()),
                                        CODELOC);
    }

    // below function call gets the lock, so unlock now
    lock.release();

    return this->getFrame(i, evaluator).subset(it.value().first, it.value().second);
}

void SystemFrames::saveFrame(const Molecules &mols,
                             const Space &space,
                             SireUnits::Dimension::Time time,
                             const Properties &props,
                             const PropertyMap &map)
{
    const bool save_coords = this->saveCoordinates();
    const bool save_vels = this->saveVelocities();
    const bool save_forces = this->saveForces();

    if (not save_coords and not save_vels and not save_forces)
    {
        return;
    }

    QVector<Vector> coordinates;
    QVector<Velocity3D> velocities;
    QVector<Force3D> forces;

    Vector *coordinates_data = 0;
    Velocity3D *velocities_data = 0;
    Force3D *forces_data = 0;

    if (save_coords)
    {
        coordinates.resize(natoms);
        coordinates_data = coordinates.data();
    }

    if (save_vels)
    {
        velocities.resize(natoms);
        velocities_data = velocities.data();
    }

    if (save_forces)
    {
        forces.resize(natoms);
        forces_data = forces.data();
    }

    auto save_frame = [&](MolNum molnum)
    {
        auto it = mol_atoms.constFind(molnum);

        if (it == mol_atoms.constEnd())
        {
            return;
        }

        auto it2 = mols.constFind(molnum);

        if (it2 == mols.constEnd())
        {
            return;
        }

        int start_atom = it.value().first;
        int nats = it.value().second;

        if (start_atom < 0 or nats == 0 or start_atom + nats > natoms)
        {
            return;
        }

        const auto &moldata = it2.value().data();

        if (save_coords)
        {
            try
            {
                const auto &coords = moldata.property(map["coordinates"]).asA<AtomCoords>();

                if (coords.nAtoms() != nats)
                {
                    return;
                }

                std::memcpy(coordinates_data + start_atom, coords.constData(CGIdx(0)),
                            nats * sizeof(Vector));
            }
            catch (...)
            {
            }
        }

        if (save_vels)
        {
            try
            {
                const auto &vels = moldata.property(map["velocities"]).asA<AtomVelocities>();

                if (vels.nAtoms() != nats)
                {
                    return;
                }

                std::memcpy(velocities_data + start_atom, vels.constData(CGIdx(0)),
                            nats * sizeof(Velocity3D));
            }
            catch (...)
            {
            }
        }

        if (save_forces)
        {
            try
            {
                const auto &frcs = moldata.property(map["forces"]).asA<AtomForces>();

                if (frcs.nAtoms() != nats)
                {
                    return;
                }

                std::memcpy(forces_data + start_atom, frcs.constData(CGIdx(0)),
                            nats * sizeof(Force3D));
            }
            catch (...)
            {
            }
        }
    };

    if (should_run_in_parallel(molnums.count(), map))
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, molnums.count()), [&](const tbb::blocked_range<int> &r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                save_frame(molnums[i]);
            } });
    }
    else
    {
        for (const auto &molnum : molnums)
        {
            save_frame(molnum);
        }
    }

    Frame frame(coordinates, velocities, forces, space, time, props);

    auto start_time = std::chrono::high_resolution_clock::now();

    QByteArray data1 = frame.toByteArray();

    auto mem_time1 = std::chrono::high_resolution_clock::now();

    QByteArray data;
    data.reserve(2 * frame.numBytes());

    auto mem_time = std::chrono::high_resolution_clock::now();

    QDataStream ds(&data, QIODevice::WriteOnly);
    ds << frame;

    auto cache_time = std::chrono::high_resolution_clock::now();

    auto handle = cache->store(data);

    auto end_time = std::chrono::high_resolution_clock::now();

    qDebug() << "Saving frame" << frames.count() << "with" << data.size() << "bytes" << frame.numBytes();
    qDebug() << "Memory time" << std::chrono::duration_cast<std::chrono::microseconds>(mem_time - start_time).count() << "microseconds";
    qDebug() << "Memory time" << std::chrono::duration_cast<std::chrono::microseconds>(mem_time - mem_time1).count() << "microseconds";
    qDebug() << "Serialization time" << std::chrono::duration_cast<std::chrono::microseconds>(cache_time - mem_time).count() << "microseconds";
    qDebug() << "Cache time" << std::chrono::duration_cast<std::chrono::microseconds>(end_time - cache_time).count() << "microseconds";
    qDebug() << "Total time" << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() << "microseconds";

    // need to hold the write lock, as we are updating global state
    // for everyone who holds this live trajectory data
    tbb::spin_mutex::scoped_lock lock(getMutex());
    frames.append(handle);
    current_frame = frame;
    current_frame_index = frames.count() - 1;
}

////////
//////// Implementation of MolSystemTrajectory
////////

static const RegisterMetaType<MolSystemTrajectory> r_moltraj;

SIRESYSTEM_EXPORT QDataStream &operator<<(QDataStream &ds, const MolSystemTrajectory &traj)
{
    writeHeader(ds, r_moltraj, 1);

    // we don't stream the trajectory as it would be too big
    ds << traj.molnum
       << static_cast<const TrajectoryData &>(traj);

    return ds;
}

SIRESYSTEM_EXPORT QDataStream &operator>>(QDataStream &ds, MolSystemTrajectory &traj)
{
    auto v = readHeader(ds, r_moltraj);

    if (v == 1)
    {
        // we don't stream the trajectory as it would be too big
        traj.clear();
        ds >> traj.molnum >> static_cast<TrajectoryData &>(traj);
    }
    else
        throw version_error(v, "1", r_moltraj, CODELOC);

    return ds;
}

MolSystemTrajectory::MolSystemTrajectory() : TrajectoryData()
{
}

MolSystemTrajectory::MolSystemTrajectory(const SystemTrajectory &trajectory,
                                         SireMol::MolNum mnum)
    : TrajectoryData(trajectory), d(trajectory.d), molnum(mnum)
{
}

MolSystemTrajectory::MolSystemTrajectory(const MolSystemTrajectory &other)
    : TrajectoryData(other), d(other.d), molnum(other.molnum)
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
        molnum = other.molnum;
    }

    return *this;
}

bool MolSystemTrajectory::operator==(const MolSystemTrajectory &other) const
{
    return TrajectoryData::operator==(other) &&
           d.get() == other.d.get() &&
           molnum == other.molnum;
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
    molnum = MolNum();
}

int MolSystemTrajectory::nFrames() const
{
    if (not d)
        return 0;
    else
        return d->nFrames();
}

int MolSystemTrajectory::nAtoms() const
{
    if (d)
    {
        return d->nAtoms(molnum);
    }
    else
    {
        return 0;
    }
}

QStringList MolSystemTrajectory::filenames() const
{
    return QStringList();
}

Frame MolSystemTrajectory::getFrame(int i) const
{
    if (d)
        return d->getFrame(molnum, i);
    else
        throw SireError::invalid_index(
            QObject::tr("Invalid frame index %1. Number of frames is 0.")
                .arg(i),
            CODELOC);
}

Frame MolSystemTrajectory::getFrame(int i, const SireBase::LazyEvaluator &evaluator) const
{
    if (d)
        return d->getFrame(molnum, i, evaluator);
    else
        throw SireError::invalid_index(
            QObject::tr("Invalid frame index %1. Number of frames is 0.")
                .arg(i),
            CODELOC);
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

SystemTrajectory::SystemTrajectory(const QList<MolNum> &molnums,
                                   const Molecules &mols,
                                   const PropertyMap &map)
    : TrajectoryData(), d(new SystemFrames(molnums, mols, map))
{
}

SystemTrajectory::SystemTrajectory(const SystemTrajectory &other)
    : TrajectoryData(other), d(other.d)
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
    }

    return *this;
}

bool SystemTrajectory::operator==(const SystemTrajectory &other) const
{
    return TrajectoryData::operator==(other) &&
           d.get() == other.d.get();
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

bool SystemTrajectory::isLive() const
{
    return true;
}

void SystemTrajectory::clear()
{
    d.reset();
}

bool SystemTrajectory::isCompatibleWith(const QList<MolNum> &molnums,
                                        const Molecules &mols,
                                        const PropertyMap &map) const
{
    if (d)
        return d->isCompatibleWith(molnums, mols, map);
    else
        return false;
}

void SystemTrajectory::saveFrame(const Molecules &mols,
                                 const Space &space,
                                 SireUnits::Dimension::Time time,
                                 const Properties &props,
                                 const PropertyMap &map)
{
    if (not d)
    {
        throw SireError::invalid_state(
            QObject::tr("The trajectory is not initialized"),
            CODELOC);
    }

    d->saveFrame(mols, space, time, props, map);
}

TrajectoryDataPtr SystemTrajectory::getTrajectory(MolNum molnum) const
{
    if (not d)
    {
        throw SireMol::missing_molecule(QObject::tr(
                                            "There is no molecule with number %1 in the system")
                                            .arg(molnum.value()),
                                        CODELOC);
    }

    return TrajectoryDataPtr(new MolSystemTrajectory(*this, molnum));
}

int SystemTrajectory::nFrames() const
{
    if (not d)
        return 0;
    else
        return d->nFrames();
}

int SystemTrajectory::nAtoms() const
{
    if (d)
    {
        return d->nAtoms();
    }
    else
    {
        return 0;
    }
}

QStringList SystemTrajectory::filenames() const
{
    return QStringList();
}

Frame SystemTrajectory::getFrame(int i) const
{
    if (d)
        return d->getFrame(i);
    else
        throw SireError::invalid_index(
            QObject::tr("Invalid frame index %1. Number of frames is 0.")
                .arg(i),
            CODELOC);
}

Frame SystemTrajectory::getFrame(int i, const LazyEvaluator &evaluator) const
{
    if (d)
        return d->getFrame(i, evaluator);
    else
        throw SireError::invalid_index(
            QObject::tr("Invalid frame index %1. Number of frames is 0.")
                .arg(i),
            CODELOC);
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
