/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#include "trajectory.h"
#include "trajectoryaligner.h"
#include "atomidxmapping.h"

#include "SireID/index.h"

#include "SireVol/space.h"

#include "SireMol/core.h"

#include "SireMaths/align.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/generalunitproperty.h"
#include "SireBase/lazyevaluator.h"
#include "SireBase/console.h"

#include "SireBase/slice.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

////////
//////// Implementation of TrajectoryData
////////

static const RegisterMetaType<TrajectoryData> r_trajdata(MAGIC_ONLY, "SireMol::TrajectoryData");

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const TrajectoryData &trajdata)
{
    writeHeader(ds, r_trajdata, 1);
    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, TrajectoryData &trajdata)
{
    VersionID v = readHeader(ds, r_trajdata);

    if (v != 1)
        throw version_error(v, "1", r_trajdata, CODELOC);

    return ds;
}

TrajectoryData::TrajectoryData() : RefCountData()
{
}

TrajectoryData::TrajectoryData(const TrajectoryData &) : RefCountData()
{
}

TrajectoryData::~TrajectoryData()
{
}

TrajectoryData &TrajectoryData::operator=(const TrajectoryData &)
{
    return *this;
}

bool TrajectoryData::operator==(const TrajectoryData &other) const
{
    return this == &other or this->_equals(other);
}

bool TrajectoryData::operator!=(const TrajectoryData &other) const
{
    return not this->operator==(other);
}

const char *TrajectoryData::typeName()
{
    return "SireMol::TrajectoryData";
}

int TrajectoryData::nFiles() const
{
    return this->filenames().count();
}

bool TrajectoryData::isEmpty() const
{
    return this->nFrames() == 0;
}

QList<Frame> TrajectoryData::getFrames() const
{
    QList<Frame> frames;

    const int n = this->nFrames();

    frames.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        frames.append(this->getFrame(i));
    }

    return frames;
}

QList<Frame> TrajectoryData::getFrames(int start_atom, int natoms) const
{
    QList<Frame> frames;

    const int n = this->nFrames();

    frames.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        frames.append(this->getFrame(i).subset(start_atom, natoms));
    }

    return frames;
}

QList<Frame> TrajectoryData::getFrames(const LazyEvaluator &evaluator) const
{
    QList<Frame> frames;

    const int n = this->nFrames();

    frames.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        frames.append(this->getFrame(i, evaluator));
    }

    return frames;
}

QList<Frame> TrajectoryData::getFrames(int start_atom, int natoms,
                                       const LazyEvaluator &evaluator) const
{
    QList<Frame> frames;

    const int n = this->nFrames();

    frames.reserve(n);

    for (int i = 0; i < n; ++i)
    {
        frames.append(this->getFrame(i, evaluator).subset(start_atom, natoms));
    }

    return frames;
}

Frame TrajectoryData::operator[](int i) const
{
    return this->getFrame(i);
}

bool TrajectoryData::isEditable() const
{
    return false;
}

void TrajectoryData::assertIsEditable() const
{
    if (not this->isEditable())
        throw SireError::incompatible_error(QObject::tr("You cannot edit a trajectory of type %1.").arg(this->what()),
                                            CODELOC);
}

TrajectoryDataPtr TrajectoryData::makeEditable() const
{
    return TrajectoryDataPtr(new MolTrajectoryData(this->getFrames()));
}

TrajectoryDataPtr TrajectoryData::makeSubsetEditable(int start_atom, int natoms) const
{
    return TrajectoryDataPtr(new MolTrajectoryData(this->getFrames(start_atom, natoms)));
}

void TrajectoryData::setFrame(int, const Frame &)
{
    this->assertIsEditable();
}

void TrajectoryData::appendFrame(const Frame &)
{
    this->assertIsEditable();
}

void TrajectoryData::insertFrame(int, const Frame &)
{
    this->assertIsEditable();
}

void TrajectoryData::deleteFrame(int)
{
    this->assertIsEditable();
}

////////
//////// Implementation of MolTrajectoryData
////////

static const RegisterMetaType<MolTrajectoryData> r_moldata;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const MolTrajectoryData &moldata)
{
    writeHeader(ds, r_moldata, 1);

    SharedDataStream sds(ds);

    sds << moldata.frames << static_cast<const TrajectoryData &>(moldata);

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, MolTrajectoryData &moldata)
{
    VersionID v = readHeader(ds, r_moldata);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> moldata.frames >> static_cast<TrajectoryData &>(moldata);
    }
    else
        throw version_error(v, "1", r_moldata, CODELOC);

    return ds;
}

MolTrajectoryData::MolTrajectoryData() : TrajectoryData()
{
}

MolTrajectoryData::MolTrajectoryData(const Frame &frame) : TrajectoryData()
{
    if (not frame.isEmpty())
        frames.append(frame);
}

MolTrajectoryData::MolTrajectoryData(const QList<Frame> &f) : TrajectoryData()
{
    bool any_empty = false;
    int natoms = 0;

    for (const auto &frame : f)
    {
        if (frame.isEmpty())
        {
            any_empty = true;
        }

        if (natoms == 0)
        {
            natoms = frame.nAtoms();
        }
        else if (natoms != frame.nAtoms())
        {
            throw SireError::incompatible_error(
                QObject::tr("Frames needs to have the same number of atoms! %1 vs %2").arg(natoms).arg(frame.nAtoms()),
                CODELOC);
        }
    }

    frames = f;

    if (any_empty)
    {
        QMutableListIterator<Frame> it(frames);

        while (it.hasNext())
        {
            auto frame = it.next();
            if (frame.isEmpty())
            {
                it.remove();
            }
        }
    }
}

MolTrajectoryData::MolTrajectoryData(const MolTrajectoryData &other) : TrajectoryData(other), frames(other.frames)
{
}

MolTrajectoryData::~MolTrajectoryData()
{
}

MolTrajectoryData &MolTrajectoryData::operator=(const MolTrajectoryData &other)
{
    frames = other.frames;
    TrajectoryData::operator=(other);
    return *this;
}

bool MolTrajectoryData::_equals(const TrajectoryData &other) const
{
    const MolTrajectoryData *ptr = dynamic_cast<const MolTrajectoryData *>(&other);

    if (ptr)
        return frames == ptr->frames;
    else
        return false;
}

const char *MolTrajectoryData::what() const
{
    return MolTrajectoryData::typeName();
}

const char *MolTrajectoryData::typeName()
{
    return QMetaType::typeName(qMetaTypeId<MolTrajectoryData>());
}

MolTrajectoryData *MolTrajectoryData::clone() const
{
    return new MolTrajectoryData(*this);
}

int MolTrajectoryData::nFrames() const
{
    return frames.count();
}

int MolTrajectoryData::nAtoms() const
{
    if (frames.isEmpty())
        return 0;
    else
        return frames.first().nAtoms();
}

QStringList MolTrajectoryData::filenames() const
{
    return QStringList();
}

Frame MolTrajectoryData::getFrame(int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    return frames[i];
}

Frame MolTrajectoryData::getFrame(int i, const LazyEvaluator &) const
{
    return MolTrajectoryData::getFrame(i);
}

bool MolTrajectoryData::isEditable() const
{
    return true;
}

TrajectoryDataPtr MolTrajectoryData::makeEditable() const
{
    return *this;
}

TrajectoryDataPtr MolTrajectoryData::makeSubsetEditable(int start_atom, int natoms) const
{
    if (start_atom == 0 and natoms == this->nAtoms())
    {
        return *this;
    }

    start_atom = SireID::Index(start_atom).map(this->nAtoms());

    if (start_atom + natoms >= this->nAtoms())
    {
        throw SireError::incompatible_error(QObject::tr("You cannot subset a trajectory as there are not enough atoms! "
                                                        "%1 versus %2")
                                                .arg(start_atom + natoms)
                                                .arg(this->nAtoms()),
                                            CODELOC);
    }

    MolTrajectoryData ret;
    ret.frames.reserve(frames.count());

    for (const auto &frame : frames)
    {
        ret.frames.append(frame.subset(start_atom, natoms));
    }

    return ret;
}

void MolTrajectoryData::setFrame(int i, const Frame &frame)
{
    i = SireID::Index(i).map(this->nFrames());

    if (frame.nAtoms() != this->nAtoms())
    {
        throw SireError::incompatible_error(
            QObject::tr("Cannot insert the frame as the number of atoms is not the same! "
                        "%1 vs %2")
                .arg(frame.nAtoms())
                .arg(this->nAtoms()),
            CODELOC);
    }

    frames[i] = frame;
}

void MolTrajectoryData::appendFrame(const Frame &frame)
{
    if (frame.isEmpty())
        return;

    if (this->nFrames() == 0)
    {
        frames.append(frame);
    }
    else
    {
        if (frame.nAtoms() != this->nAtoms())
        {
            throw SireError::incompatible_error(
                QObject::tr("Cannot append the frame as the number of atoms is not the same! "
                            "%1 vs %2")
                    .arg(frame.nAtoms())
                    .arg(this->nAtoms()),
                CODELOC);
        }

        frames.append(frame);
    }
}

void MolTrajectoryData::insertFrame(int i, const Frame &frame)
{
    if (frame.isEmpty())
        return;

    if (i == this->nFrames())
    {
        this->appendFrame(frame);
    }
    else
    {
        i = SireID::Index(i).map(this->nFrames());

        if (frame.nAtoms() != this->nAtoms())
        {
            throw SireError::incompatible_error(
                QObject::tr("Cannot insert the frame as the number of atoms is not the same! "
                            "%1 vs %2")
                    .arg(frame.nAtoms())
                    .arg(this->nAtoms()),
                CODELOC);
        }

        frames.insert(i, frame);
    }
}

void MolTrajectoryData::deleteFrame(int i)
{
    i = SireID::Index(i).map(this->nFrames());

    frames.removeAt(i);
}

/////////////
///////////// Implementation of Trajectory
/////////////

static const RegisterMetaType<Trajectory> r_traj;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const Trajectory &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);
    sds << traj.d << traj.start_atom << traj.natoms;

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, Trajectory &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> traj.d >> traj.start_atom >> traj.natoms;
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

Trajectory::Trajectory() : ConcreteProperty<Trajectory, MoleculeProperty>(), start_atom(0), natoms(0)
{
}

Trajectory::Trajectory(const TrajectoryData &data)
    : ConcreteProperty<Trajectory, MoleculeProperty>(), start_atom(0), natoms(0)
{
    start_atom = 0;
    natoms = data.nAtoms();

    if (data.nFrames() > 0)
        d.append(TrajectoryDataPtr(data));
}

Trajectory::Trajectory(const QList<TrajectoryDataPtr> &data)
    : ConcreteProperty<Trajectory, MoleculeProperty>(), start_atom(0), natoms(0)
{
    start_atom = 0;
    natoms = 0;

    for (const auto &ptr : data)
    {
        if (ptr.constData() != 0)
        {
            if (ptr->nFrames() > 0)
            {
                if (natoms == 0)
                {
                    natoms = ptr->nAtoms();
                }
                else if (natoms != ptr->nAtoms())
                {
                    throw SireError::incompatible_error(QObject::tr("Different number of atoms in different "
                                                                    "trajectories: %1 versus %2.")
                                                            .arg(natoms)
                                                            .arg(ptr->nAtoms()),
                                                        CODELOC);
                }

                d.append(ptr);
            }
        }
    }
}

Trajectory::Trajectory(const TrajectoryData &data, int s, int n)
    : ConcreteProperty<Trajectory, MoleculeProperty>(), start_atom(s), natoms(n)
{
    if (natoms <= 0 or data.nFrames() <= 0)
        return;

    if (start_atom < 0 or (start_atom + natoms) > data.nAtoms())
        throw SireError::incompatible_error(
            QObject::tr("Cannot use start_atom %1 and natoms %2 for a trajectory with %3 atoms.")
                .arg(start_atom)
                .arg(natoms)
                .arg(data.nAtoms()),
            CODELOC);

    d.append(TrajectoryDataPtr(data));
}

Trajectory::Trajectory(const QList<TrajectoryDataPtr> &data, int s, int n)
    : ConcreteProperty<Trajectory, MoleculeProperty>(), start_atom(s), natoms(n)
{
    if (natoms <= 0 or data.count() <= 0)
        return;

    n = 0;

    for (const auto &ptr : data)
    {
        if (ptr.constData() != 0)
        {
            if (ptr->nFrames() > 0)
            {
                if (n == 0)
                {
                    n = ptr->nAtoms();
                }
                else if (n != ptr->nAtoms())
                {
                    throw SireError::incompatible_error(QObject::tr("Different number of atoms in different "
                                                                    "trajectories: %1 versus %2.")
                                                            .arg(natoms)
                                                            .arg(ptr->nAtoms()),
                                                        CODELOC);
                }

                if (start_atom + natoms > ptr->nAtoms())
                {
                    throw SireError::incompatible_error(QObject::tr("Cannot use start_atom %1 and natoms %2 "
                                                                    "for a trajectory with %3 atoms.")
                                                            .arg(start_atom)
                                                            .arg(natoms)
                                                            .arg(ptr->nAtoms()),
                                                        CODELOC);
                }

                d.append(ptr);
            }
        }
    }
}

Trajectory::Trajectory(const Trajectory &other)
    : ConcreteProperty<Trajectory, MoleculeProperty>(other), d(other.d), start_atom(other.start_atom),
      natoms(other.natoms)
{
}

Trajectory::~Trajectory()
{
}

Trajectory &Trajectory::operator=(const Trajectory &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        d = other.d;
        start_atom = other.start_atom;
        natoms = other.natoms;
    }

    return *this;
}

bool Trajectory::operator==(const Trajectory &other) const
{
    return start_atom == other.start_atom and natoms == other.natoms and d == other.d;
}

bool Trajectory::operator!=(const Trajectory &other) const
{
    return not Trajectory::operator==(other);
}

const char *Trajectory::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Trajectory>());
}

const char *Trajectory::what() const
{
    return Trajectory::typeName();
}

Trajectory *Trajectory::clone() const
{
    return new Trajectory(*this);
}

QString Trajectory::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("Trajectory::null");
    }
    else
    {
        return QObject::tr("Trajectory( num_frames=%1, num_atoms=%2 )").arg(this->nFrames()).arg(this->nAtoms());
    }
}

bool Trajectory::isEmpty() const
{
    return d.isEmpty();
}

int Trajectory::nFrames() const
{
    int nframes = 0;

    for (const auto &ptr : d)
    {
        nframes += ptr->nFrames();
    }

    return nframes;
}

int Trajectory::count() const
{
    return this->nFrames();
}

int Trajectory::size() const
{
    return this->nFrames();
}

int Trajectory::nAtoms() const
{
    return natoms;
}

// THIS WILL CHANGE 'frame' SO IT HAS THE RIGHT VALUE
// FOR THE SPECIFIED TRAJECTORY
int Trajectory::_getIndexForFrame(int &frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    for (int i = 0; i < d.count(); ++i)
    {
        int n = d[i]->nFrames();

        if (frame < n)
            return i;
        else
            frame -= n;
    }

    // we should never get here
    throw SireError::program_bug(QObject::tr("How did we get here? %1 %2").arg(frame).arg(this->nFrames()), CODELOC);

    return 0;
}

Frame Trajectory::operator[](int i) const
{
    return this->getFrame(i);
}

QList<Frame> Trajectory::operator[](const QList<qint64> &idxs) const
{
    QList<Frame> frames;
    frames.reserve(idxs.count());

    for (const auto &idx : idxs)
    {
        frames.append(this->getFrame(idx));
    }

    return frames;
}

QList<Frame> Trajectory::operator[](const SireBase::Slice &slice) const
{
    QList<Frame> ret;

    for (auto it = slice.begin(this->nFrames()); not it.atEnd(); it.next())
    {
        ret.append(this->getFrame(it.value()));
    }

    return ret;
}

Frame Trajectory::getFrame(int i) const
{
    i = Index(i).map(this->nFrames());

    int idx = _getIndexForFrame(i);

    auto frame = d[idx]->getFrame(i);

    if (frame.nAtoms() == natoms)
        return frame;
    else
        return frame.subset(start_atom, natoms);
}

Frame Trajectory::getFrame(int i, const LazyEvaluator &evaluator) const
{
    i = Index(i).map(this->nFrames());

    int idx = _getIndexForFrame(i);

    auto frame = d[idx]->getFrame(i, evaluator);

    if (frame.nAtoms() == natoms)
        return frame;
    else
        return frame.subset(start_atom, natoms);
}

Frame Trajectory::getFrame(int i, const FrameTransform &transform) const
{
    const auto nframes = this->nFrames();

    i = Index(i).map(nframes);

    auto smooth = transform.nSmooth();

    if (smooth > nframes)
        smooth = nframes;

    Frame frame = this->getFrame(i);

    if (smooth > 1)
    {
        int half = int(smooth / 2);
        int start_frame = i - half;
        int end_frame = i + half + 1;

        if (smooth % 2 == 1)
            end_frame += 1;

        if (start_frame < 0)
        {
            start_frame = 0;
        }

        if (end_frame > nframes)
        {
            end_frame = nframes;
        }

        QList<Frame> frames;

        for (int j = start_frame; j < end_frame; ++j)
        {
            if (i != j)
                frames.append(this->getFrame(j));
        }

        frame = frame.smooth(frames);
    }

    return frame.transform(transform);
}

Frame Trajectory::getFrame(int i, const FrameTransform &transform,
                           const LazyEvaluator &evaluator) const
{
    const auto nframes = this->nFrames();

    i = Index(i).map(nframes);

    auto smooth = transform.nSmooth();

    if (smooth > nframes)
        smooth = nframes;

    Frame frame = this->getFrame(i, evaluator);

    if (smooth > 1)
    {
        int half = int(smooth / 2);
        int start_frame = i - half;
        int end_frame = i + half + 1;

        if (smooth % 2 == 1)
            end_frame += 1;

        if (start_frame < 0)
        {
            start_frame = 0;
        }

        if (end_frame > nframes)
        {
            end_frame = nframes;
        }

        QList<Frame> frames;

        for (int j = start_frame; j < end_frame; ++j)
        {
            if (i != j)
                frames.append(this->getFrame(j, evaluator));
        }

        frame = frame.smooth(frames);
    }

    return frame.transform(transform);
}

TrajectoryData &Trajectory::_makeEditable(int &frame)
{
    frame = Index(frame).map(this->nFrames());

    int idx = _getIndexForFrame(frame);

    if (d[idx]->nAtoms() != natoms or (not d[idx]->isEditable()))
    {
        if (d[idx]->nAtoms() == natoms)
        {
            d[idx] = d[idx]->makeEditable();
        }
        else
        {
            d[idx] = d[idx]->makeSubsetEditable(start_atom, natoms);

            if (d[idx]->nAtoms() != natoms)
                throw SireError::program_bug(QObject::tr("The result from making the subset editable should "
                                                         "be a data with natoms = %1 (not %2)")
                                                 .arg(natoms)
                                                 .arg(d[idx]->nAtoms()),
                                             CODELOC);
        }
    }

    return *(d[idx]);
}

Frame Trajectory::_subset(const Frame &frame) const
{
    if (this->isEmpty() or frame.nAtoms() == natoms)
        return frame;
    else
        return frame.subset(start_atom, natoms);
}

void Trajectory::setFrame(int i, const Frame &frame)
{
    if (frame.isEmpty())
    {
        this->deleteFrame(i);
    }
    else
    {
        this->_makeEditable(i).setFrame(i, this->_subset(frame));
    }
}

void Trajectory::appendFrame(const Frame &frame)
{
    if (frame.isEmpty())
        return;

    auto subset = this->_subset(frame);

    if (subset.nAtoms() == 0)
    {
        throw SireError::program_bug(
            QObject::tr("This should not be empty? %1 %2 %3").arg(subset.nAtoms()).arg(start_atom).arg(natoms),
            CODELOC);
    }

    if (d.isEmpty())
    {
        d.append(TrajectoryDataPtr(new MolTrajectoryData(subset)));
        start_atom = 0;
        natoms = subset.nAtoms();
    }
    else if (d.last()->isEditable() and d.last()->nAtoms() == natoms)
    {
        d.last()->appendFrame(frame);
    }
    else
    {
        d.append(TrajectoryDataPtr(new MolTrajectoryData(subset)));
    }
}

void Trajectory::insertFrame(int i, const Frame &frame)
{
    if (frame.isEmpty())
        return;

    auto subset = this->_subset(frame);

    this->_makeEditable(i).insertFrame(i, subset);
}

void Trajectory::deleteFrame(int i)
{
    int frame = Index(i).map(this->nFrames());

    int idx = _getIndexForFrame(frame);

    this->_makeEditable(i).deleteFrame(i);

    if (d[idx]->isEmpty())
    {
        d.removeAt(idx);
    }
}

bool Trajectory::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return this->nAtoms() == molinfo.nAtoms();
}

/** Merge this property with another property */
PropertyList Trajectory::merge(const MolViewProperty &other,
                               const AtomIdxMapping &mapping,
                               const QString &ghost,
                               const SireBase::PropertyMap &map) const
{
    if (not other.isA<Trajectory>())
    {
        throw SireError::incompatible_error(QObject::tr("Cannot merge %1 with %2 as they are different types.")
                                                .arg(this->what())
                                                .arg(other.what()),
                                            CODELOC);
    }

    SireBase::Console::warning(QObject::tr("Merging %1 properties is not yet implemented. Returning two copies of the original property.")
                                   .arg(this->what()));

    SireBase::PropertyList ret;

    ret.append(*this);
    ret.append(*this);

    return ret;
}

///////
/////// Implementation of Frame
///////

static const RegisterMetaType<Frame> r_frame;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const Frame &frame)
{
    writeHeader(ds, r_frame, 2);

    SharedDataStream sds(ds);

    sds << frame.coordinates() << frame.velocities() << frame.forces() << frame.spc
        << frame.t.to(picosecond) << frame.props;

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, Frame &frame)
{
    VersionID v = readHeader(ds, r_frame);

    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> frame.coords >> frame.vels >> frame.frcs >> frame.spc;

        double time;
        sds >> time;

        frame.t = time * picosecond;
        frame.start_atom = 0;
        frame.num_atoms = std::max(frame.coords.count(),
                                   std::max(frame.vels.count(),
                                            frame.frcs.count()));

        sds >> frame.props;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> frame.coords >> frame.vels >> frame.frcs >> frame.spc;

        double time;
        sds >> time;

        frame.t = time * picosecond;
        frame.start_atom = 0;
        frame.num_atoms = std::max(frame.coords.count(),
                                   std::max(frame.vels.count(),
                                            frame.frcs.count()));

        frame.props = Properties();
    }
    else
        throw version_error(v, "1", r_frame, CODELOC);

    return ds;
}

Frame::Frame()
    : ConcreteProperty<Frame, MoleculeProperty>(),
      spc(Space::null()), t(0), start_atom(-1), num_atoms(0)
{
}

Frame::Frame(const Molecule &mol, const PropertyMap &map)
    : ConcreteProperty<Frame, MoleculeProperty>(),
      spc(Space::null()), t(0), start_atom(-1), num_atoms(0)
{
    this->operator=(Frame(mol.data(), map));
}

Frame::Frame(const MoleculeData &mol, const PropertyMap &map)
    : ConcreteProperty<Frame, MoleculeProperty>(),
      spc(Space::null()), t(0), start_atom(-1), num_atoms(0)
{
    int natoms = 0;
    int broken = false;

    if (mol.hasProperty(map["coordinates"]))
    {
        try
        {
            coords = mol.property(map["coordinates"]).asA<AtomCoords>().toVector();
            natoms = coords.count();
        }
        catch (const SireError::exception &e)
        {
            qDebug() << "WARNING: Problem saving coordinates!";
            qDebug() << e.what() << e.error();
        }
    }

    if (mol.hasProperty(map["velocities"]))
    {
        try
        {
            vels = mol.property(map["velocities"]).asA<AtomVelocities>().toVector();

            if (natoms == 0)
                natoms = vels.count();
            else if (natoms != vels.count())
            {
                broken = true;
            }
        }
        catch (const SireError::exception &e)
        {
            qDebug() << "WARNING: Problem saving velocities!";
            qDebug() << e.what() << e.error();
        }
    }

    if (mol.hasProperty(map["forces"]))
    {
        try
        {
            frcs = mol.property(map["forces"]).asA<AtomForces>().toVector();

            if (natoms == 0)
                natoms = frcs.count();
            else if (natoms != frcs.count())
            {
                broken = true;
            }
        }
        catch (const SireError::exception &e)
        {
            qDebug() << "WARNING: Problem saving forces!";
            qDebug() << e.what() << e.error();
        }
    }

    if (broken)
    {
        throw SireError::incompatible_error(
            QObject::tr("Incompatibility between the coordinates, velocities and forces? "
                        "%1 | %2 | %3")
                .arg(coords.count())
                .arg(vels.count())
                .arg(frcs.count()),
            CODELOC);
    }

    if (mol.hasProperty(map["time"]))
    {
        try
        {
            t = mol.property(map["time"]).asA<GeneralUnitProperty>();
        }
        catch (...)
        {
        }
    }

    if (mol.hasProperty(map["space"]))
    {
        try
        {
            spc = mol.property(map["space"]).asA<Space>();
        }
        catch (...)
        {
        }
    }

    num_atoms = std::max(coords.count(), std::max(vels.count(),
                                                  frcs.count()));

    this->assertSane();
}

Frame::Frame(const QVector<Vector> &coordinates, const Space &space, Time time)
    : ConcreteProperty<Frame, MoleculeProperty>(),
      coords(coordinates), spc(space), t(time),
      start_atom(-1), num_atoms(coordinates.count())
{
    this->assertSane();
}

Frame::Frame(const QVector<Vector> &coordinates, const QVector<Velocity3D> &velocities, const Space &space, Time time)
    : ConcreteProperty<Frame, MoleculeProperty>(),
      coords(coordinates), vels(velocities), spc(space), t(time),
      start_atom(-1)
{
    num_atoms = std::max(coords.count(), vels.count());

    this->assertSane();
}

Frame::Frame(const QVector<Vector> &coordinates, const QVector<Velocity3D> &velocities, const QVector<Force3D> &forces,
             const Space &space, Time time)
    : ConcreteProperty<Frame, MoleculeProperty>(), coords(coordinates), vels(velocities), frcs(forces), spc(space),
      t(time), start_atom(-1)
{
    num_atoms = std::max(coords.count(), std::max(vels.count(),
                                                  frcs.count()));

    this->assertSane();
}

Frame::Frame(const QVector<Vector> &coordinates, const QVector<Velocity3D> &velocities, const QVector<Force3D> &forces,
             const Space &space, Time time, const Properties &properties)
    : ConcreteProperty<Frame, MoleculeProperty>(), coords(coordinates), vels(velocities),
      frcs(forces), spc(space), t(time), props(properties),
      start_atom(-1)
{
    num_atoms = std::max(coords.count(), std::max(vels.count(),
                                                  frcs.count()));

    this->assertSane();
}

Frame::Frame(const Frame &other)
    : ConcreteProperty<Frame, MoleculeProperty>(), coords(other.coords), vels(other.vels), frcs(other.frcs),
      spc(other.spc), t(other.t), props(other.props),
      start_atom(other.start_atom), num_atoms(other.num_atoms)
{
}

Frame::~Frame()
{
}

Frame &Frame::operator=(const Frame &other)
{
    if (this != &other)
    {
        coords = other.coords;
        vels = other.vels;
        frcs = other.frcs;
        spc = other.spc;
        t = other.t;
        start_atom = other.start_atom;
        num_atoms = other.num_atoms;
        props = other.props;
    }

    return *this;
}

bool Frame::operator==(const Frame &other) const
{
    return (this == &other) or
           (start_atom == other.start_atom and
            num_atoms == other.num_atoms and
            coords == other.coords and
            vels == other.vels and
            frcs == other.frcs and
            spc == other.spc and
            t == other.t and props == other.props);
}

bool Frame::operator!=(const Frame &other) const
{
    return not operator==(other);
}

const char *Frame::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Frame>());
}

const char *Frame::what() const
{
    return Frame::typeName();
}

Frame *Frame::clone() const
{
    return new Frame(*this);
}

void Frame::assertSane() const
{
    int natoms = coords.count();
    bool broken = false;

    if (natoms == 0)
    {
        natoms = vels.count();
    }
    else if (vels.count() != 0 and natoms != vels.count())
    {
        broken = true;
    }

    if (natoms == 0)
    {
        natoms = frcs.count();
    }
    else if (frcs.count() != 0 and natoms != frcs.count())
    {
        broken = true;
    }

    if (natoms != num_atoms)
        broken = true;

    if (broken)
    {
        throw SireError::program_bug(QObject::tr("Somehow the Frame has got into an invalid state! "
                                                 "%1 | %2 | %3 | %4")
                                         .arg(coords.count())
                                         .arg(vels.count())
                                         .arg(frcs.count())
                                         .arg(num_atoms),
                                     CODELOC);
    }

    if (t.value() < 0)
    {
        throw SireError::program_bug(
            QObject::tr("The time of this frame is negative! This is not supported. %1").arg(t.toString()), CODELOC);
    }
}

QString Frame::toString() const
{
    QStringList c;

    if (this->hasCoordinates())
        c.append("C");

    if (this->hasVelocities())
        c.append("V");

    if (this->hasForces())
        c.append("F");

    if (this->hasProperties())
        c.append("P");

    if (c.isEmpty())
        return QObject::tr("Frame::null");

    return QObject::tr("Frame( nAtoms=%1, time=%2, {%3} )")
        .arg(this->nAtoms())
        .arg(this->time().toString())
        .arg(c.join(","));
}

bool Frame::isEmpty() const
{
    return num_atoms == 0 and coords.isEmpty() and
           vels.isEmpty() and frcs.isEmpty() and props.isEmpty();
}

template <class T>
QVector<T> _reorder(const QVector<T> &orig, const QVector<qint64> &order)
{
    if (orig.isEmpty() or order.isEmpty())
        return orig;

    QVector<T> ret = orig;

    for (int i = 0; i < order.count(); ++i)
    {
        int j = order[i];

        if (j >= 0 and j < orig.count())
            ret[i] = orig[j];
    }

    return ret;
}

/** Return a copy of this frame which has been reordered according
 *  to 'order' (i.e. atom i is moved to order[i]). This does nothing
 *  if the order is empty. It silently ignores invalid orders, and
 *  will leave atoms that aren't referenced in their original
 *  positions
 */
Frame Frame::reorder(const QVector<qint64> &order) const
{
    Frame ret(*this);

    if (order.isEmpty())
        return ret;

    ret.coords = _reorder(ret.coords, order);
    ret.vels = _reorder(ret.vels, order);
    ret.frcs = _reorder(ret.frcs, order);

    return ret;
}

Frame Frame::extract() const
{
    if (start_atom == -1)
        return *this;
    else
    {
        Frame ret(*this);
        ret.coords = this->coordinates();
        ret.vels = this->velocities();
        ret.frcs = this->forces();
        ret.start_atom = -1;
        ret.num_atoms = std::max(ret.coords.count(),
                                 std::max(ret.vels.count(),
                                          ret.frcs.count()));

        return ret;
    }
}

Frame Frame::transform(const FrameTransform &tform) const
{
    Frame ret = this->extract();
    ret.coords = tform.apply(ret.coords, ret.spc.read());
    ret.spc = tform.apply(ret.spc.read());
    return ret;
}

Frame Frame::reverse(const FrameTransform &tform) const
{
    Frame ret = this->extract();
    ret.coords = tform.reverse(ret.coords);
    ret.spc = tform.reverse(ret.spc.read());
    return ret;
}

/** Return the frame which has this frame smoothed with the coordinates
 *  from all of the other frames. The result will be returned
 */
Frame Frame::smooth(const QList<Frame> &frames) const
{
    if (frames.isEmpty() or this->coords.isEmpty())
    {
        return *this;
    }

    Frame ret(*this);
    ret.vels = this->velocities();
    ret.frcs = this->forces();
    ret.start_atom = -1;
    ret.num_atoms = this->num_atoms;

    const double weight = 1.0 / (1 + frames.count());

    QVector<Vector> smoothed(this->num_atoms);
    auto smoothed_data = smoothed.data();

    Vector mincoords = this->coords[0];
    Vector maxcoords = this->coords[0];

    const Vector *coords_data = this->coordinatesData();

    for (int i = 0; i < this->num_atoms; ++i)
    {
        const auto &c = coords_data[i];

        mincoords.setMin(c);
        maxcoords.setMax(c);

        smoothed_data[i] = weight * c;
    }

    Vector center = mincoords + (0.5 * (maxcoords - mincoords));

    for (int i = 0; i < frames.count(); ++i)
    {
        const auto frame_coords_data = frames.at(i).coordinatesData();

        if (frames.at(i).nAtoms() != num_atoms)
        {
            throw SireError::incompatible_error(
                QObject::tr(
                    "Cannot smooth a trajectory if the number of "
                    "atoms / coordinates per frame changes."),
                CODELOC);
        }

        Vector delta = frame_coords_data[0] - coords_data[0];

        if (delta.length2() > 30)
        {
            // the first atom has jumped by 5-6 A. This suggests that this
            // frame may have been wrapped. We will now check the delta
            // between the centers of the two frames...
            mincoords = frame_coords_data[0];
            maxcoords = mincoords;

            for (int j = 1; j < num_atoms; ++j)
            {
                mincoords.setMin(frame_coords_data[j]);
                maxcoords.setMax(frame_coords_data[j]);
            }

            Vector frame_center = mincoords + (0.5 * (maxcoords - mincoords));

            delta = center - frame_center;

            if (delta.length2() > 30)
            {
                // Yes - the center has moved too - suggests a frame offset
                // Translate to remove the center offset
                for (int j = 0; j < num_atoms; ++j)
                {
                    smoothed_data[j] += weight * (frame_coords_data[j] + delta);
                }
            }
            else
            {
                // No - not enough of a change - just use the
                // original coordinates
                for (int j = 0; j < num_atoms; ++j)
                {
                    smoothed_data[j] += weight * frame_coords_data[j];
                }
            }
        }
        else
        {
            for (int j = 0; j < num_atoms; ++j)
            {
                smoothed_data[j] += weight * frame_coords_data[j];
            }
        }
    }

    ret.coords = smoothed;
    return ret;
}

bool Frame::hasCoordinates() const
{
    return not coords.isEmpty();
}

bool Frame::hasVelocities() const
{
    return not vels.isEmpty();
}

bool Frame::hasForces() const
{
    return not frcs.isEmpty();
}

bool Frame::hasProperties() const
{
    return not props.isEmpty();
}

template <class T>
QVector<T> _mid(const QVector<T> &v, int start, int count)
{
    if (v.count() == 0)
        return v;

    if (start < 0 or count < 0 or start + count > v.count())
        throw SireError::incompatible_error(QObject::tr("Cannot subset a vector of length %1 using start %2, count %3.")
                                                .arg(v.count())
                                                .arg(start)
                                                .arg(count),
                                            CODELOC);

    return v.mid(start, count);
}

const Vector *Frame::coordinatesData() const
{
    if (coords.isEmpty())
        return 0;
    else if (start_atom == -1)
        return coords.constData();
    else
        return coords.constData() + start_atom;
}

const Velocity3D *Frame::velocitiesData() const
{
    if (vels.isEmpty())
        return 0;
    else if (start_atom == -1)
        return vels.constData();
    else
        return vels.constData() + start_atom;
}

const Force3D *Frame::forcesData() const
{
    if (frcs.isEmpty())
        return 0;
    else if (start_atom == -1)
        return frcs.constData();
    else
        return frcs.constData() + start_atom;
}

QVector<Vector> Frame::coordinates() const
{
    if (start_atom == -1)
        return coords;
    else
        return _mid(coords, start_atom, num_atoms);
}

QVector<Velocity3D> Frame::velocities() const
{
    if (start_atom == -1)
        return vels;
    else
        return _mid(vels, start_atom, num_atoms);
}

QVector<Force3D> Frame::forces() const
{
    if (start_atom == -1)
        return frcs;
    else
        return _mid(frcs, start_atom, num_atoms);
}

SireUnits::Dimension::Time Frame::time() const
{
    return t;
}

const Space &Frame::space() const
{
    return *spc;
}

const Properties &Frame::properties() const
{
    return props;
}

bool Frame::hasProperty(const PropertyName &key) const
{
    return props.hasProperty(key);
}

const Property &Frame::property(const PropertyName &key) const
{
    return props.property(key);
}

const Property &Frame::property(const SireBase::PropertyName &key,
                                const SireBase::Property &default_value) const
{
    return props.property(key, default_value);
}

int Frame::numBytes() const
{
    // 3 doubles times the number of coordinates, velocities and forces
    return 24 * (coords.count() + vels.count() + frcs.count());
}

int Frame::nAtoms() const
{
    return num_atoms;
}

Frame Frame::subset(int start, int natoms) const
{
    if (start == 0 and natoms == this->nAtoms())
        return *this;

    Frame ret(*this);

    if (start < 0 or natoms < 0 or start + natoms > num_atoms)
        throw SireError::incompatible_error(QObject::tr("Cannot subset a Frame of length %1 using start %2, count %3.")
                                                .arg(num_atoms)
                                                .arg(start)
                                                .arg(natoms),
                                            CODELOC);

    if (start_atom == -1)
    {
        ret.start_atom = start;
        ret.num_atoms = natoms;
    }
    else
    {
        ret.start_atom = start_atom + start;
        ret.num_atoms = natoms;
    }

    return ret;
}

bool Frame::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    return this->nAtoms() == molinfo.nAtoms();
}

/** Join the vector of passed frames into a single frame. The
 *  frames are joined in order, e.g. from the first atom in
 *  the first frame to the last atom in the last frame
 */
Frame Frame::join(const QVector<Frame> &frames,
                  bool use_parallel)
{
    if (frames.isEmpty())
    {
        return Frame();
    }
    else if (frames.count() == 1)
    {
        return frames.at(0);
    }

    const int nframes = frames.count();
    const Frame *frames_data = frames.constData();

    QVector<int> start_idxs;
    start_idxs.reserve(nframes);

    bool have_coords = false;
    bool have_vels = false;
    bool have_frcs = false;

    int natoms = 0;

    for (const auto &frame : frames)
    {
        start_idxs.append(natoms);

        natoms += frame.nAtoms();

        if (frame.hasCoordinates())
            have_coords = true;

        if (frame.hasVelocities())
            have_vels = true;

        if (frame.hasForces())
            have_frcs = true;
    }

    const int *start_idxs_data = start_idxs.constData();

    // take all of the global data from the values
    // of the first frame
    Frame ret(frames.at(0));
    ret = ret.extract();

    Vector *coords_data = 0;
    Velocity3D *vels_data = 0;
    Force3D *frcs_data = 0;

    if (have_coords)
    {
        ret.coords = QVector<Vector>(natoms, Vector(0));
        coords_data = ret.coords.data();
    }
    else
    {
        ret.coords = QVector<Vector>();
    }

    if (have_vels)
    {
        ret.vels = QVector<Velocity3D>(natoms, Velocity3D(0));
        vels_data = ret.vels.data();
    }
    else
    {
        ret.vels = QVector<Velocity3D>();
    }

    if (have_frcs)
    {
        ret.frcs = QVector<Force3D>(natoms, Force3D(0));
        frcs_data = ret.frcs.data();
    }
    else
    {
        ret.frcs = QVector<Force3D>();
    }

    if (use_parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nframes), [&](tbb::blocked_range<int> r)
                          {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                const Frame frame = frames_data[i].extract();
                const int start_idx = start_idxs_data[i];

                if (have_coords and frame.hasCoordinates())
                    std::memcpy(coords_data + start_idx,
                                frame.coordinates().constData(),
                                frame.nAtoms() * sizeof(Vector));

                if (have_vels and frame.hasVelocities())
                    std::memcpy(vels_data + start_idx,
                                frame.velocities().constData(),
                                frame.nAtoms() * sizeof(Velocity3D));

                if (have_frcs and frame.hasForces())
                    std::memcpy(frcs_data + start_idx,
                                frame.forces().constData(),
                                frame.nAtoms() * sizeof(Force3D));
            } });
    }
    else
    {
        for (int i = 0; i < nframes; ++i)
        {
            const Frame frame = frames_data[i].extract();
            const int start_idx = start_idxs_data[i];

            if (have_coords and frame.hasCoordinates())
                std::memcpy(coords_data + start_idx,
                            frame.coordinates().constData(),
                            frame.nAtoms() * sizeof(Vector));

            if (have_vels and frame.hasVelocities())
                std::memcpy(vels_data + start_idx,
                            frame.velocities().constData(),
                            frame.nAtoms() * sizeof(Velocity3D));

            if (have_frcs and frame.hasForces())
                std::memcpy(frcs_data + start_idx,
                            frame.forces().constData(),
                            frame.nAtoms() * sizeof(Force3D));
        }
    }

    return ret;
}

/** Merge this property with another property */
PropertyList Frame::merge(const MolViewProperty &other,
                          const AtomIdxMapping &mapping,
                          const QString &ghost,
                          const SireBase::PropertyMap &map) const
{
    if (not other.isA<Frame>())
    {
        throw SireError::incompatible_error(QObject::tr("Cannot merge %1 with %2 as they are different types.")
                                                .arg(this->what())
                                                .arg(other.what()),
                                            CODELOC);
    }

    SireBase::Console::warning(QObject::tr("Merging %1 properties is not yet implemented. Returning two copies of the original property.")
                                   .arg(this->what()));

    SireBase::PropertyList ret;

    ret.append(*this);
    ret.append(*this);

    return ret;
}
