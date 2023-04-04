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

#include "trajectoryaligner.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<TrajectoryAligner> r_ta;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds,
                                       const TrajectoryAligner &ta)
{
    writeHeader(ds, r_ta, 1);

    SharedDataStream sds(ds);
    sds << ta.atms << ta.refcoords;

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, TrajectoryAligner &ta)
{
    VersionID v = readHeader(ds, r_ta);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> ta.atms >> ta.refcoords;
    }
    else
        throw version_error(v, "1", r_ta, CODELOC);

    return ds;
}

TrajectoryAligner::TrajectoryAligner()
{
}

TrajectoryAligner::TrajectoryAligner(const SelectorM<Atom> &atoms,
                                     const PropertyMap &m)
    : ConcreteProperty<TrajectoryAligner, Property>(),
      atms(atoms), map(m)
{
    if (atms.count() > 0)
    {
        const auto c = atms.property<Vector>(map["coordinates"]);
        refcoords = QVector<Vector>(c.begin(), c.end());
    }
}

TrajectoryAligner::TrajectoryAligner(const TrajectoryAligner &other)
    : ConcreteProperty<TrajectoryAligner, Property>(other),
      atms(other.atms), refcoords(other.refcoords), map(other.map)
{
}

TrajectoryAligner::~TrajectoryAligner()
{
}

TrajectoryAligner &TrajectoryAligner::operator=(const TrajectoryAligner &other)
{
    if (this != &other)
    {
        atms = other.atms;
        refcoords = other.refcoords;
        map = other.map;
    }

    return *this;
}

bool TrajectoryAligner::operator==(const TrajectoryAligner &other) const
{
    return atms == other.atms and refcoords == other.refcoords and
           map == other.map;
}

bool TrajectoryAligner::operator!=(const TrajectoryAligner &other) const
{
    return not this->operator==(other);
}

TrajectoryAligner *TrajectoryAligner::clone() const
{
    return new TrajectoryAligner(*this);
}

const char *TrajectoryAligner::what() const
{
    return TrajectoryAligner::typeName();
}

const char *TrajectoryAligner::typeName()
{
    return QMetaType::typeName(qMetaTypeId<TrajectoryAligner>());
}

QString TrajectoryAligner::toString() const
{
    return QObject::tr("TrajectoryAligner()");
}

const SelectorM<Atom> &TrajectoryAligner::atoms() const
{
    return atms;
}

Transform TrajectoryAligner::operator[](int i) const
{
    // get the coordinates for this frame
    const auto nframes = this->count();

    if (atms.count() == 0 or nframes <= 1)
        // no alignment possible
        return Transform();

    // make sure that we cap if we run out of frames
    if (i >= nframes)
        i = nframes - 1;
    else if (i < -nframes)
        i = -nframes;

    i = Index(i).map(nframes);

    // the below part could be cached :-)
    auto c = SelectorM<Atom>(atms);
    c.loadFrame(i, map);

    const auto coords = c.property<Vector>(map["coordinates"]);

    return SireMaths::getAlignment(refcoords,
                                   QVector<Vector>(coords.begin(), coords.end()),
                                   true);
}

QList<Transform> TrajectoryAligner::operator[](const QList<qint64> &idxs) const
{
    QList<Transform> ret;

    for (const auto &idx : idxs)
    {
        ret.append(this->operator[](idx));
    }

    return ret;
}

QList<Transform> TrajectoryAligner::operator[](const SireBase::Slice &slice) const
{
    QList<Transform> ret;

    for (auto it = slice.begin(this->nFrames()); not it.atEnd(); it.next())
    {
        ret.append(this->operator[](it.value()));
    }

    return ret;
}

int TrajectoryAligner::count() const
{
    return atms.nFrames(map);
}

int TrajectoryAligner::size() const
{
    return this->count();
}

int TrajectoryAligner::nFrames() const
{
    return this->count();
}

bool TrajectoryAligner::isCompatibleWith(const MoleculeInfoData &) const
{
    return true;
}
