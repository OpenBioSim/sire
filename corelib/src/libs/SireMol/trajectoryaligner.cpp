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
#include "atommatch.h"

#include "SireMaths/vector.h"

#include "SireBase/numberproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of FrameTransform
////////

static const RegisterMetaType<FrameTransform> r_ft;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds,
                                       const FrameTransform &ft)
{
    writeHeader(ds, r_ft, 1);

    SharedDataStream sds(ds);

    sds << ft.tform << ft.cent << ft.old_space << ft.new_space
        << ft.smooth << ft.autowrap;

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds,
                                       FrameTransform &ft)
{
    VersionID v = readHeader(ds, r_ft);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> ft.tform >> ft.cent >> ft.old_space >> ft.new_space >> ft.smooth >> ft.autowrap;
    }
    else
        throw version_error(v, "1", r_ft, CODELOC);

    return ds;
}

FrameTransform::FrameTransform()
    : ConcreteProperty<FrameTransform, Property>(),
      cent(0), smooth(1), autowrap(false)
{
}

FrameTransform::FrameTransform(const SireBase::PropertyMap &map)
    : ConcreteProperty<FrameTransform, Property>(),
      cent(0), smooth(1), autowrap(false)
{
    if (map.specified("wrap"))
    {
        autowrap = map["wrap"].value().asABoolean();
    }

    if (map.specified("smooth"))
    {
        smooth = map["smooth"].value().asAnInteger();
    }
}

FrameTransform::FrameTransform(const SireMaths::Transform &transform,
                               const SireMaths::Vector &center,
                               const SireVol::Space &space,
                               int sm,
                               bool wrap)
    : ConcreteProperty<FrameTransform, Property>(),
      tform(transform), cent(center), old_space(space),
      smooth(sm), autowrap(wrap)
{
    if (tform.isNull())
        new_space = old_space;
    else
        new_space = old_space.read().transform(tform);

    if (smooth < 1)
        smooth = 1;
}

FrameTransform::FrameTransform(const FrameTransform &other)
    : ConcreteProperty<FrameTransform, Property>(),
      tform(other.tform), cent(other.cent),
      old_space(other.old_space), new_space(other.new_space),
      smooth(other.smooth), autowrap(other.autowrap)
{
}

FrameTransform::~FrameTransform()
{
}

FrameTransform &FrameTransform::operator=(const FrameTransform &other)
{
    if (this != &other)
    {
        tform = other.tform;
        cent = other.cent;
        old_space = other.old_space;
        new_space = other.new_space;
        smooth = other.smooth;
        autowrap = other.autowrap;
    }

    return *this;
}

bool FrameTransform::operator==(const FrameTransform &other) const
{
    return this == &other or
           (tform == other.tform and cent == other.cent and
            old_space == other.old_space and new_space == other.new_space and
            smooth == other.smooth and autowrap == other.autowrap);
}

bool FrameTransform::operator!=(const FrameTransform &other) const
{
    return not this->operator==(other);
}

FrameTransform *FrameTransform::clone() const
{
    return new FrameTransform(*this);
}

const char *FrameTransform::what() const
{
    return FrameTransform::typeName();
}

const char *FrameTransform::typeName()
{
    return QMetaType::typeName(qMetaTypeId<FrameTransform>());
}

QString FrameTransform::toString() const
{
    return QObject::tr("FrameTransform( %1, nSmooth=%2, wrap=%3 )")
        .arg(this->transform().toString())
        .arg(this->nSmooth())
        .arg(this->wrap());
}

SireVol::SpacePtr FrameTransform::apply(const SireVol::Space &space) const
{
    if (this->transform().isNull())
    {
        return space;
    }
    else if (this->old_space.read().equals(space))
    {
        return this->new_space.read();
    }
    else
    {
        return space.transform(this->transform(), true);
    }
}

SireMaths::Vector FrameTransform::apply(const SireMaths::Vector &coords) const
{
    Vector ret = coords;

    if (this->wrap())
    {
        ret = this->space().getMinimumImage(ret, this->center());
    }

    return this->transform().apply(ret);
}

QVector<SireMaths::Vector> FrameTransform::apply(const QVector<Vector> &coords) const
{
    QVector<Vector> ret = coords;

    if (this->wrap())
    {
        // wrap these coordinates into the old space before we
        // perform the transfrom
        ret = this->space().getMinimumImage(ret, this->center());
    }

    // now perform the transformation
    this->transform().apply(ret.data(), ret.count());

    return ret;
}

SireMaths::Vector FrameTransform::apply(const SireMaths::Vector &coords,
                                        const SireVol::Space &spc) const
{
    Vector ret = coords;

    if (this->wrap())
    {
        ret = spc.getMinimumImage(ret, this->center());
    }

    return this->transform().apply(ret);
}

QVector<SireMaths::Vector> FrameTransform::apply(const QVector<Vector> &coords,
                                                 const SireVol::Space &spc) const
{
    QVector<Vector> ret = coords;

    if (this->wrap())
    {
        // wrap these coordinates into the old space before we
        // perform the transfrom
        ret = spc.getMinimumImage(ret, this->center());
    }

    // now perform the transformation
    this->transform().apply(ret.data(), ret.count());

    return ret;
}

Frame FrameTransform::apply(const Frame &frame) const
{
    return frame.transform(*this);
}

SireVol::SpacePtr FrameTransform::reverse(const SireVol::Space &space) const
{
    if (this->transform().isNull())
    {
        return space;
    }
    else if (this->new_space.read().equals(space))
    {
        return this->old_space.read();
    }
    else
    {
        return space.transform(this->transform(), false);
    }
}

SireMaths::Vector FrameTransform::reverse(const SireMaths::Vector &coords) const
{
    Vector ret = coords;

    if (this->wrap())
    {
        ret = this->new_space.read().getMinimumImage(ret, this->center());
    }

    return this->transform().reverse(ret);
}

QVector<SireMaths::Vector> FrameTransform::reverse(const QVector<SireMaths::Vector> &coords) const
{
    QVector<Vector> ret = coords;

    if (this->wrap())
    {
        ret = this->new_space.read().getMinimumImage(ret, this->center());
    }

    this->transform().reverse(ret.data(), ret.count());

    return ret;
}

Frame FrameTransform::reverse(const Frame &frame) const
{
    return frame.reverse(*this);
}

const SireMaths::Transform &FrameTransform::transform() const
{
    return this->tform;
}

const SireMaths::Vector &FrameTransform::center() const
{
    return this->cent;
}

const SireVol::Space &FrameTransform::space() const
{
    return this->old_space.read();
}

int FrameTransform::nSmooth() const
{
    return this->smooth;
}

bool FrameTransform::wrap() const
{
    return this->autowrap;
}

bool FrameTransform::isCompatibleWith(const MoleculeInfoData &) const
{
    return true;
}

/////////
///////// Implementation of TrajectoryAligner
/////////

static const RegisterMetaType<TrajectoryAligner> r_ta;

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds,
                                       const TrajectoryAligner &ta)
{
    writeHeader(ds, r_ta, 1);

    SharedDataStream sds(ds);
    sds << ta.atms << ta.refcoords << ta.cent << ta.nsmooth << ta.map;

    return ds;
}

SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, TrajectoryAligner &ta)
{
    VersionID v = readHeader(ds, r_ta);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> ta.atms >> ta.refcoords >> ta.cent >> ta.nsmooth >> ta.map;
    }
    else
        throw version_error(v, "1", r_ta, CODELOC);

    return ds;
}

TrajectoryAligner::TrajectoryAligner() : nsmooth(0), autowrap(false)
{
}

void TrajectoryAligner::_populate(const PropertyMap &m)
{
    this->map = m;

    if (this->map.specified("smooth"))
    {
        this->nsmooth = this->map["smooth"].value().asAnInteger();

        if (this->nsmooth < 1)
            this->nsmooth = 1;
    }

    if (this->map.specified("wrap"))
    {
        this->autowrap = this->map["wrap"].value().asABoolean();
    }

    // now modify the map that is used to get the reference frame
    // to save the number of frames to smooth over, and to turn off
    // wrapping (as we don't want to wrap when we get the reference
    // coordinates)
    this->map.set("smooth", NumberProperty(qint64(this->nsmooth)));
    this->map.set("wrap", BooleanProperty(false));
}

TrajectoryAligner::TrajectoryAligner(const SireMaths::Vector &center,
                                     const PropertyMap &m)
    : ConcreteProperty<TrajectoryAligner, Property>(),
      cent(0), nsmooth(1), autowrap(true)
{
    cent = center;
    this->_populate(m);
}

TrajectoryAligner::TrajectoryAligner(const SelectorM<Atom> &atoms,
                                     const PropertyMap &m)
    : ConcreteProperty<TrajectoryAligner, Property>(),
      atms(atoms), cent(0), nsmooth(1), map(m), autowrap(true)
{
    if (atms.count() > 0)
    {
        const auto c = atms.property<Vector>(map["coordinates"]);
        refcoords = QVector<Vector>(c.begin(), c.end());

        cent = refcoords[0];

        if (refcoords.count() > 1)
        {
            Vector mincoords = cent;
            Vector maxcoords = cent;

            for (int i = 1; i < refcoords.count(); ++i)
            {
                mincoords.setMin(refcoords[i]);
                maxcoords.setMax(refcoords[i]);
            }

            cent = mincoords + 0.5 * (maxcoords - mincoords);
        }
    }

    this->_populate(m);
}

TrajectoryAligner::TrajectoryAligner(const SelectorM<Atom> &atoms,
                                     const QVector<Vector> &points,
                                     const PropertyMap &m)
    : ConcreteProperty<TrajectoryAligner, Property>(),
      atms(atoms), cent(0), nsmooth(1), map(m), autowrap(true)
{
    if (atms.count() > 0)
    {
        if (points.count() != atms.count())
            throw SireError::incompatible_error(QObject::tr(
                                                    "Not enought reference points (%1) to align this number "
                                                    "of atoms (%2)")
                                                    .arg(points.count())
                                                    .arg(atoms.count()),
                                                CODELOC);

        refcoords = points;

        cent = refcoords[0];

        if (refcoords.count() > 1)
        {
            Vector mincoords = cent;
            Vector maxcoords = cent;

            for (int i = 1; i < refcoords.count(); ++i)
            {
                mincoords.setMin(refcoords[i]);
                maxcoords.setMax(refcoords[i]);
            }

            cent = mincoords + 0.5 * (maxcoords - mincoords);
        }
    }

    this->_populate(m);
}

TrajectoryAligner::TrajectoryAligner(const TrajectoryAligner &other)
    : ConcreteProperty<TrajectoryAligner, Property>(other),
      atms(other.atms), refcoords(other.refcoords),
      cent(other.cent), nsmooth(other.nsmooth), map(other.map),
      autowrap(other.autowrap)
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
        cent = other.cent;
        nsmooth = other.nsmooth;
        map = other.map;
        autowrap = other.autowrap;
    }

    return *this;
}

bool TrajectoryAligner::operator==(const TrajectoryAligner &other) const
{
    return atms == other.atms and refcoords == other.refcoords and
           cent == other.cent and nsmooth == other.nsmooth and
           map == other.map and autowrap == other.autowrap;
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
    if (this->nSmooth() > 1)
    {
        if (this->wrap())
            return QObject::tr("TrajectoryAligner(center=%1, nSmooth=%2, wrapped)")
                .arg(this->center().toString())
                .arg(this->nSmooth());
        else
            return QObject::tr("TrajectoryAligner(center=%1, nSmooth=%2)")
                .arg(this->center().toString())
                .arg(this->nSmooth());
    }
    else
    {
        if (this->wrap())
            return QObject::tr("TrajectoryAligner(center=%1, wrapped)")
                .arg(this->center().toString());
        else
            return QObject::tr("TrajectoryAligner(center=%1)")
                .arg(this->center().toString());
    }
}

int TrajectoryAligner::nSmooth() const
{
    return this->nsmooth;
}

Vector TrajectoryAligner::center() const
{
    return this->cent;
}

bool TrajectoryAligner::wrap() const
{
    return this->autowrap;
}

const SelectorM<Atom> &TrajectoryAligner::atoms() const
{
    return atms;
}

FrameTransform TrajectoryAligner::operator[](int i) const
{
    // get the coordinates for this frame
    const auto nframes = this->count();

    const auto space_property = map["space"];

    SireVol::SpacePtr space = SireVol::Cartesian();

    if (atms.count() == 0 or nframes <= 1)
    {
        // use the first space that we find
        if (this->wrap() and (atms.count() > 0))
        {
            const auto mols = atms.molecules();

            for (int i = 0; i < mols.count(); ++i)
            {
                const auto &mol = mols[i];

                if (mol.data().hasProperty(space_property))
                {
                    space = mol.data().property(space_property).asA<Space>();
                    break;
                }
            }
        }

        // no alignment possible - just potential for mapping into
        // the space
        return FrameTransform(SireMaths::Transform(),
                              this->center(),
                              space.read(),
                              this->nSmooth(),
                              this->wrap());
    }

    // make sure that we cap if we run out of frames
    if (i >= nframes)
        i = nframes - 1;
    else if (i < -nframes)
        i = -nframes;

    i = Index(i).map(nframes);

    // the below part could be cached :-)
    auto c = SelectorM<Atom>(atms);

    // this will load the frame - smoothing the coordinates
    // if needed
    c.loadFrame(i, this->map);

    // use the first space that we find
    if (this->wrap() and (c.count() > 0))
    {
        const auto mols = c.molecules();

        for (int i = 0; i < mols.count(); ++i)
        {
            const auto &mol = mols[i];

            if (mol.data().hasProperty(space_property))
            {
                space = mol.data().property(space_property).asA<Space>();
                break;
            }
        }
    }

    const auto coords = c.property<Vector>(map["coordinates"]);

    return FrameTransform(SireMaths::getAlignment(refcoords,
                                                  QVector<Vector>(coords.begin(), coords.end()),
                                                  true),
                          this->center(),
                          space.read(),
                          this->nSmooth(),
                          this->wrap());
}

QList<FrameTransform> TrajectoryAligner::operator[](const QList<qint64> &idxs) const
{
    QList<FrameTransform> ret;

    for (const auto &idx : idxs)
    {
        ret.append(this->operator[](idx));
    }

    return ret;
}

QList<FrameTransform> TrajectoryAligner::operator[](const SireBase::Slice &slice) const
{
    QList<FrameTransform> ret;

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
