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

#include <QMutex>
#include <limits>

#include "cartesian.h"
#include "space.h"
#include "transformedspace.h"

#include "SireMaths/align.h"
#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireVol;
using namespace SireBase;
using namespace SireStream;

///////////////
/////////////// Implementation of Space
///////////////

static const RegisterMetaType<Space> r_space(MAGIC_ONLY, "SireVol::Space");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Space &space)
{
    writeHeader(ds, r_space, 1) << static_cast<const Property &>(space);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Space &space)
{
    VersionID v = readHeader(ds, r_space);

    if (v == 1)
    {
        ds >> static_cast<Property &>(space);
    }
    else
        throw version_error(v, "1", r_space, CODELOC);

    return ds;
}

/** Construct a Space. */
Space::Space() : Property()
{
}

/** Copy constructor */
Space::Space(const Space &other) : Property(other)
{
}

/** Destructor */
Space::~Space()
{
}

/** Return the maximum cutoff that can be used in this space
    so that calculations obey the minimum image convention.
    This is the largest possible double for non-periodic spaces,
    and half the smallest box side for periodic spaces
*/
SireUnits::Dimension::Length Space::maximumCutoff() const
{
    return SireUnits::Dimension::Length(std::numeric_limits<double>::max());
}

/** Return the matrix that describes this box - these are the
 *  three box vectors. This raises an exception if this is
 *  not a periodic space
 */
Matrix Space::boxMatrix() const
{
    throw SireError::incompatible_error(QObject::tr(
                                            "The space '%1' is not a periodic space, and so does not have "
                                            "a box matrix")
                                            .arg(this->toString()),
                                        CODELOC);
}

/** Change the volume of this space by 'delta' */
SpacePtr Space::changeVolume(SireUnits::Dimension::Volume delta) const
{
    return this->setVolume(this->volume() + delta);
}

/** Return a random point within this space using the passed
    random number generator to generate any necessary random
    numbers, and centering the box at the origin */
Vector Space::getRandomPoint(const RanGenerator &generator) const
{
    return this->getRandomPoint(Vector(0, 0, 0), generator);
}

/** Return a random point within this space using the global
    random number generator and with the box centered at 'center' */
Vector Space::getRandomPoint(const Vector &center) const
{
    RanGenerator generator;
    return this->getRandomPoint(center, generator);
}

/** Return a random point within this space using the global
    random number generator, and with the box centered at the origin */
Vector Space::getRandomPoint() const
{
    RanGenerator generator;
    return this->getRandomPoint(Vector(0, 0, 0), generator);
}

double Space::minimumDistance(const Vector &point, const CoordGroup &group) const
{
    QVector<Vector> coords(1, point);
    return this->minimumDistance(CoordGroup(coords), group);
}

/** Assert that 'other' is of the same type as this space

    \throw SireError::incompatible_error
*/
void Space::assertCompatible(const Space &other) const
{
    if (QLatin1String(this->what()) != QLatin1String(other.what()))
        throw SireError::incompatible_error(QObject::tr("This space (of type \"%1\") is incompatible with "
                                                        "a space of type \"%2\".")
                                                .arg(this->what())
                                                .arg(other.what()),
                                            CODELOC);
}

/** Return the minimum image copy of 'coords' with respect to 'center',
    where the coordinates are "made whole". This means that they are
    translated as a single group, but the group as a whole will not
    be split across a periodic boundary. Use this function if you want
    to restore a molecule that has been split over a space into a single,
    coherent entity (all of the coordinates physically close to
    one another)
*/
QVector<Vector> Space::makeWhole(const QVector<Vector> &coords, const Vector &center) const
{
    if (coords.isEmpty() or (not this->isPeriodic()))
        return coords;
    else if (coords.count() == 1)
    {
        Vector image = this->getMinimumImage(coords[0], center);

        if (image != coords[0])
            // this has changed
            return QVector<Vector>({image});
        else
            return coords;
    }

    // build this up, getting the minimum image for each point as we
    // walk along the points
    QVector<Vector> mapped = coords;
    Vector *mapped_data = 0;

    const auto coords_data = coords.constData();
    const int n = coords.count();

    // the starting point in the center
    Vector ref_point = center;

    for (int i = 0; i < n; ++i)
    {
        auto image = this->getMinimumImage(coords_data[i], ref_point);

        if (image != coords_data[i])
        {
            if (mapped_data == 0)
                mapped_data = mapped.data();

            mapped_data[i] = image;
        }

        ref_point = image;
    }

    // this will only be different to `coords` if any
    // of the atoms have been changed
    return mapped;
}

/** Make the passed group of coordinates 'whole'. This will make sure
 *  that they are all next to each other, and aren't split across a
 *  periodic image boundary. The box that will be chosen will be the
 *  one that contains the center of the points, with the points mapped
 *  from the first to the last
 */
QVector<Vector> Space::makeWhole(const QVector<Vector> &coords) const
{
    if (coords.count() < 2 or (not this->isPeriodic()))
        return coords;

    // calculate the spatial central coordinates
    const Vector center = AABox(coords).center();

    return this->makeWhole(coords, center);
}

/** Return the minimum image copy of 'coords' with respect to 'center',
    where the coordinates are "made whole". This means that they are
    translated as a single group, but the group as a whole will not
    be split across a periodic boundary. Use this function if you want
    to restore a molecule that has been split over a space into a single,
    coherent entity (all of the coordinates physically close to
    one another). This treats all of the passed arrays
    of coordinates as a single unit that should not be split
*/
QVector<QVector<Vector>> Space::makeWhole(const QVector<QVector<Vector>> &coords,
                                          const Vector &center) const
{
    if (coords.isEmpty() or (not this->isPeriodic()))
        return coords;

    QVector<QVector<Vector>> new_coords(coords);
    QVector<Vector> *new_coords_data = 0;

    const QVector<Vector> *coords_data = coords.constData();

    const int n = coords.count();

    for (int i = 0; i < n; ++i)
    {
        const auto whole_coords = this->makeWhole(coords_data[i], center);

        if (whole_coords.constData() != coords_data[i].constData())
        {
            // the coordinates changed
            if (new_coords_data == 0)
            {
                new_coords_data = new_coords.data();
            }

            new_coords_data[i] = whole_coords;
        }
    }

    return new_coords;
}

/** Make the passed group of coordinates 'whole'. This will make sure
 *  that they are all next to each other, and aren't split across a
 *  periodic image boundary. The box that will be chosen will be the
 *  one that contains the center of the points, with the points mapped
 *  from the first to the last. This treats all of the passed arrays
 *  of coordinates as a single unit that should not be split
 */
QVector<QVector<Vector>> Space::makeWhole(const QVector<QVector<Vector>> &coords) const
{
    if (coords.isEmpty() or (not this->isPeriodic()))
        return coords;

    const Vector center = AABox(coords).center();

    return this->makeWhole(coords, center);
}

SpacePtr Space::transform(const Transform &tform, bool forwards) const
{
    if (forwards)
        return SpacePtr(new TransformedSpace(*this, tform));
    else
        return SpacePtr(new TransformedSpace(*this, tform.inverse()));
}

Q_GLOBAL_STATIC(Cartesian, nullCartesian)

/** Return the default space (Cartesian infinite box) */
const Cartesian &Space::null()
{
    return *(nullCartesian());
}
