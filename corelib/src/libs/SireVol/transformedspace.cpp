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
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#include "transformedspace.h"
#include "coordgroup.h"

#include "SireBase/countflops.h"

#include "SireMaths/rangenerator.h"
#include "SireMaths/align.h"

#include "SireError/errors.h"
#include "SireVol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireVol;
using namespace SireBase;
using namespace SireMaths;
using namespace SireUnits::Dimension;
using namespace SireStream;

using boost::tuple;

static const RegisterMetaType<TransformedSpace> r_ts;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const TransformedSpace &ts)
{
    writeHeader(ds, r_ts, 1);

    SharedDataStream sds(ds);

    sds << ts.tform << ts.spc << static_cast<const Space &>(ts);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, TransformedSpace &ts)
{
    VersionID v = readHeader(ds, r_ts);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> ts.tform >> ts.spc >> static_cast<Space &>(ts);
        ts.inv_tform = ts.tform.inverse();
    }
    else
        throw version_error(v, "1", r_ts, CODELOC);

    return ds;
}

/** Construct a default TransformedSpace volume */
TransformedSpace::TransformedSpace() : ConcreteProperty<TransformedSpace, Space>()
{
}

/** Construct to represent 'space' transformed by the passed transform */
TransformedSpace::TransformedSpace(const Space &space, const Transform &transform)
    : ConcreteProperty<TransformedSpace, Space>(),
      tform(transform), inv_tform(transform.inverse()),
      spc(space)
{
}

/** Copy constructor */
TransformedSpace::TransformedSpace(const TransformedSpace &other)
    : ConcreteProperty<TransformedSpace, Space>(other),
      tform(other.tform), inv_tform(other.inv_tform), spc(other.spc)
{
}

/** Destructor */
TransformedSpace::~TransformedSpace()
{
}

/** Copy assignment operator */
TransformedSpace &TransformedSpace::operator=(const TransformedSpace &other)
{
    if (this != &other)
    {
        tform = other.tform;
        inv_tform = other.inv_tform;
        spc = other.spc;
        Space::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool TransformedSpace::operator==(const TransformedSpace &other) const
{
    return tform == other.tform and spc.read().equals(other.spc.read());
}

/** Comparison operator */
bool TransformedSpace::operator!=(const TransformedSpace &other) const
{
    return not this->operator==(other);
}

/** Return a string representation of this space */
QString TransformedSpace::toString() const
{
    return QObject::tr("TransformedSpace( %1 by %2 )")
        .arg(spc.read().toString())
        .arg(tform.toString());
}

bool TransformedSpace::isPeriodic() const
{
    return spc.read().isPeriodic();
}

bool TransformedSpace::isCartesian() const
{
    return spc.read().isCartesian();
}

Matrix TransformedSpace::boxMatrix() const
{
    // return the transformed value
    return tform.rotationMatrix() * spc.read().boxMatrix();
}

SireUnits::Dimension::Volume TransformedSpace::volume() const
{
    return spc.read().volume();
}

SpacePtr TransformedSpace::setVolume(SireUnits::Dimension::Volume volume) const
{
    TransformedSpace ret(*this);
    ret.spc = this->spc.read().setVolume(volume);
    return ret;
}

/** Calculate the distance between two points */
double TransformedSpace::calcDist(const Vector &point0, const Vector &point1) const
{
    return spc.read().calcDist(inv_tform(point0), inv_tform(point1));
}

/** Calculate the distance squared between two points */
double TransformedSpace::calcDist2(const Vector &point0, const Vector &point1) const
{
    return spc.read().calcDist2(inv_tform(point0), inv_tform(point1));
}

/** Populate the matrix 'mat' with the distances between all points in
    the group 'group'. Return the shortest distance between points. */
double TransformedSpace::calcDist(const CoordGroup &group, DistMatrix &mat) const
{
    CoordGroup t = group.edit().transform(inv_tform).commit();

    return spc.read().calcDist(t, mat);
}

/** Populate the matrix 'mat' with the distances^2 between all points in
    the group 'group'. Return the shortest distance between points. */
double TransformedSpace::calcDist2(const CoordGroup &group, DistMatrix &mat) const
{
    CoordGroup t = group.edit().transform(inv_tform).commit();

    return spc.read().calcDist2(t, mat);
}

/** Populate the matrix 'mat' with the inverse distances between all points in
    the group 'group'. Return the smallest distance between points. */
double TransformedSpace::calcInvDist(const CoordGroup &group, DistMatrix &mat) const
{
    CoordGroup t = group.edit().transform(inv_tform).commit();

    return spc.read().calcInvDist(t, mat);
}

/** Populate the matrix 'mat' with the inverse distances^2 between all points in
    the group 'group'. Return the smallest distance between points. */
double TransformedSpace::calcInvDist2(const CoordGroup &group, DistMatrix &mat) const
{
    CoordGroup t = group.edit().transform(inv_tform).commit();

    return spc.read().calcInvDist2(t, mat);
}

/** Populate the matrix 'mat' with the distances between all of the
    points of the two CoordGroups. Return the shortest distance between the two
    CoordGroups. */
double TransformedSpace::calcDist(const CoordGroup &group0, const CoordGroup &group1, DistMatrix &mat) const
{
    CoordGroup g0 = group0.edit().transform(inv_tform).commit();
    CoordGroup g1 = group1.edit().transform(inv_tform).commit();

    return spc.read().calcDist(g0, g1, mat);
}

/** Populate the matrix 'mat' with the distances between all of the
    points of the passed CoordGroup and 'point'. Returns the shortest distance. */
double TransformedSpace::calcDist(const CoordGroup &group, const Vector &point, DistMatrix &mat) const
{
    CoordGroup g = group.edit().transform(inv_tform).commit();

    return spc.read().calcDist(g, inv_tform(point), mat);
}

/** Populate the matrix 'mat' with the distances^2 between all of the
    points of the two CoordGroups. Return the shortest distance between the
    two CoordGroups. */
double TransformedSpace::calcDist2(const CoordGroup &group0, const CoordGroup &group1, DistMatrix &mat) const
{
    CoordGroup g0 = group0.edit().transform(inv_tform).commit();
    CoordGroup g1 = group1.edit().transform(inv_tform).commit();

    return spc.read().calcDist2(g0, g1, mat);
}

/** Populate the matrix 'mat' with the distances between all of the
    points of the passed CoordGroup and 'point'. Returns the shortest distance. */
double TransformedSpace::calcDist2(const CoordGroup &group, const Vector &point, DistMatrix &mat) const
{
    CoordGroup g = group.edit().transform(inv_tform).commit();

    return spc.read().calcDist2(g, inv_tform(point), mat);
}

/** Populate the matrix 'mat' with the inverse distances between all of the
    points of the two CoordGroups. Return the shortest distance between
    the two CoordGroups. */
double TransformedSpace::calcInvDist(const CoordGroup &group0, const CoordGroup &group1, DistMatrix &mat) const
{
    CoordGroup g0 = group0.edit().transform(inv_tform).commit();
    CoordGroup g1 = group1.edit().transform(inv_tform).commit();

    return spc.read().calcInvDist(g0, g1, mat);
}

/** Populate the matrix 'mat' with the inverse distances^2 between all of the
    points of the two CoordGroups. Return the shortest distance between
    the two CoordGroups. */
double TransformedSpace::calcInvDist2(const CoordGroup &group0, const CoordGroup &group1, DistMatrix &mat) const
{
    CoordGroup g0 = group0.edit().transform(inv_tform).commit();
    CoordGroup g1 = group1.edit().transform(inv_tform).commit();

    return spc.read().calcInvDist2(g0, g1, mat);
}

/** Calculate the distance vector between two points */
DistVector TransformedSpace::calcDistVector(const Vector &point0, const Vector &point1) const
{
    return spc.read().calcDistVector(inv_tform(point0), inv_tform(point1));
}

/** Populate the matrix 'distmat' with all of the interpoint distance vectors
    between all points within the CoordGroup. This is *not* a symmetrical matrix,
    as the direction from point A to point B is the negative of the
    direction from point B to point A. This returns the shortest distance
    between two points in the group (that is not the self-self distance) */
double TransformedSpace::calcDistVectors(const CoordGroup &group, DistVectorMatrix &mat) const
{
    CoordGroup g = group.edit().transform(inv_tform).commit();

    return spc.read().calcDistVectors(g, mat);
}

/** Populate the matrix 'distmat' between all the points of the two CoordGroups
    'group1' and 'group2' - the returned matrix has the vectors pointing
    from each point in 'group1' to each point in 'group2'. This returns
    the shortest distance between two points in the group */
double TransformedSpace::calcDistVectors(const CoordGroup &group0, const CoordGroup &group1, DistVectorMatrix &mat) const
{
    CoordGroup g0 = group0.edit().transform(inv_tform).commit();
    CoordGroup g1 = group1.edit().transform(inv_tform).commit();

    return spc.read().calcDistVectors(g0, g1, mat);
}

/** Populate the matrix 'distmat' between all the points of the passed
    CoordGroup with 'point' - the returned matrix has the vectors pointing
    from the point, to each point in 'group'. This returns the shortest distance. */
double TransformedSpace::calcDistVectors(const CoordGroup &group, const Vector &point, DistVectorMatrix &mat) const
{
    CoordGroup g = group.edit().transform(inv_tform).commit();

    return spc.read().calcDistVectors(g, inv_tform(point), mat);
}

/** Calculate the angle between the passed three points. This should return
    the acute angle between the points, which should lie between 0 and 180 degrees */
Angle TransformedSpace::calcAngle(const Vector &point0, const Vector &point1, const Vector &point2) const
{
    return spc.read().calcAngle(inv_tform(point0), inv_tform(point1), inv_tform(point2));
}

/** Calculate the torsion angle between the passed four points. This should
    return the torsion angle measured clockwise when looking down the
    torsion from point0-point1-point2-point3. This will lie between 0 and 360
    degrees */
Angle TransformedSpace::calcDihedral(const Vector &point0, const Vector &point1, const Vector &point2,
                                     const Vector &point3) const
{
    return spc.read().calcDihedral(inv_tform(point0), inv_tform(point1),
                                   inv_tform(point2), inv_tform(point3));
}

/** Return whether or not two groups enclosed by the AABoxes 'aabox0'
    and 'aabox1' are definitely beyond the cutoff distance */
bool TransformedSpace::beyond(double dist, const AABox &aabox0, const AABox &aabox1) const
{
    AABox box0(aabox0);
    AABox box1(aabox1);

    box0.transform(inv_tform);
    box1.transform(inv_tform);

    return spc.read().beyond(dist, box0, box1);
}

/** Return whether or not these two groups are definitely beyond the cutoff distance. */
bool TransformedSpace::beyond(double dist, const CoordGroup &group0, const CoordGroup &group1) const
{
    return TransformedSpace::beyond(dist, group0.aaBox(), group1.aaBox());
}

/** Return the minimum distance between the points in 'group0' and 'group1'. */
double TransformedSpace::minimumDistance(const CoordGroup &group0, const CoordGroup &group1) const
{
    CoordGroup g0 = group0.edit().transform(inv_tform).commit();
    CoordGroup g1 = group1.edit().transform(inv_tform).commit();

    return spc.read().minimumDistance(g0, g1);
}

/** Return the minimum distance between the two passed boxes */
double TransformedSpace::minimumDistance(const AABox &box0, const AABox &box1) const
{
    AABox b0(box0);
    AABox b1(box1);

    b0.transform(inv_tform);
    b1.transform(inv_tform);

    return spc.read().minimumDistance(b0, b1);
}

/** Return the minimum distance between a point and a box */
double TransformedSpace::minimumDistance(const Vector &point, const AABox &box) const
{
    AABox tbox(box);
    tbox.transform(inv_tform);

    return spc.read().minimumDistance(inv_tform(point), tbox);
}

/** Return the minimum distance between points within the group 'group'. */
double TransformedSpace::minimumDistance(const CoordGroup &group) const
{
    CoordGroup g = group.edit().transform(inv_tform).commit();
    return spc.read().minimumDistance(g);
}

QVector<Vector> TransformedSpace::getMinimumImage(const QVector<Vector> &coords,
                                                  const Vector &point) const
{
    return tform(spc.read().getMinimumImage(inv_tform(coords), inv_tform(point)));
}

CoordGroup TransformedSpace::getMinimumImage(const CoordGroup &group, const Vector &point) const
{
    CoordGroup g = group.edit().transform(inv_tform).commit();

    g = spc.read().getMinimumImage(group, inv_tform(point));

    return g.edit().transform(tform).commit();
}

CoordGroupArray TransformedSpace::getMinimumImage(const CoordGroupArray &groups,
                                                  const Vector &point, bool translate_as_one) const
{
    if (translate_as_one)
    {
        QVector<CoordGroup> tgroups;

        for (int i = 0; i < groups.count(); ++i)
        {
            tgroups.append(groups.at(i).edit().transform(inv_tform).commit());
        }

        CoordGroupArray tarray(tgroups);

        tarray = spc.read().getMinimumImage(tarray, inv_tform(point), true);

        tgroups.clear();

        for (int i = 0; i < tarray.count(); ++i)
        {
            tgroups.append(tarray.at(i).edit().transform(tform).commit());
        }

        return CoordGroupArray(tgroups);
    }
    else
    {
        QVector<CoordGroup> tgroups;

        for (int i = 0; i < groups.count(); ++i)
        {
            tgroups.append(this->getMinimumImage(groups.at(i), point));
        }

        return CoordGroupArray(tgroups);
    }
}

/** A TransformedSpace space is not periodic, so this just returns the input 'aabox' */
AABox TransformedSpace::getMinimumImage(const AABox &aabox, const Vector &point) const
{
    AABox box(aabox);
    box.transform(inv_tform);
    return spc.read().getMinimumImage(box, inv_tform(point));
}

/** A TransformedSpace space is not periodic, so this just returns the input 'point' */
Vector TransformedSpace::getMinimumImage(const Vector &point, const Vector &center) const
{
    return spc.read().getMinimumImage(inv_tform(point), inv_tform(center));
}

/** Return all periodic images of 'point' with respect to 'center' within
    'dist' distance of 'center' */
QVector<Vector> TransformedSpace::getImagesWithin(const Vector &point, const Vector &center, double dist) const
{
    return tform(spc.read().getImagesWithin(inv_tform(point), inv_tform(center), dist));
}

/** Return a list of copies of CoordGroup 'group' that are within
    'distance' of the CoordGroup 'center', translating 'group' so that
    it has the right coordinates to be around 'center'. As this is not
    a periodic space, this will merely return a copy of 'group' if
    it is within the specified distance. */
QList<tuple<double, CoordGroup>> TransformedSpace::getCopiesWithin(const CoordGroup &group, const CoordGroup &center,
                                                                   double dist) const
{
    CoordGroup g = group.edit().transform(inv_tform).commit();
    CoordGroup c = center.edit().transform(inv_tform).commit();

    auto result = spc.read().getCopiesWithin(g, c, dist);

    for (auto &cgroup : result)
    {
        cgroup = tuple<double, CoordGroup>(cgroup.get<0>(),
                                           cgroup.get<1>().edit().transform(tform).commit());
    }

    return result;
}

/** Return a random point in this space - this can be truly anywhere!!!

    (well, it is limited to within -10^20 and 10^20 angstroms)
*/
Vector TransformedSpace::getRandomPoint(const Vector &center, const RanGenerator &generator) const
{
    return tform(spc.read().getRandomPoint(center, generator));
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at the origin */
Vector TransformedSpace::getBoxCenter(const Vector &point) const
{
    return tform(spc.read().getBoxCenter(inv_tform(point)));
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at 'center' */
Vector TransformedSpace::getBoxCenter(const Vector &p, const Vector &center) const
{
    return tform(spc.read().getBoxCenter(inv_tform(p), inv_tform(center)));
}

const char *TransformedSpace::typeName()
{
    return QMetaType::typeName(qMetaTypeId<TransformedSpace>());
}
