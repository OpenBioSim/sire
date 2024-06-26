/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREVOL_TRICLINICBOX_H
#define SIREVOL_TRICLINICBOX_H

#include "cartesian.h"

#include "SireMaths/matrix.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
    class TriclinicBox;
}

SIREVOL_EXPORT QDataStream &operator<<(QDataStream &, const SireVol::TriclinicBox &);
SIREVOL_EXPORT QDataStream &operator>>(QDataStream &, SireVol::TriclinicBox &);

namespace SireVol
{

    using SireMaths::Matrix;
    using SireMaths::Vector;

    /**
    A TriclinicBox is a volume that represents standard periodic boundary conditions
    (a 3D box replicated to infinity along all three dimensions).

    To support triclinic boxes that work across a range of molecular simulation
    engines, e.g. AMBER, GROMACS, OpenMM, we represent the triclinic space in
    reduced form, using the approach documented in Appendix A of Chapter 3 from
    "Molecular dynamics of sense and sensibility in processing and analysis of data"
    by Tsjerk A. Wassenaar.

    @author Lester Hedges
    */
    class SIREVOL_EXPORT TriclinicBox : public SireBase::ConcreteProperty<TriclinicBox, Cartesian>
    {

        friend SIREVOL_EXPORT QDataStream & ::operator<<(QDataStream &, const TriclinicBox &);
        friend SIREVOL_EXPORT QDataStream & ::operator>>(QDataStream &, TriclinicBox &);

    public:
        TriclinicBox();

        /** Construct a triclinic box from lattice vectors.

            @param v0
                The first lattice vector.

            @param v1
                The second lattice vector.

            @param v2
                The third lattice vector.

            @param auto_rotate
                Whether to automatically rotate the box to comply with the
                constraints of molecular dynamics engines, i.e. vector0 aligned
                with x axis, vector1 in x-y plane, and vector2 with positive
                z component.

            @param auto_reduce
                Whether to automatically perform a lattice reduction on the
                box.
         */
        TriclinicBox(
            const Vector &v0,
            const Vector &v1,
            const Vector &v2,
            bool auto_rotate=false,
            bool auto_reduce=false
        );

        /** Construct a triclinic box from box lengths and angles.

            @param a
                The length of the first box vector.

            @param b
                The length of the second box vector.

            @param c
                The length of the third box vector.

            @param alpha
                The angle between the second and third box vectors.

            @param beta
                The angle between the first and third box vectors.

            @param gamma
                The angle between the second and first box vectors.

            @param auto_rotate
                Whether to automatically rotate the box to comply with the
                constraints of molecular dynamics engines, i.e. vector0 aligned
                with x axis, vector1 in x-y plane, and vector2 with positive
                z component:w.

            @param auto_reduce
                Whether to automatically perform a lattice reduction on the
                box.
         */
        TriclinicBox(
            double a,
            double b,
            double c,
            const SireUnits::Dimension::Angle &alpha,
            const SireUnits::Dimension::Angle &beta,
            const SireUnits::Dimension::Angle &gamma,
            bool auto_rotate=false,
            bool auto_reduce=false
        );

        TriclinicBox(const TriclinicBox &other);

        ~TriclinicBox();

        TriclinicBox &operator=(const TriclinicBox &other);

        bool operator==(const TriclinicBox &other) const;
        bool operator!=(const TriclinicBox &other) const;

        bool isPeriodic() const;
        bool isCartesian() const;

        QString toString() const;

        /** Get the maximum cutoff distance for the triclinic box. */
        SireUnits::Dimension::Length maximumCutoff() const;

        /** Get the volume of the triclinic box. */
        SireUnits::Dimension::Volume volume() const;

        /** Set the volume of the triclinic box. */
        SpacePtr setVolume(SireUnits::Dimension::Volume volume) const;

        Matrix boxMatrix() const;

        static const char *typeName();

        /** Whether the triclinic cell has been rotated to comply with the constraints
            of molecular dynamics engines, i.e. vector0 aligned with x axis, vector1
            in x-y plane, and vector2 with positive z component.
         */
        bool isRotated() const;

        /** Whether an automatic lattice reduction has been performed. */
        bool isReduced() const;

        double calcDist(const Vector &point0, const Vector &point1) const;

        double calcDist2(const Vector &point0, const Vector &point1) const;

        double calcDist(const CoordGroup &group1, const CoordGroup &group2, DistMatrix &distmat) const;

        double calcDist(const CoordGroup &group, const Vector &point, DistMatrix &mat) const;

        double calcDist2(const CoordGroup &group, const Vector &point, DistMatrix &mat) const;

        double calcDist2(const CoordGroup &group1, const CoordGroup &group2, DistMatrix &distmat) const;

        double calcInvDist(const CoordGroup &group1, const CoordGroup &group2, DistMatrix &distmat) const;

        double calcInvDist2(const CoordGroup &group1, const CoordGroup &group2, DistMatrix &distmat) const;

        DistVector calcDistVector(const Vector &point0, const Vector &point1) const;

        double calcDistVectors(const CoordGroup &group1, const CoordGroup &group2, DistVectorMatrix &distmat) const;

        double calcDistVectors(const CoordGroup &group, const Vector &point, DistVectorMatrix &distmat) const;

        SireUnits::Dimension::Angle calcAngle(const Vector &point0, const Vector &point1, const Vector &point2) const;

        SireUnits::Dimension::Angle calcDihedral(const Vector &point0, const Vector &point1, const Vector &point2,
                                                 const Vector &point3) const;

        bool beyond(double dist, const AABox &aabox0, const AABox &aabox1) const;

        bool beyond(double dist, const CoordGroup &group0, const CoordGroup &group1) const;

        double minimumDistance(const CoordGroup &group0, const CoordGroup &group1) const;

        double minimumDistance(const AABox &box0, const AABox &box1) const;

        QVector<Vector> getMinimumImage(const QVector<Vector> &coords,
                                        const Vector &center) const;

        Vector getMinimumImage(const Vector &point, const Vector &center) const;

        CoordGroup getMinimumImage(const CoordGroup &group, const Vector &center) const;

        CoordGroupArray getMinimumImage(const CoordGroupArray &groups, const Vector &center,
                                        bool translate_as_one = false) const;

        AABox getMinimumImage(const AABox &aabox, const Vector &center) const;

        QVector<Vector> getImagesWithin(const Vector &point, const Vector &center, double dist) const;

        QList<boost::tuple<double, CoordGroup>> getCopiesWithin(const CoordGroup &group, const CoordGroup &center,
                                                                double dist) const;

        Vector getRandomPoint(const Vector &center, const RanGenerator &generator) const;

        Vector getBoxCenter(const Vector &p) const;
        Vector getBoxCenter(const Vector &p, const Vector &center) const;

        /** Return the first box vector. */
        const Vector &vector0() const;

        /** Return the second box vector. */
        const Vector &vector1() const;

        /** Return the third box vector. */
        const Vector &vector2() const;

        /** Return the angle between v1 and v2 in degrees. */
        double alpha() const;

        /** Return the angle between v0 and v2 in degrees. */
        double beta() const;

        /** Return the angle between v1 and v0 in degrees. */
        double gamma() const;

        /** Return the rotation matrix. */
        const Matrix &rotationMatrix() const;

        /** Return the cell matrix. */
        Matrix cellMatrix() const;

        /** Return a cubic TriclinicBox with image distance d. */
        static TriclinicBox cubic(double d);

        /** Return a square rhombic dodecahedron TriclinicBox with image distance d. */
        static TriclinicBox rhombicDodecahedronSquare(double d, bool auto_rotate=true, bool auto_reduce=true);

        /** Return a hexagonal rhombic dodecahedron TriclinicBox with image distance d. */
        static TriclinicBox rhombicDodecahedronHexagon(double d, bool auto_rotate=true, bool auto_reduce=true);

        /** Return a truncated octahedron with image distance d. */
        static TriclinicBox truncatedOctahedron(double d, bool auto_rotate=true, bool auto_reduce=true);

        /** Rotate the triclinic cell to comply with the constraints of certain
            molecular dynamics engines, i.e. such that vector0 is aligned with
            the x axis, vector1, lies in the x-y plane, and vector2 has a positive
            z component.

            @param precision
                The precision to use when sorting the lattice vectors based on
                their magnitude. This can be used to prevent unwanted rotation
                when using input fixed-precision ascii molecular input files.
         */
        void rotate(double precision=0.0);

        /** Perform a lattice reduction on the triclinic cell.

            @param bias
                The bias to use when rounding during the lattice reduction.
                Negative values biases towards left-tilting boxes, whereas
                positive values biases towards right-tilting boxes. This can
                be used to ensure that rounding is performed in a consistent
                direction, avoiding oscillation when the TriclinicBox is
                instantiated from box vectors, or dimensions and angles, that
                have been read from fixed-precision input files.
         */
        void reduce(double bias=0.0);

    protected:
        void construct(const Vector &v0, const Vector &v1, const Vector &v2, bool auto_rotate=false, bool auto_reduce=false);

        void setAttributes();

        Vector wrapDelta(const Vector &v0, const Vector &v1) const;

        CoordGroupArray _pvt_getMinimumImage(const CoordGroupArray &groups, const Vector &point) const;

        /** The first box vector */
        Vector v0;

        /** The second box vector */
        Vector v1;

        /** The third box vector */
        Vector v2;

        /** The rotation matrix used to transform the box to meet the requirements
            of molecular dynamics engines.
          */
        Matrix rotation_matrix;

        /** The cell matrix. */
        Matrix cell_matrix;

        /** The inverse of the cell matrix. */
        Matrix cell_matrix_inverse;

        /** The matrix product of cell_matrix_inverse and cell_matrix. */
        Matrix M;

        /** The maximum distance within which a point will always be closer to the
            origin than any of its images.
          */
        double dist_max;

        /** The maximum axis length of the cell. */
        double max_length;

        /** The angle between vectors v1 and v2. */
        double _alpha;

        /** The angle between vectors v0 and v2. */
        double _beta;

        /** The angle between vectors v1 and v2. */
        double _gamma;

        /** The volume of the triclinic cell. */
        double vol;

        /** Whether the triclinic cell has been rotated. */
        bool is_rotated;

        /** Whether the triclinic cell has been reduced. */
        bool is_reduced;

        /** The inverse of the lengths of each side of the box */
        Vector invlength;
    };

} // namespace SireVol

Q_DECLARE_METATYPE(SireVol::TriclinicBox)

SIRE_EXPOSE_CLASS(SireVol::TriclinicBox)

SIRE_END_HEADER

#endif
