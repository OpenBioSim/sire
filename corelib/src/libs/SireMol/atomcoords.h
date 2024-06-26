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

#ifndef SIREMOL_ATOMCOORDS_H
#define SIREMOL_ATOMCOORDS_H

#include "atomproperty.hpp"

#include "SireVol/coordgroup.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    template <>
    class AtomProperty<SireMaths::Vector>;
}

SIREMOL_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::AtomProperty<SireMaths::Vector> &);
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &, SireMol::AtomProperty<SireMaths::Vector> &);

namespace SireMaths
{
    class Vector;
    class Quaternion;
    class AxisSet;
    class Matrix;
    class Transform;
} // namespace SireMaths

namespace SireBase
{
    class Slice;
}

namespace SireVol
{
    class Space;
}

namespace SireMol
{

    using SireMaths::AxisSet;
    using SireMaths::Matrix;
    using SireMaths::Quaternion;
    using SireMaths::Transform;
    using SireMaths::Vector;

    using SireVol::CoordGroup;
    using SireVol::CoordGroupArray;
    using SireVol::Space;

    /** This is an explicit specialisation of AtomProperty<T> for the Vector
        class, as the Vector implies coordinates, which are arranged into
        CoordGroups (so that bounding boxes are calculated automatically)

        @author Christopher Woods
    */
    template <>
    class SIREMOL_EXPORT AtomProperty<Vector> : public SireBase::ConcreteProperty<AtomProperty<Vector>, AtomProp>
    {

        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const AtomProperty<SireMaths::Vector> &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, AtomProperty<SireMaths::Vector> &);

    public:
        AtomProperty();

        AtomProperty(const MoleculeInfoData &molinfo);

        AtomProperty(const CoordGroup &cgroup);
        AtomProperty(const CoordGroupArray &cgroups);

        AtomProperty(const AtomProperty<Vector> &other);

        ~AtomProperty();

        AtomProperty<Vector> &operator=(const AtomProperty<Vector> &other);

        static const char *typeName();

        AtomProperty<Vector> *clone() const;

        bool isEmpty() const;

        QString toString() const;

        bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

        bool operator==(const AtomProperty<Vector> &other) const;
        bool operator!=(const AtomProperty<Vector> &other) const;

        bool canConvert(const QVariant &value) const;

        void assignFrom(const AtomProperty<QVariant> &values);

        AtomProperty<QVariant> toVariant() const;

        static AtomProperty<Vector> fromVariant(const AtomProperty<QVariant> &variant);

        const CoordGroup &operator[](CGIdx cgidx) const;

        const CoordGroup &at(CGIdx cgidx) const;
        const CoordGroup &get(CGIdx cgidx) const;

        const Vector &operator[](int i) const;
        const Vector &operator[](const CGAtomIdx &cgatomidx) const;
        const Vector &at(int i) const;
        const Vector &at(const CGAtomIdx &cgatomidx) const;
        const Vector &get(int i) const;
        const Vector &get(const CGAtomIdx &cgatomidx) const;

        QList<Vector> operator[](const QList<qint64> &idxs) const;
        QList<Vector> operator[](const SireBase::Slice &slice) const;

        QVariant getAsVariant(const CGAtomIdx &cgatomidx) const;
        PropertyPtr getAsProperty(const CGAtomIdx &cgatomidx) const;

        AtomProperty<Vector> &set(const CGAtomIdx &cgatomidx, const Vector &value);

        AtomProperty<Vector> &set(CGIdx cgidx, const QVector<Vector> &values);
        AtomProperty<Vector> &set(CGIdx cgidx, const CoordGroup &cgroup);
        AtomProperty<Vector> &set(int i, const Vector &value);

        void translate(const Vector &delta);
        void translate(CGIdx cgidx, const Vector &delta);

        void rotate(const Quaternion &quat, const Vector &point);
        void rotate(const Matrix &rotmat, const Vector &point);

        void rotate(CGIdx cgidx, const Quaternion &quat, const Vector &point);
        void rotate(CGIdx cgidx, const Matrix &rotmat, const Vector &point);

        void transform(const Transform &t);
        void transform(CGIdx cgidx, const Transform &t);

        void mapInto(const AxisSet &axes);
        void mapInto(CGIdx cgidx, const AxisSet &axes);

        void changeFrame(const AxisSet &from_frame, const AxisSet &to_frame);
        void changeFrame(CGIdx cgidx, const AxisSet &from_frame, const AxisSet &to_frame);

        const CoordGroup *data() const;
        const CoordGroup *constData() const;

        const Vector *data(CGIdx cgidx) const;
        const Vector *constData(CGIdx cgidx) const;

        int size() const;
        int count() const;

        int nCutGroups() const;

        int nAtoms() const;
        int nAtoms(CGIdx cgidx) const;

        const CoordGroupArray &array() const;

        QVector<Vector> toVector() const;
        QVector<Vector> toVector(const AtomSelection &selection) const;

        QList<Vector> toList() const;
        QList<Vector> toList(const AtomSelection &selection) const;

        void copyFrom(const QVector<Vector> &values);
        void copyFrom(const QVector<Vector> &values, const AtomSelection &selection);

        PropertyPtr merge(const MoleculeInfoData &molinfo) const;
        PropertyPtr divide(const QVector<AtomSelection> &beads) const;
        PropertyPtr divideByResidue(const MoleculeInfoData &molinfo) const;

        void assertCanConvert(const QVariant &value) const;

        SireBase::PropertyList merge(const MolViewProperty &other,
                                     const AtomIdxMapping &mapping,
                                     const QString &ghost = QString(),
                                     const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    private:
        /** The actual atomic coordinates, arranged into CoordGroups */
        CoordGroupArray coords;
    };

    /** Return the raw array that hold the coordinates */
    SIRE_ALWAYS_INLINE const CoordGroupArray &AtomProperty<Vector>::array() const
    {
        return coords;
    }

    typedef AtomProperty<Vector> AtomCoords;

} // namespace SireMol

Q_DECLARE_METATYPE(SireMol::AtomCoords);

SIRE_EXPOSE_ATOM_PROPERTY(SireMaths::Vector, SireMol::AtomCoords)

SIRE_END_HEADER

#endif
