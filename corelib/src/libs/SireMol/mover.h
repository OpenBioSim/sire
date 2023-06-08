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

#ifndef SIREMOL_MOVER_H
#define SIREMOL_MOVER_H

#include "atomselection.h"
#include "moleculeview.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
    class AxisSet;
    class Transform;
} // namespace SireMaths

namespace SireVol
{
    class Space;
}

namespace SireMol
{

    template <>
    class AtomProperty<SireMaths::Vector>;

    typedef AtomProperty<SireMaths::Vector> AtomCoords;

    class BondID;
    class AngleID;
    class DihedralID;
    class ImproperID;

    using SireMaths::AxisSet;
    using SireMaths::Matrix;
    using SireMaths::Quaternion;
    using SireMaths::Transform;
    using SireMaths::Vector;

    using SireBase::PropertyMap;

    using SireVol::Space;

    SIREMOL_EXPORT Transform getAlignment(const MoleculeView &view0, const MoleculeView &view1, bool fit = true);
    SIREMOL_EXPORT Transform getAlignment(const MoleculeView &view0, const MoleculeView &view1, const PropertyMap &map,
                                          bool fit = true);
    SIREMOL_EXPORT Transform getAlignment(const MoleculeView &view0, const PropertyMap &map0, const MoleculeView &view1,
                                          const PropertyMap &map1, bool fit = true);
    SIREMOL_EXPORT Transform getAlignment(const MoleculeView &view0, const MoleculeView &view1, const AtomMatcher &matcher,
                                          bool fit = true);
    SIREMOL_EXPORT Transform getAlignment(const MoleculeView &view0, const MoleculeView &view1, const AtomMatcher &matcher,
                                          const PropertyMap &map, bool fit = true);
    SIREMOL_EXPORT Transform getAlignment(const MoleculeView &view0, const PropertyMap &map0, const MoleculeView &view1,
                                          const PropertyMap &map1, const AtomMatcher &matcher, bool fit = true);

    /** This class provides the template-independent part
        of Mover<T>. This class is not designed to be used
        on its own!

        @author Christopher Woods
    */
    class SIREMOL_EXPORT MoverBase
    {
    public:
        MoverBase();

        MoverBase(const MoverBase &other);

        ~MoverBase();

        MoverBase &operator=(const MoverBase &other);

    protected:
        MoverBase(const AtomSelection &selected_atom);

        void setMovableAtoms(const AtomSelection &movable_atoms);

        void makeWhole(MoleculeData &view, const PropertyMap &map) const;
        void makeWhole(MoleculeData &view, const SireMaths::Vector &center, const PropertyMap &map) const;

        void translate(MoleculeData &data, const SireMaths::Vector &delta, const PropertyMap &map) const;

        void rotate(MoleculeData &data, const Quaternion &quat, const SireMaths::Vector &point, const PropertyMap &map) const;

        void rotate(MoleculeData &data, const Matrix &rotmat, const SireMaths::Vector &point, const PropertyMap &map) const;

        void transform(MoleculeData &data, const Transform &transform, const PropertyMap &map) const;

        void mapInto(MoleculeData &data, const AxisSet &axes, const PropertyMap &map) const;

        void changeFrame(MoleculeData &data, const AxisSet &from_frame, const AxisSet &to_frame,
                         const PropertyMap &map) const;

        void change(MoleculeData &data, const BondID &bond, SireUnits::Dimension::Length delta,
                    const PropertyMap &map) const;

        void change(MoleculeData &data, const AngleID &angle, SireUnits::Dimension::Angle delta,
                    const PropertyMap &map) const;

        void change(MoleculeData &data, const DihedralID &dihedral, SireUnits::Dimension::Angle delta,
                    const PropertyMap &map) const;

        void change(MoleculeData &data, const BondID &bond, SireUnits::Dimension::Angle delta,
                    const PropertyMap &map) const;

        void change(MoleculeData &data, const ImproperID &improper, SireUnits::Dimension::Angle delta,
                    const PropertyMap &map) const;

        void set(MoleculeData &data, const BondID &bond, SireUnits::Dimension::Length value, const PropertyMap &map) const;

        void set(MoleculeData &data, const AngleID &angle, SireUnits::Dimension::Angle value, const PropertyMap &map) const;

        void set(MoleculeData &data, const DihedralID &dihedral, SireUnits::Dimension::Angle value,
                 const PropertyMap &map) const;

        void setAll(MoleculeData &data, const DihedralID &dihedral, SireUnits::Dimension::Angle value,
                    const PropertyMap &map) const;

        void set(MoleculeData &data, const ImproperID &improper, SireUnits::Dimension::Angle value,
                 const PropertyMap &map) const;

        static void mapInto(MoleculeData &data, const AtomSelection &selected_atoms, const AxisSet &axes,
                            const PropertyMap &map);

        static void changeFrame(MoleculeData &data, const AtomSelection &selected_atoms, const AxisSet &from_axes,
                                const AxisSet &to_axes, const PropertyMap &map);

        static void makeWhole(MoleculeData &view, const AtomSelection &selected_atoms,
                              const PropertyMap &map);

        static void makeWhole(MoleculeData &view, const AtomSelection &selected_atoms,
                              const SireMaths::Vector &center, const PropertyMap &map);

        static void translate(MoleculeData &view, const AtomSelection &selected_atoms, const SireMaths::Vector &delta,
                              const PropertyMap &map);

        static void rotate(MoleculeData &view, const AtomSelection &selected_atoms, const Matrix &rotmat,
                           const SireMaths::Vector &point, const PropertyMap &map);

        static void transform(MoleculeData &view, const AtomSelection &selected_atoms, const Transform &t,
                              const PropertyMap &map);

        static void mapInto(AtomCoords &coords, const AtomSelection &selected_atoms, const AxisSet &axes);

        static void changeFrame(AtomCoords &coords, const AtomSelection &selected_atoms, const AxisSet &from_frame,
                                const AxisSet &to_frame);

        static void makeWhole(AtomCoords &coords, const AtomSelection &selected_atoms,
                              const SireVol::Space &space);
        static void makeWhole(AtomCoords &coords, const AtomSelection &selected_atoms,
                              const SireMaths::Vector &center, const SireVol::Space &space);

        static void translate(AtomCoords &coords, const AtomSelection &selected_atoms, const SireMaths::Vector &delta);

        static void rotate(AtomCoords &coords, const AtomSelection &selected_atoms, const Matrix &rotmat,
                           const SireMaths::Vector &point);

        static void rotate(AtomCoords &coords, const AtomSelection &selected_atoms, const Quaternion &quat,
                           const SireMaths::Vector &point);

        static void transform(AtomCoords &coords, const AtomSelection &selected_atoms, const Transform &t);

        /** The only atoms that can be moved by this Mover */
        AtomSelection movable_atoms;
    };

} // namespace SireMol

SIRE_EXPOSE_CLASS(SireMol::MoverBase)
SIRE_EXPOSE_FUNCTION(SireMol::getAlignment)

SIRE_END_HEADER

#endif
