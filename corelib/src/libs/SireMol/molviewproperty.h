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

#ifndef SIREMOL_MOLVIEWPROPERTY_H
#define SIREMOL_MOLVIEWPROPERTY_H

#include <QHash>

#include "SireBase/property.h"
#include "SireBase/propertylist.h"

#include "moleculeinfo.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class MolViewProperty;
}

namespace SireMol
{

    class MoleculeInfoData;
    class MoleculeInfo;
    class MoleculeView;
    class AtomSelection;
    class AtomIdx;
    class AtomMatcher;
    class AtomIdxMapping;

    /** This is the base class of all properties that are specifically
        attached to views of a molecule (e.g. AtomProperty, ResProperty,
        SegProperty etc.)

        @author Christopher Woods
    */
    class SIREMOL_EXPORT MolViewProperty : public SireBase::Property
    {

    public:
        MolViewProperty();

        MolViewProperty(const MolViewProperty &other);

        ~MolViewProperty();

        static const char *typeName()
        {
            return "SireMol::MolViewProperty";
        }

        virtual bool isCompatibleWith(const MoleculeInfoData &molinfo) const = 0;
        virtual bool isCompatibleWith(const MoleculeInfo &molinfo) const;

        SireBase::PropertyPtr makeCompatibleWith(const MoleculeInfoData &molinfo, const AtomMatcher &atommatcher) const;
        SireBase::PropertyPtr makeCompatibleWith(const MoleculeInfoData &molinfo, const QHash<AtomIdx, AtomIdx> &map) const;

        SireBase::PropertyPtr makeCompatibleWith(const MoleculeInfoData &molinfo) const;

        SireBase::PropertyPtr makeCompatibleWith(const MoleculeView &mol) const;
        SireBase::PropertyPtr makeCompatibleWith(const MoleculeView &mol, const AtomMatcher &atommatcher) const;
        SireBase::PropertyPtr makeCompatibleWith(const MoleculeView &mol, const QHash<AtomIdx, AtomIdx> &map) const;

        virtual SireBase::PropertyList merge(const MolViewProperty &other,
                                             const AtomIdxMapping &mapping,
                                             const QString &ghost = QString(),
                                             const SireBase::PropertyMap &map = SireBase::PropertyMap()) const = 0;

        void assertCompatibleWith(const MoleculeInfoData &molinfo) const;
        void assertCompatibleWith(const MoleculeInfo &molinfo) const;

    protected:
        virtual SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                              const AtomMatcher &atommatcher) const;
        virtual SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                              const QHash<AtomIdx, AtomIdx> &map) const;
    };

    /** This specifically is a property that pertains to an entire
        molecule (e.g. a selection of atoms from that molecule)

        @author Christopher Woods
    */
    class SIREMOL_EXPORT MoleculeProperty : public MolViewProperty
    {
    public:
        MoleculeProperty();

        MoleculeProperty(const MoleculeProperty &other);

        ~MoleculeProperty();

        static const char *typeName()
        {
            return "SireMol::MoleculeProperty";
        }
    };

    SIRE_EXPOSE_CLASS(SireMol::MolViewProperty)
    SIRE_EXPOSE_CLASS(SireMol::MoleculeProperty)

} // namespace SireMol

SIRE_END_HEADER

#endif
