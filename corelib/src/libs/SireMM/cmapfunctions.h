/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2025  Christopher Woods
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

#ifndef SIREMM_CMAPFUNCTIONS_H
#define SIREMM_CMAPFUNCTIONS_H

#include "SireMol/molviewproperty.h"
#include "SireMol/atomid.h"
#include "SireMol/atomidx.h"
#include "SireMol/cgatomidx.h"

#include "cmapparameter.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class CMAPFunction;
    class CMAPFunctions;
} // namespace SireMM

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::CMAPFunction &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::CMAPFunction &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::CMAPFunctions &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::CMAPFunctions &);

namespace SireMol
{
    class AtomSelection;
}

namespace SireMM
{
    using SireMol::AtomID;
    using SireMol::AtomIdx;
    using SireMol::CGAtomIdx;

    /** This class holds a single CMAP function for a single set
     *  of five atoms
     */
    class SIREMM_EXPORT CMAPFunction
    {

        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const CMAPFunction &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, CMAPFunction &);

    public:
        CMAPFunction();
        CMAPFunction(const CGAtomIdx &atm0, const CGAtomIdx &atm1, const CGAtomIdx &atm2,
                     const CGAtomIdx &atm3, const CGAtomIdx &atm4, const CMAPParameter &param);

        CMAPFunction(const CMAPFunction &other);
        ~CMAPFunction();

        CMAPFunction &operator=(const CMAPFunction &other);

        bool operator==(const CMAPFunction &other) const;
        bool operator!=(const CMAPFunction &other) const;

        static const char *typeName();
        const char *what() const;

        QString toString() const;

        CMAPFunction *clone() const;

        const CGAtomIdx &atom0() const;
        const CGAtomIdx &atom1() const;
        const CGAtomIdx &atom2() const;
        const CGAtomIdx &atom3() const;
        const CGAtomIdx &atom4() const;

        const CMAPParameter &parameter() const;

    private:
        /** The indicies of the five atoms */
        CGAtomIdx atm0;
        CGAtomIdx atm1;
        CGAtomIdx atm2;
        CGAtomIdx atm3;
        CGAtomIdx atm4;

        /** The parameter for the CMAP function */
        CMAPParameter param;
    };

    namespace detail
    {
        class IDQuint
        {
        public:
            IDQuint(quint32 atm0 = 0, quint32 atm1 = 0, quint32 atm2 = 0,
                    quint32 atm3 = 0, quint32 atm4 = 0);

            IDQuint(const IDQuint &other);

            ~IDQuint();

            IDQuint &operator=(const IDQuint &other);

            bool operator==(const IDQuint &other) const;
            bool operator!=(const IDQuint &other) const;

            bool operator<(const IDQuint &other) const;
            bool operator<=(const IDQuint &other) const;
            bool operator>(const IDQuint &other) const;
            bool operator>=(const IDQuint &other) const;

            quint32 operator[](int i) const
            {
                switch (i)
                {
                case 0:
                    return atom0;
                case 1:
                    return atom1;
                case 2:
                    return atom2;
                case 3:
                    return atom3;
                default:
                    return atom4;
                }
            }

            quint32 atom0;
            quint32 atom1;
            quint32 atom2;
            quint32 atom3;
            quint32 atom4;
        };

        SIRE_ALWAYS_INLINE uint qHash(const IDQuint &idquint)
        {
            return (idquint.atom0 << 28) | ((idquint.atom1 << 24) & 0x0F000000) |
                   ((idquint.atom2 << 16) & 0x00FF0000) |
                   ((idquint.atom3 << 8) & 0x0000FF00) |
                   (idquint.atom4 & 0x000000FF);
        }
    }

    /** This class holds all of the CMAP parameters for a single molecule */
    class SIREMM_EXPORT CMAPFunctions : public SireBase::ConcreteProperty<CMAPFunctions, SireMol::MoleculeProperty>
    {
        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const CMAPFunctions &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, CMAPFunctions &);

    public:
        CMAPFunctions();

        CMAPFunctions(const SireMol::MoleculeData &moldata);
        CMAPFunctions(const SireMol::MoleculeInfoData &molinfo);

        CMAPFunctions(const CMAPFunctions &other);

        ~CMAPFunctions();

        static const char *typeName();
        const char *what() const;

        CMAPFunctions &operator=(const CMAPFunctions &other);

        bool operator==(const CMAPFunctions &other) const;
        bool operator!=(const CMAPFunctions &other) const;

        QString toString() const;

        virtual CMAPFunctions *clone() const;

        void set(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2, const AtomID &atom3,
                 const AtomID &atom4, const CMAPParameter &parameter);

        void set(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3,
                 AtomIdx atom4, const CMAPParameter &parameter);

        void clear();

        void clear(AtomIdx atom);
        void clear(const AtomID &atom);

        void clear(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3,
                   AtomIdx atom4);
        void clear(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2, const AtomID &atom3,
                   const AtomID &atom4);

        void clear(const QList<AtomIdx> &atoms, bool exclusive = true);

        bool isEmpty() const;

        int nFunctions() const;

        CMAPParameter parameter(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3,
                                AtomIdx atom4) const;

        CMAPParameter parameter(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2,
                                const AtomID &atom3, const AtomID &atom4) const;

        QVector<CMAPFunction> parameters() const;

        QVector<CMAPFunction> parameters(const QList<AtomIdx> &atoms, bool exclusive = true) const;

        CMAPFunctions includeOnly(const SireMol::AtomSelection &selected_atoms, bool isstrict = true) const;

        bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

        SireBase::PropertyList merge(const MolViewProperty &other,
                                     const SireMol::AtomIdxMapping &mapping,
                                     const QString &ghost = QString(),
                                     const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    protected:
        SireBase::PropertyPtr _pvt_makeCompatibleWith(const SireMol::MoleculeInfoData &molinfo,
                                                      const SireMol::AtomMatcher &atommatcher) const;

        SireBase::PropertyPtr _pvt_makeCompatibleWith(const SireMol::MoleculeInfoData &molinfo,
                                                      const QHash<AtomIdx, AtomIdx> &map) const;

        const SireMol::MoleculeInfoData &info() const;

    private:
        /** The metadata for the molecule that contains these functions */
        SireBase::SharedDataPointer<SireMol::MoleculeInfoData> molinfo;

        /** All of the CMAP parameters, identified by the IDQuint of
         *  that contains the parameters */
        QHash<detail::IDQuint, CMAPParameter> parameters_by_atoms;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

    //////
    ////// Inline functions of CMAPFunction
    //////

    /** Return the first atom of the quint */
    SIRE_ALWAYS_INLINE const CGAtomIdx &CMAPFunction::atom0() const
    {
        return atm0;
    }

    /** Return the second atom of the quint */
    SIRE_ALWAYS_INLINE const CGAtomIdx &CMAPFunction::atom1() const
    {
        return atm1;
    }

    /** Return the third atom of the quint */
    SIRE_ALWAYS_INLINE const CGAtomIdx &CMAPFunction::atom2() const
    {
        return atm2;
    }

    /** Return the fourth atom of the quint */
    SIRE_ALWAYS_INLINE const CGAtomIdx &CMAPFunction::atom3() const
    {
        return atm3;
    }

    /** Return the fifth atom of the quint */
    SIRE_ALWAYS_INLINE const CGAtomIdx &CMAPFunction::atom4() const
    {
        return atm4;
    }

    /** Return the parameter itself */
    SIRE_ALWAYS_INLINE const CMAPParameter &CMAPFunction::parameter() const
    {
        return param;
    }

#endif // SIRE_SKIP_INLINE_FUNCTIONS
}

Q_DECLARE_METATYPE(SireMM::CMAPFunction);
Q_DECLARE_METATYPE(SireMM::CMAPFunctions);

SIRE_EXPOSE_CLASS(SireMM::CMAPFunction)
SIRE_EXPOSE_CLASS(SireMM::CMAPFunctions)

SIRE_END_HEADER

#endif