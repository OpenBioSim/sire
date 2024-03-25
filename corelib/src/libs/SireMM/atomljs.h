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

#ifndef SIREMM_ATOMLJS_H
#define SIREMM_ATOMLJS_H

#include "SireMol/atomproperty.hpp"

#include "ljparameter.h"
#include "lj1264parameter.h"

#include <QHash>

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
    template <>
    class AtomProperty<SireMM::LJParameter>;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMol::AtomProperty<SireMM::LJParameter> &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMol::AtomProperty<SireMM::LJParameter> &);

namespace SireMM
{
    class LJExceptionID;
    class LJException;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::LJExceptionID &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::LJExceptionID &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::LJException &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::LJException &);

namespace SireMM
{
    /** This is a simple ID number that is used to match up LJ exceptions
     *  across different molecules. Each pair of atom-pair exceptions will
     *  have thier own ID number, and will recognise each other by this
     *  number. Note that this number is per-session, and will be changed
     *  when a molecule is streamed and re-loaded
     */
    class SIREMM_EXPORT LJExceptionID
    {
        friend QDataStream & ::operator<<(QDataStream &, const LJExceptionID &);
        friend QDataStream & ::operator>>(QDataStream &, LJExceptionID &);

        // friend so that is can see the 'is_first' flag
        friend class LJException;

    public:
        LJExceptionID();
        LJExceptionID(const LJExceptionID &other);
        ~LJExceptionID();

        static LJExceptionID generate();

        LJExceptionID &operator=(const LJExceptionID &other);

        bool operator==(const LJExceptionID &other) const;
        bool operator!=(const LJExceptionID &other) const;

        const char *what() const;
        static const char *typeName();

        bool pairsWith(const LJExceptionID &other) const;

        LJExceptionID getPair() const;

        QString toString() const;

    private:
        quint32 id;
        bool is_first;
    };

    /** This class represents a single LJException. It combines
     *  the LJExceptionID with the LJ1264Parameter value that
     *  is used to represent the exception
     */
    class SIREMM_EXPORT LJException
    {
        friend QDataStream & ::operator<<(QDataStream &, const LJException &);
        friend QDataStream & ::operator>>(QDataStream &, LJException &);

    public:
        LJException();
        LJException(const LJ1264Parameter &value);
        LJException(const LJExceptionID &id, const LJ1264Parameter &value);
        LJException(const LJException &other);
        ~LJException();

        LJException &operator=(const LJException &other);

        bool operator==(const LJException &other) const;
        bool operator!=(const LJException &other) const;

        const char *what() const;
        static const char *typeName();

        bool pairsWith(const LJException &other) const;

        LJExceptionID ID() const;
        LJ1264Parameter value() const;

        QString toString() const;

        LJException getPair() const;

    private:
        LJExceptionID id;
        LJ1264Parameter val;
    };
}

namespace SireMol
{
    /** This is an explicit specialisation of AtomProperty<T> for the LJParameter
        class, as the AtomLJs also has to store atom-pair exceptions

        @author Christopher Woods
    */
    template <>
    class SIREMM_EXPORT AtomProperty<SireMM::LJParameter> : public SireBase::ConcreteProperty<AtomProperty<SireMM::LJParameter>, AtomProp>
    {
        friend SIREMOL_EXPORT QDataStream & ::operator<<(QDataStream &, const AtomProperty<SireMM::LJParameter> &);
        friend SIREMOL_EXPORT QDataStream & ::operator>>(QDataStream &, AtomProperty<SireMM::LJParameter> &);

    public:
        typedef typename PackedArray2D<SireMM::LJParameter>::Array Array;

        AtomProperty();

        AtomProperty(const MoleculeInfo &molinfo);
        AtomProperty(const MoleculeInfo &molinfo, const SireMM::LJParameter &default_value);

        AtomProperty(const MoleculeView &molview);
        AtomProperty(const MoleculeView &molview, const SireMM::LJParameter &default_value);

        AtomProperty(const MoleculeInfoData &molinfo);
        AtomProperty(const MoleculeInfoData &molinfo, const SireMM::LJParameter &default_value);

        AtomProperty(const SireMM::LJParameter &value);
        AtomProperty(const PackedArray2D<SireMM::LJParameter> &values);

        AtomProperty(const AtomProperty<SireMM::LJParameter> &other);

        ~AtomProperty();

        AtomProperty<SireMM::LJParameter> &operator=(const AtomProperty<SireMM::LJParameter> &other);

        static const char *typeName();

        AtomProperty<SireMM::LJParameter> *clone() const;

        bool operator==(const AtomProperty<SireMM::LJParameter> &other) const;
        bool operator!=(const AtomProperty<SireMM::LJParameter> &other) const;

        bool isEmpty() const;

        QString toString() const;

        AtomProperty<QVariant> toVariant() const;
        static AtomProperty<SireMM::LJParameter> fromVariant(const AtomProperty<QVariant> &variant);

        void assignFrom(const AtomProperty<QVariant> &values);

        const typename PackedArray2D<SireMM::LJParameter>::Array &operator[](CGIdx cgidx) const;
        const typename PackedArray2D<SireMM::LJParameter>::Array &at(CGIdx cgidx) const;
        const typename PackedArray2D<SireMM::LJParameter>::Array &get(CGIdx cgidx) const;

        const SireMM::LJParameter &operator[](int i) const;
        const SireMM::LJParameter &operator[](const CGAtomIdx &cgatomidx) const;
        QList<SireMM::LJParameter> operator[](const QList<qint64> &idxs) const;
        QList<SireMM::LJParameter> operator[](const SireBase::Slice &slice) const;

        const SireMM::LJParameter &at(int i) const;
        const SireMM::LJParameter &at(const CGAtomIdx &cgatomidx) const;
        const SireMM::LJParameter &get(int i) const;
        const SireMM::LJParameter &get(const CGAtomIdx &cgatomidx) const;

        SireMM::LJ1264Parameter get(int i, int j) const;
        SireMM::LJ1264Parameter get(int i, int j, const AtomProperty<SireMM::LJParameter> &other) const;

        QVariant getAsVariant(const CGAtomIdx &cgatomidx) const;
        SireBase::PropertyPtr getAsProperty(const CGAtomIdx &cgatomidx) const;

        AtomProperty<SireMM::LJParameter> &set(int i, const SireMM::LJParameter &value);
        AtomProperty<SireMM::LJParameter> &set(const CGAtomIdx &cgatomidx, const SireMM::LJParameter &value);

        AtomProperty<SireMM::LJParameter> &set(CGIdx cgidx, const QVector<SireMM::LJParameter> &values);

        AtomProperty<SireMM::LJParameter> &set(int i, const QList<SireMM::LJException> &values);

        AtomProperty<SireMM::LJParameter> &set(int i, int j, const SireMM::LJ1264Parameter &value);
        AtomProperty<SireMM::LJParameter> &set(int i, int j,
                                               AtomProperty<SireMM::LJParameter> &other,
                                               const SireMM::LJ1264Parameter &value);

        void clearExceptions();
        void clearExceptions(int i);

        const PackedArray2D<SireMM::LJParameter> &array() const;

        const typename PackedArray2D<SireMM::LJParameter>::Array *data() const;
        const typename PackedArray2D<SireMM::LJParameter>::Array *constData() const;

        const SireMM::LJParameter *data(CGIdx cgidx) const;
        const SireMM::LJParameter *constData(CGIdx cgidx) const;

        int size() const;
        int count() const;

        int nCutGroups() const;

        int nAtoms() const;
        int nAtoms(CGIdx cgidx) const;

        bool isCompatibleWith(const MoleculeInfoData &molinfo) const;
        bool isCompatibleWith(const MoleculeInfo &molinfo) const;

        AtomProperty<SireMM::LJParameter> matchToSelection(const AtomSelection &selection) const;

        QVector<SireMM::LJParameter> toVector() const;
        QVector<SireMM::LJParameter> toVector(const AtomSelection &selection) const;

        QList<SireMM::LJParameter> toList() const;
        QList<SireMM::LJParameter> toList(const AtomSelection &selection) const;

        PropertyPtr merge(const MoleculeInfoData &molinfo) const;
        PropertyPtr divide(const QVector<AtomSelection> &beads) const;
        PropertyPtr divideByResidue(const MoleculeInfoData &molinfo) const;

        void copyFrom(const QVector<SireMM::LJParameter> &values);
        void copyFrom(const QVector<SireMM::LJParameter> &values, const AtomSelection &selection);

        bool canConvert(const QVariant &value) const;

        void assertCanConvert(const QVariant &value) const;

        bool hasExceptions() const;
        bool hasExceptions(int i) const;

        QList<SireMM::LJException> getExceptions(int i) const;

        QList<boost::tuple<int, int, SireMM::LJ1264Parameter>> getExceptions() const;
        QList<boost::tuple<int, int, SireMM::LJ1264Parameter>> getExceptions(const AtomProperty<SireMM::LJParameter> &other) const;

        SireBase::PropertyList merge(const MolViewProperty &other,
                                     const AtomIdxMapping &mapping,
                                     const QString &ghost = QString(),
                                     const SireBase::PropertyMap &map = SireBase::PropertyMap()) const;

    private:
        /** The actual atomic property values */
        SireBase::PackedArray2D<SireMM::LJParameter> props;

        /** The LJExceptionIDs of any atoms that have them - note that
         *  a single atom may have multiple exceptions, so this is a
         *  multi-hash
         */
        QHash<quint32, QList<SireMM::LJException>> lj_exceptions;
    };

}

namespace SireMM
{
    typedef SireMol::AtomProperty<LJParameter> AtomLJs;
}

Q_DECLARE_METATYPE(SireMM::AtomLJs);

Q_DECLARE_METATYPE(SireMM::LJException);
Q_DECLARE_METATYPE(SireMM::LJExceptionID);

Q_DECLARE_TYPEINFO(SireMM::LJException, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(SireMM::LJExceptionID, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS(SireMM::LJExceptionID)
SIRE_EXPOSE_CLASS(SireMM::LJException)

SIRE_EXPOSE_ATOM_PROPERTY(SireMM::LJParameter, SireMM::AtomLJs)

SIRE_END_HEADER

#endif
