/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "atomljs.h"

#include "SireBase/quickcopy.hpp"
#include "SireBase/incremint.h"

#include "SireStream/magic_error.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QUuid>
#include <QMutex>

using namespace SireMM;
using namespace SireMol;
using namespace SireStream;

///////
/////// Implementation of LJExceptionID
///////

static const RegisterMetaType<LJExceptionID> r_ljexceptionid(NO_ROOT);

// we need to transparently map unique IDs on streaming. We do this
// by using a UUID mapping on streaming and loading

static QMutex lj_exception_mutex;
static QHash<quint64, QUuid> uuid_to_id;
static QHash<QUuid, LJExceptionID> uuid_to_ljexceptionid;

QUuid get_uid(quint64 id)
{
    QMutexLocker lkr(&lj_exception_mutex);

    if (id_to_uuid.contains(id))
    {
        return id_to_uuid[id];
    }
    else
    {
        QUuid ret = QUuid::createUuid();
        id_to_uuid[id] = ret;
        return ret;
    }
}

LJExceptionID get_ljexceptionid(const QUuid &uuid)
{
    QMutexLocker lkr(&lj_exception_mutex);

    if (uuid_to_ljexceptionid.contains(uuid))
    {
        return uuid_to_ljexceptionid[uuid];
    }
    else
    {
        LJExceptionID ret = LJExceptionID::generate();
        uuid_to_ljexceptionid[uuid] = ret;
        return ret;
    }
}

QDataStream &operator<<(QDataStream &ds, const LJExceptionID &id)
{
    writeHeader(ds, r_ljexceptionid, 1);

    auto uid = get_uid(id.id);

    ds << uid;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, LJExceptionID &id)
{
    VersionID v = readHeader(ds, r_ljexceptionid);

    if (v == 1)
    {
        QUuid uid;
        ds >> uid;

        id = get_ljexceptionid(uid);
    }
    else
        throw SireStream::version_error(v, "1", r_ljexceptionid, CODELOC);

    return ds;
}

LJExceptionID::LJExceptionID() : id(0)
{
}

LJExceptionID::LJExceptionID(const LJExceptionID &other) : id(other.id)
{
}

LJExceptionID::~LJExceptionID()
{
}

static SireBase::Incremint global_ljexceptionid;

LJExceptionID LJExceptionID::generate()
{
    LJExceptionID ret;
    ret.id = global_ljexceptionid.increment();
    return ret;
}

LJExceptionID &LJExceptionID::operator=(const LJExceptionID &other)
{
    id = other.id;
    return *this;
}

bool LJExceptionID::operator==(const LJExceptionID &other) const
{
    return id == other.id;
}

bool LJExceptionID::operator!=(const LJExceptionID &other) const
{
    return id != other.id;
}

QString LJExceptionID::toString() const
{
    return QString("LJExceptionID(%1)").arg(id);
}

///////
/////// Implementation of LJParameter
///////

static const RegisterMetaType<AtomLJs> r_atomljs;

QDataStream &operator<<(QDataStream &ds, const AtomProperty<LJParameter> &atomljs)
{
    writeHeader(ds, r_atomljs, 2);

    SharedDataStream sds(ds);

    sds << static_cast<const SireMol::AtomProp &>(atomljs);
    sds << atomljs.props;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AtomProperty<LJParameter> &atomljs)
{
    VersionID v;

    try
    {
        v = readHeader(ds, r_atomljs);
    }
    catch (const SireStream::magic_error &)
    {
        // original version of AtomLJs didn't have a header!
        v = 1;
    }

    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> static_cast<SireMol::AtomProp &>(atomljs);
        sds >> atomljs.props;
    }
    else if (v == 1)
    {
        // only read the underlying AtomProperty
        ds >> static_cast<SireMol::AtomProp &>(atomljs);
        ds >> atomljs.props;
    }

    return ds;
}

/** Return whether or not the variant 'value' can be converted to be
    held as an AtomProperty<LJParameter> */
bool AtomProperty<LJParameter>::canConvert(const QVariant &value) const
{
    return value.isNull() or value.canConvert<LJParameter>();
}

/** Assert that the passed variant value can be converted to be
    held within this property

    \throw SireError::invalid_cast
*/
void AtomProperty<LJParameter>::assertCanConvert(const QVariant &value) const
{
    if (not(value.isNull() or value.canConvert<LJParameter>()))
    {
        throw SireError::invalid_cast(QObject::tr("It is not possible to convert the value of type %1 to "
                                                  "type %2, as is required for storing in the AtomProperty %3.")
                                          .arg(value.typeName())
                                          .arg(QMetaType::typeName(qMetaTypeId<LJParameter>()))
                                          .arg(this->what()),
                                      CODELOC);
    }
}

/** Null constructor */
AtomProperty<LJParameter>::AtomProperty() : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
}

/** Create an AtomProperty that holds one value for each
    atom described in 'molinfo'. Each atom starts with
    a default-constructed value of the property */
AtomProperty<LJParameter>::AtomProperty(const MoleculeInfoData &molinfo)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
    int ncg = molinfo.nCutGroups();

    if (ncg > 0)
    {
        // create space for each CutGroup
        QVector<QVector<LJParameter>> tmp_props = QVector<QVector<LJParameter>>(ncg);
        QVector<LJParameter> *tmp_props_array = tmp_props.data();

        for (CGIdx i(0); i < ncg; ++i)
        {
            // now create space for all of the atoms
            tmp_props_array[i] = QVector<LJParameter>(molinfo.nAtoms(i));
        }

        // now copy this into the PackedArray
        props = PackedArray2D<LJParameter>(tmp_props);
    }
}

/** Create an AtomProperty that holds one value for each
    atom described in 'molinfo'. Each atom starts with
    a default-constructed value of the property */
AtomProperty<LJParameter>::AtomProperty(const MoleculeInfo &molinfo)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
    this->operator=(AtomProperty<LJParameter>(molinfo.data()));
}

AtomProperty<LJParameter>::AtomProperty(const MoleculeInfo &molinfo, const LJParameter &default_value)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
    this->operator=(AtomProperty<LJParameter>(molinfo.data(), default_value));
}

AtomProperty<LJParameter>::AtomProperty(const MoleculeView &molview)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
    this->operator=(AtomProperty<LJParameter>(MoleculeInfo(molview)));
}

AtomProperty<LJParameter>::AtomProperty(const MoleculeView &molview, const LJParameter &default_value)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
    this->operator=(AtomProperty<LJParameter>(MoleculeInfo(molview), default_value));
}

/** Create an AtomProperty that holds one value for each
    atom described in 'molinfo'. Each atom starts with
    the value 'default_value' */
AtomProperty<LJParameter>::AtomProperty(const MoleculeInfoData &molinfo, const LJParameter &default_value)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
    int ncg = molinfo.nCutGroups();

    if (ncg > 0)
    {
        // create space for each CutGroup
        QVector<QVector<LJParameter>> tmp_props = QVector<QVector<LJParameter>>(ncg);
        QVector<LJParameter> *tmp_props_array = tmp_props.data();

        for (CGIdx i(0); i < ncg; ++i)
        {
            // now create space for all of the atoms
            tmp_props_array[i] = QVector<LJParameter>(molinfo.nAtoms(i), default_value);
        }

        // now copy this into the PackedArray
        props = PackedArray2D<LJParameter>(tmp_props);
    }
}

/** Construct an Atom property that holds a single value (only
    suitable for a molecule that has just one atom in just one CutGroup) */
AtomProperty<LJParameter>::AtomProperty(const LJParameter &value)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>()
{
    QVector<LJParameter> tmp_props(1, value);
    props = PackedArray2D<LJParameter>(tmp_props);
}

/** Construct the Atom property from the PackedArray2D of values */
AtomProperty<LJParameter>::AtomProperty(const PackedArray2D<LJParameter> &values)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>(), props(values)
{
}

/** Copy constructor */
AtomProperty<LJParameter>::AtomProperty(const AtomProperty<LJParameter> &other)
    : SireBase::ConcreteProperty<AtomProperty<LJParameter>, AtomProp>(), props(other.props)
{
}

/** Destructor */
AtomProperty<LJParameter>::~AtomProperty()
{
}

/** Copy assignment operator */
AtomProperty<LJParameter> &AtomProperty<LJParameter>::operator=(const AtomProperty<LJParameter> &other)
{
    props = other.props;
    return *this;
}

/** Comparison operator */
bool AtomProperty<LJParameter>::operator==(const AtomProperty<LJParameter> &other) const
{
    return props == other.props;
}

/** Comparison operator */
bool AtomProperty<LJParameter>::operator!=(const AtomProperty<LJParameter> &other) const
{
    return props != other.props;
}

AtomProperty<LJParameter> *AtomProperty<LJParameter>::clone() const
{
    return new AtomProperty<LJParameter>(*this);
}

const char *AtomProperty<LJParameter>::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AtomProperty<LJParameter>>());
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'

    \throw SireError::invalid_index
*/
const typename PackedArray2D<LJParameter>::Array &AtomProperty<LJParameter>::operator[](CGIdx cgidx) const
{
    return props.constData()[cgidx.map(props.count())];
}

/** Convert the contained properties into an array of arrays of QVariants.
    There is one array per CutGroup */
AtomProperty<QVariant> AtomProperty<LJParameter>::toVariant() const
{
    return AtomProperty<QVariant>(props.toVariant());
}

/** Assign the values of the properties from the array of QVariants

    \throw SireError::invalid_cast
*/
void AtomProperty<LJParameter>::assignFrom(const AtomProperty<QVariant> &variant)
{
    props = SireBase::PackedArray2D<LJParameter>::fromVariant(variant.array());
}

/** Return an AtomProperty constructed from an array of QVariants */
AtomProperty<LJParameter> AtomProperty<LJParameter>::fromVariant(const AtomProperty<QVariant> &variant)
{
    return AtomProperty<LJParameter>(SireBase::PackedArray2D<LJParameter>::fromVariant(variant.array()));
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'

    \throw SireError::invalid_index
*/
const typename PackedArray2D<LJParameter>::Array &AtomProperty<LJParameter>::at(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the array of properties for the atoms in the CutGroup
    identified by index 'cgidx'

    \throw SireError::invalid_index
*/
const typename PackedArray2D<LJParameter>::Array &AtomProperty<LJParameter>::get(CGIdx cgidx) const
{
    return this->operator[](cgidx);
}

/** Return the property at the specified index */
const LJParameter &AtomProperty<LJParameter>::operator[](int i) const
{
    i = SireID::Index(i).map(this->nAtoms());

    return props.constValueData()[i];
}

const LJParameter &AtomProperty<LJParameter>::at(int i) const
{
    return this->operator[](i);
}

const LJParameter &AtomProperty<LJParameter>::get(int i) const
{
    return this->operator[](i);
}

QList<LJParameter> AtomProperty<LJParameter>::operator[](const QList<qint64> &idxs) const
{
    QList<LJParameter> ret;

    for (auto idx : idxs)
    {
        ret.append(this->operator[](idx));
    }

    return ret;
}

QList<LJParameter> AtomProperty<LJParameter>::operator[](const SireBase::Slice &slice) const
{
    QList<LJParameter> ret;

    for (auto it = slice.begin(this->nAtoms()); not it.atEnd(); it.next())
    {
        ret.append(this->operator[](it.value()));
    }

    return ret;
}

/** Return the property for the atom at index 'cgatomidx'

    \throw SireError::invalid_index
*/
const LJParameter &AtomProperty<LJParameter>::operator[](const CGAtomIdx &cgatomidx) const
{
    const typename PackedArray2D<LJParameter>::Array &group_props = this->operator[](cgatomidx.cutGroup());

    return group_props.constData()[cgatomidx.atom().map(group_props.count())];
}

/** Return the property for the atom at index 'cgatomidx'

    \throw SireError::invalid_index
*/
const LJParameter &AtomProperty<LJParameter>::at(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Return the property for the atom at index 'cgatomidx'

    \throw SireError::invalid_index
*/
const LJParameter &AtomProperty<LJParameter>::get(const CGAtomIdx &cgatomidx) const
{
    return this->operator[](cgatomidx);
}

/** Return the value for the atom at index 'cgatomidx', as
    a QVariant. This lets you get the value without knowing the
    actual type of this AtomProperty<LJParameter>

    \throw SireError::invalid_index
*/
QVariant AtomProperty<LJParameter>::getAsVariant(const CGAtomIdx &cgatomidx) const
{
    const auto &value = this->get(cgatomidx);

    return QVariant::fromValue(value);
}

/** Return the value for the atom at index 'cgatomidx' as a
    Property. This lets you get the value without knowing the
    actual type of this AtomProperty<LJParameter>

   \throw SireError::invalid_index
*/
SireBase::PropertyPtr AtomProperty<LJParameter>::getAsProperty(const CGAtomIdx &cgatomidx) const
{
    return SireBase::convert_property(this->get(cgatomidx));
}

/** Set the value of the property of the ith atoms to 'value'
 *
 * \throw SireError::invalid_index
 */
AtomProperty<LJParameter> &AtomProperty<LJParameter>::set(int i, const LJParameter &value)
{
    i = SireID::Index(i).map(this->nAtoms());

    props.valueData()[i] = value;

    return *this;
}

/** Set the value of the property for the atom at index 'cgatomidx'

    \throw SireError::invalid_index
*/
AtomProperty<LJParameter> &AtomProperty<LJParameter>::set(const CGAtomIdx &cgatomidx, const LJParameter &value)
{
    quint32 cgidx = cgatomidx.cutGroup().map(props.count());
    quint32 atomidx = cgatomidx.atom().map(props.at(cgidx).count());

    props(cgidx, atomidx) = value;

    return *this;
}

/** Set the values for all atoms in the CutGroup at index 'cgidx'

    \throw SireError::incompatible_error
    \throw SireError::invalid_index
*/
AtomProperty<LJParameter> &AtomProperty<LJParameter>::set(CGIdx cgidx, const QVector<LJParameter> &values)
{
    props.update(cgidx.map(props.count()), values);

    return *this;
}

/** Return a const-reference to the PackedArray2D used to store
    all of the atom properties */
const PackedArray2D<LJParameter> &AtomProperty<LJParameter>::array() const
{
    return props;
}

/** Return a raw pointer to the array of arrays */
const typename PackedArray2D<LJParameter>::Array *AtomProperty<LJParameter>::data() const
{
    return props.constData();
}

/** Return a raw pointer to the array of arrays */
const typename PackedArray2D<LJParameter>::Array *AtomProperty<LJParameter>::constData() const
{
    return props.constData();
}

/** Return a raw pointer to the array of properties for
    the atoms in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
const LJParameter *AtomProperty<LJParameter>::data(CGIdx cgidx) const
{
    return this->at(cgidx).constData();
}

/** Return a raw pointer to the array of properties for
    the atoms in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
const LJParameter *AtomProperty<LJParameter>::constData(CGIdx cgidx) const
{
    return this->at(cgidx).constData();
}

/** Return the number of atoms in the molecule */
int AtomProperty<LJParameter>::size() const
{
    return this->nAtoms();
}

/** Return the number of atoms in the molecule */
int AtomProperty<LJParameter>::count() const
{
    return this->nAtoms();
}

/** Return the number of CutGroups in the molecule */
int AtomProperty<LJParameter>::nCutGroups() const
{
    return props.count();
}

/** Return whether or not this is empty */
bool AtomProperty<LJParameter>::isEmpty() const
{
    return this->nAtoms() == 0;
}

QString AtomProperty<LJParameter>::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("AtomLJs::empty");
    }
    else
    {
        QStringList parts;

        const auto n = this->count();

        if (n <= 10)
        {
            for (int i = 0; i < n; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->operator[](i).toString()));
            }
        }
        else
        {
            for (int i = 0; i < 5; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->operator[](i).toString()));
            }

            parts.append("...");

            for (int i = n - 5; i < n; ++i)
            {
                parts.append(QObject::tr("%1: %2").arg(i).arg(this->operator[](i).toString()));
            }
        }

        return QObject::tr("AtomLJs( size=%2\n%3\n)").arg(n).arg(parts.join("\n"));
    }
}

/** Return the total number of atoms in the molecule */
int AtomProperty<LJParameter>::nAtoms() const
{
    return props.nValues();
}

/** Return the number of atoms in the CutGroup at index 'cgidx'

    \throw SireError::invalid_index
*/
int AtomProperty<LJParameter>::nAtoms(CGIdx cgidx) const
{
    return this->at(cgidx).count();
}

/** Return whether or not this property is compatible with the
    Molecule whose layout is described in 'molinfo'
*/
bool AtomProperty<LJParameter>::isCompatibleWith(const MoleculeInfoData &molinfo) const
{
    if (props.nValues() != molinfo.nAtoms())
        return false;

    int ncg = molinfo.nCutGroups();

    if (ncg != props.count())
        return false;

    const typename PackedArray2D<LJParameter>::Array *props_array = props.constData();

    for (CGIdx i(0); i < ncg; ++i)
    {
        if (molinfo.nAtoms(i) != props_array[i].count())
            return false;
    }

    return true;
}

/** Return whether or not this property is compatible with the
    Molecule whose layout is described in 'molinfo'
*/
bool AtomProperty<LJParameter>::isCompatibleWith(const MoleculeInfo &molinfo) const
{
    if (props.nValues() != molinfo.nAtoms())
        return false;

    int ncg = molinfo.nCutGroups();

    if (ncg != props.count())
        return false;

    const typename PackedArray2D<LJParameter>::Array *props_array = props.constData();

    for (CGIdx i(0); i < ncg; ++i)
    {
        if (molinfo.nAtoms(i) != props_array[i].count())
            return false;
    }

    return true;
}

/** Convert this atom property to an array of values. The values
    are written in CGAtomIdx order */
QVector<LJParameter> AtomProperty<LJParameter>::toVector() const
{
    if (this->nAtoms() == 0)
        return QVector<LJParameter>();

    QVector<LJParameter> ret(this->nAtoms());

    SireBase::quickCopy<LJParameter>(ret.data(), props.constValueData(), this->nAtoms());

    return ret;
}

QList<LJParameter> AtomProperty<LJParameter>::toList() const
{
    const int nats = this->nAtoms();

    if (nats == 0)
        return QList<LJParameter>();

    QList<LJParameter> ret;
    ret.reserve(nats);

    for (int i = 0; i < nats; ++i)
    {
        ret.append(props.constValueData()[i]);
    }

    return ret;
}

/** Convert the properties of the atoms selected in 'selection' to an
    array of values. The values are written in CGAtomIdx order

    \throw SireError::incompatible_error
*/
QVector<LJParameter> AtomProperty<LJParameter>::toVector(const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*this);

    if (selected_atoms.selectedAll())
        return this->toVector();

    else if (selected_atoms.selectedNone())
        return QVector<LJParameter>();

    else if (selected_atoms.selectedAllCutGroups())
    {
        QVector<LJParameter> vals(selected_atoms.nSelected());
        LJParameter *value = vals.data();

        const int ncg = selected_atoms.nCutGroups();
        const typename PackedArray2D<LJParameter>::Array *props_array = props.constData();

        for (CGIdx i(0); i < ncg; ++i)
        {
            const LJParameter *group_props = props_array[i].constData();

            if (selected_atoms.selectedAll(i))
            {
                const int nats = props_array[i].nValues();

                SireBase::quickCopy<LJParameter>(value, group_props, nats);

                value += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (SireID::Index idx, idxs)
                {
                    *value = group_props[idx];
                    ++value;
                }
            }
        }

        return vals;
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        std::sort(cgidxs.begin(), cgidxs.end());

        QVector<LJParameter> vals(selected_atoms.nSelected());
        LJParameter *value = vals.data();

        const typename PackedArray2D<LJParameter>::Array *props_array = props.constData();

        foreach (CGIdx i, cgidxs)
        {
            const LJParameter *group_props = props_array[i].constData();

            if (selected_atoms.selectedAll(i))
            {
                const int nats = props_array[i].nValues();

                SireBase::quickCopy<LJParameter>(value, group_props, nats);

                value += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (SireID::Index idx, idxs)
                {
                    *value = group_props[idx];
                    ++value;
                }
            }
        }

        return vals;
    }
}

QList<LJParameter> AtomProperty<LJParameter>::toList(const AtomSelection &selected_atoms) const
{
    return this->toVector().toList();
}

/** Copy into this atom property set the values from 'values'. The values
    are copied in CGAtomIdx order, and there must be as many values
    as there are atoms

    \throw SireError::incompatible_error
*/
void AtomProperty<LJParameter>::copyFrom(const QVector<LJParameter> &values)
{
    if (values.count() != this->nAtoms())
        this->throwIncorrectNumberOfAtoms(values.count(), this->nAtoms());

    SireBase::quickCopy<LJParameter>(props.valueData(), values.constData(), values.count());
}

/** Copy into this atom property set the values from 'values', but only
    for the atoms selected in 'selection'. This copies the properties
    in in CGAtomIdx order, and there must be the same number of values
    as there are selected atoms

    \throw SireError::incompatible_error
*/
void AtomProperty<LJParameter>::copyFrom(const QVector<LJParameter> &values, const AtomSelection &selected_atoms)
{
    selected_atoms.assertCompatibleWith(*this);

    if (selected_atoms.selectedAll())
    {
        this->copyFrom(values);
        return;
    }

    if (values.count() != selected_atoms.nSelected())
        this->throwIncorrectNumberOfSelectedAtoms(values.count(), selected_atoms.nSelected());

    const LJParameter *values_array = values.constData();

    if (selected_atoms.selectedAllCutGroups())
    {
        const int ncg = selected_atoms.nCutGroups();

        for (CGIdx i(0); i < ncg; ++i)
        {
            LJParameter *group_props = props.data(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = props.nValues(i);
                SireBase::quickCopy<LJParameter>(group_props, values_array, nats);

                values_array += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (SireID::Index idx, idxs)
                {
                    group_props[idx] = *values_array;
                    ++values_array;
                }
            }
        }
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();
        std::sort(cgidxs.begin(), cgidxs.end());

        foreach (CGIdx i, cgidxs)
        {
            LJParameter *group_props = props.data(i);

            if (selected_atoms.selectedAll(i))
            {
                const int nats = props.nValues(i);
                SireBase::quickCopy<LJParameter>(group_props, values_array, nats);

                values_array += nats;
            }
            else
            {
                QList<SireID::Index> idxs = selected_atoms.selectedAtoms(i).values();
                std::sort(idxs.begin(), idxs.end());

                foreach (SireID::Index idx, idxs)
                {
                    group_props[idx] = *values_array;
                    ++values_array;
                }
            }
        }
    }
}

/** Match this property to the passed selection. This returns
    the property only for the CutGroups that have been selected,
    and with default values for any atoms in those CutGroups that
    have not been selected. This is useful, e.g. for the forcefield
    classes, as this allows an AtomProperty<LJParameter> to be returned
    for only the atoms that are selected as part of the forcefield.

    \throw SireError::incompatible_error
*/
AtomProperty<LJParameter> AtomProperty<LJParameter>::matchToSelection(const AtomSelection &selected_atoms) const
{
    selected_atoms.assertCompatibleWith(*this);

    if (selected_atoms.selectedAll())
        return *this;
    else if (selected_atoms.selectedNone())
        return AtomProperty<LJParameter>();

    else if (selected_atoms.selectedAllCutGroups())
    {
        PackedArray2D<LJParameter> new_props = props;

        int ncg = props.count();

        for (CGIdx i(0); i < ncg; ++i)
        {
            if (not selected_atoms.selectedAll(i))
            {
                int nats = new_props.at(i).count();
                LJParameter *atom_props_array = new_props.data(i);

                const QSet<Index> &atoms = selected_atoms.selectedAtoms(i);

                for (Index j(0); j < nats; ++j)
                {
                    if (not atoms.contains(j))
                        atom_props_array[j] = LJParameter::dummy();
                }
            }
        }

        return AtomProperty<LJParameter>(new_props);
    }
    else
    {
        QList<CGIdx> cgidxs = selected_atoms.selectedCutGroups();

        QVector<QVector<LJParameter>> new_props = QVector<QVector<LJParameter>>(cgidxs.count());
        QVector<LJParameter> *new_props_array = new_props.data();

        const typename PackedArray2D<LJParameter>::Array *props_array = props.constData();

        int n = 0;

        foreach (CGIdx i, cgidxs)
        {
            if (selected_atoms.selectedAll(i))
            {
                new_props_array[n] = props_array[i].toQVector();
                ++n;
            }
            else
            {
                const QSet<Index> &atoms = selected_atoms.selectedAtoms(i);

                QVector<LJParameter> atom_props = props_array[i].toQVector();
                int nats = atom_props.count();

                LJParameter *atom_props_array = atom_props.data();

                for (Index j(0); j < nats; ++j)
                {
                    if (not atoms.contains(j))
                        atom_props_array[j] = LJParameter::dummy();
                }

                new_props_array[n] = atom_props;
                ++n;
            }
        }

        return AtomProperty<LJParameter>(new_props);
    }
}

/** Merge all of the atomic properties into a single array, with
    the properties arranged in AtomIdx order */
PropertyPtr AtomProperty<LJParameter>::merge(const MoleculeInfoData &moldata) const
{
    this->assertCompatibleWith(moldata);

    QVector<LJParameter> vals(moldata.nAtoms());

    LJParameter *vals_array = vals.data();

    for (AtomIdx i(0); i < moldata.nAtoms(); ++i)
    {
        vals_array[i] = this->at(moldata.cgAtomIdx(i));
    }

    return AtomProperty<LJParameter>(vals);
}

/** Divide the AtomProperty into beads according to the passed atom selections,
    and returning the properties in AtomIdx order within each bead

    \throw SireError::incompatible_error
*/
PropertyPtr AtomProperty<LJParameter>::divide(const QVector<AtomSelection> &beads) const
{
    if (beads.isEmpty())
        return PropertyPtr();

    const int nbeads = beads.count();
    const AtomSelection *beads_array = beads.constData();

    QVector<QVector<LJParameter>> bead_vals(nbeads);
    QVector<LJParameter> *bead_vals_array = bead_vals.data();

    for (int i = 0; i < nbeads; ++i)
    {
        const AtomSelection &bead = beads_array[i];

        bead.assertCompatibleWith<LJParameter>(*this);

        QVector<LJParameter> vals(bead.nSelected());
        LJParameter *vals_array = vals.data();

        if (bead.selectedAll())
        {
            for (AtomIdx j(0); j < bead.nSelected(); ++j)
            {
                vals_array[j] = this->at(bead.info().cgAtomIdx(j));
            }
        }
        else
        {
            foreach (const AtomIdx &j, bead.selectedAtoms())
            {
                *vals_array = this->at(bead.info().cgAtomIdx(j));
                ++vals_array;
            }
        }

        bead_vals_array[i] = vals;
    }

    return AtomProperty<LJParameter>(bead_vals);
}

/** Divide the properties into residues. This returns the values in
    Residue/Index order

    \throw SireError::incompatible_error
*/
PropertyPtr AtomProperty<LJParameter>::divideByResidue(const MoleculeInfoData &molinfo) const
{
    this->assertCompatibleWith(molinfo);

    QVector<QVector<LJParameter>> res_vals(molinfo.nResidues());
    QVector<LJParameter> *res_vals_array = res_vals.data();

    for (ResIdx i(0); i < molinfo.nResidues(); ++i)
    {
        const int nats = molinfo.nAtoms(i);

        QVector<LJParameter> vals(nats);
        LJParameter *vals_array = vals.data();

        for (int j = 0; j < nats; ++j)
        {
            vals_array[j] = this->at(molinfo.cgAtomIdx(molinfo.getAtom(i, j)));
        }

        res_vals_array[i] = vals;
    }

    return AtomProperty<LJParameter>(res_vals);
}

/** Return the combined LJ parameter for atoms i and j, assuming either
 *  standard combining rules (arithmetic) or returning the exception
 *  if one has been set
 */
LJ1264Parameter AtomProperty<LJParameter>::get(int i, int j) const
{
    return LJ1264Parameter(this->get(i).combine(this->get(j)));
}

/** Return the LJ exception between atom 'i' and the exception representing
 *  atom 'j'. If no exception has been set, then the standard LJ parameter
 *  for atom 'i' is returned
 */
LJ1264Parameter AtomProperty<LJParameter>::get(int i, LJExceptionID j) const
{
    return LJ1264Parameter(this->get(i));
}

/** Return the LJ exception between the exception representing atom 'i' and
 *  atom 'j'. If no exception has been set, then the standard LJ parameter
 * for atom 'j' is returned
 */
LJ1264Parameter AtomProperty<LJParameter>::get(LJExceptionID i, int j) const
{
    return LJ1264Parameter(this->get(j));
}

/** Return the LJ exception between atom 'i' in this set and atom 'j'
 *  in other. If no exception has been set, then the standard LJ parameter
 *  combination using arithmetic combining rules is returned.
 */
LJ1264Parameter AtomProperty<LJParameter>::get(int i, int j, const AtomProperty<LJParameter> &other) const
{
    return LJ1264Parameter(this->get(i).combine(other.get(j)));
}

/** Set the exception for atoms i and j equal to 'value' */
AtomProperty<LJParameter> &AtomProperty<LJParameter>::set(int i, int j, const LJ1264Parameter &value)
{
    return *this;
}

/** Set the exception for atom i and the exception representing 'j' to 'value' */
AtomProperty<LJParameter> &AtomProperty<LJParameter>::set(int i, LJExceptionID j, const LJ1264Parameter &value)
{
    return *this;
}

/** Set the exception for the exception representing 'i' and atom j to 'value' */
AtomProperty<LJParameter> &AtomProperty<LJParameter>::set(LJExceptionID i, int j, const LJ1264Parameter &value)
{
    return *this;
}

/** Return whether or not there are any exceptions in this AtomProperty */
bool AtomProperty<LJParameter>::hasExceptions() const
{
    return false;
}

/** Return all of the exceptions between this LJ set and 'other' */
QList<std::tuple<int, int, SireMM::LJ1264Parameter>> AtomProperty<LJParameter>::getExceptions(const AtomProperty<SireMM::LJParameter> &other) const
{
    return QList<std::tuple<int, int, SireMM::LJ1264Parameter>>();
}
