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

#include "energytrajectory.h"

#include "SireBase/console.h"

#include "SireUnits/units.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;

static const RegisterMetaType<EnergyTrajectory> r_etraj;

QDataStream &operator<<(QDataStream &ds, const EnergyTrajectory &etraj)
{
    writeHeader(ds, r_etraj, 1);

    SharedDataStream sds(ds);

    sds << etraj.time_values << etraj.energy_values
        << etraj.label_values << etraj.props
        << static_cast<const Property &>(etraj);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, EnergyTrajectory &etraj)
{
    VersionID v = readHeader(ds, r_etraj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> etraj.time_values >> etraj.energy_values >> etraj.label_values >> etraj.props >> static_cast<Property &>(etraj);
    }
    else
        throw version_error(v, "1", r_etraj, CODELOC);

    return ds;
}

EnergyTrajectory::EnergyTrajectory()
    : ConcreteProperty<EnergyTrajectory, Property>()
{
}

EnergyTrajectory::EnergyTrajectory(const EnergyTrajectory &other)
    : ConcreteProperty<EnergyTrajectory, Property>(other),
      time_values(other.time_values), energy_values(other.energy_values),
      label_values(other.label_values), props(other.props)
{
}

EnergyTrajectory::~EnergyTrajectory()
{
}

EnergyTrajectory &EnergyTrajectory::operator=(const EnergyTrajectory &other)
{
    if (this != &other)
    {
        time_values = other.time_values;
        energy_values = other.energy_values;
        label_values = other.label_values;
        props = other.props;
        Property::operator=(other);
    }

    return *this;
}

bool EnergyTrajectory::operator==(const EnergyTrajectory &other) const
{
    return time_values == other.time_values and
           energy_values == other.energy_values and
           label_values == other.label_values and
           props == other.props;
}

bool EnergyTrajectory::operator!=(const EnergyTrajectory &other) const
{
    return not this->operator==(other);
}

EnergyTrajectory *EnergyTrajectory::clone() const
{
    return new EnergyTrajectory(*this);
}

const char *EnergyTrajectory::typeName()
{
    return QMetaType::typeName(qMetaTypeId<EnergyTrajectory>());
}

const char *EnergyTrajectory::what() const
{
    return EnergyTrajectory::typeName();
}

QString EnergyTrajectory::toString() const
{
    if (this->isEmpty())
        return QObject::tr("EnergyTrajectory::empty");

    QStringList parts;

    const auto n = this->count();

    auto keys = this->energy_values.keys();
    keys.sort();

    auto label_keys = this->label_values.keys();
    label_keys.sort();

    auto headers = label_keys + keys;

    headers.insert(0, "time");

    parts.append(headers.join("\t"));

    if (n <= 10)
    {
        for (int i = 0; i < n; ++i)
        {
            const auto vals = this->get(i, GeneralUnit(picosecond), GeneralUnit(kcal_per_mol));
            const auto label_vals = this->getLabels(i);

            QStringList v;

            v.append(QString::number(vals["time"]));

            for (const auto &key : label_keys)
            {
                v.append(label_vals[key]);
            }

            for (const auto &key : keys)
            {
                v.append(QString::number(vals[key]));
            }

            parts.append(v.join("\t"));
        }
    }
    else
    {
        for (int i = 0; i < 5; ++i)
        {
            const auto vals = this->get(i, GeneralUnit(picosecond), GeneralUnit(kcal_per_mol));
            const auto label_vals = this->getLabels(i);

            QStringList v;

            v.append(QString::number(vals["time"]));

            for (const auto &key : label_keys)
            {
                v.append(label_vals[key]);
            }

            for (const auto &key : keys)
            {
                v.append(QString::number(vals[key]));
            }

            parts.append(v.join("\t"));
        }

        parts.append("...");

        for (int i = n - 5; i < n; ++i)
        {
            const auto vals = this->get(i, GeneralUnit(picosecond), GeneralUnit(kcal_per_mol));
            const auto label_vals = this->getLabels(i);

            QStringList v;

            v.append(QString::number(vals["time"]));

            for (const auto &key : label_keys)
            {
                v.append(label_vals[key]);
            }

            for (const auto &key : keys)
            {
                v.append(QString::number(vals[key]));
            }

            parts.append(v.join("\t"));
        }
    }

    return QObject::tr("EnergyTrajectory( size=%1\n%2\n)").arg(n).arg(parts.join("\n"));
}

/** Set the energies at time 'time' to the components contained
 *  in 'energies', and the labels to those contained in 'labels'
 */
void EnergyTrajectory::set(const GeneralUnit &time,
                           const QHash<QString, GeneralUnit> &energies,
                           const QHash<QString, QString> &labels)
{
    if (energies.isEmpty() and labels.isEmpty())
        return;

    auto t = Time(time);

    // round the time to 5 decimal places
    auto t_round = qRound(t.value() * 100000.0) / 100000.0;

    for (auto it = energies.constBegin(); it != energies.constEnd(); ++it)
    {
        // make sure all of the values are valid
        if (it.key() == "time")
        {
            throw SireError::invalid_key(QObject::tr(
                                             "You cannot call an energy component 'time'. This name "
                                             "is reserved for the time component."),
                                         CODELOC);
        }

        auto nrg = MolarEnergy(it.value());
    }

    for (auto it = labels.constBegin(); it != labels.constEnd(); ++it)
    {
        // make sure all of the values are valid
        if (it.key() == "time")
        {
            throw SireError::invalid_key(QObject::tr(
                                             "You cannot call a label component 'time'. This name "
                                             "is reserved for the time component."),
                                         CODELOC);
        }
    }

    // we won't check for duplication in the energies and labels as these
    // can be resolved on printout

    // make sure we have populated all of the columns - this
    // populates with NaNs if the column doesn't exist
    QVector<double> null_column;

    for (const auto &key : energies.keys())
    {
        if (not this->energy_values.contains(key))
        {
            if (null_column.isEmpty() and not time_values.isEmpty())
            {
                null_column = QVector<double>(time_values.count(), NAN);
            }

            this->energy_values.insert(key, null_column);
        }
    }

    QVector<QString> null_label_column;

    for (const auto &key : labels.keys())
    {
        if (not this->label_values.contains(key))
        {
            if (null_label_column.isEmpty() and not time_values.isEmpty())
            {
                null_label_column = QVector<QString>(time_values.count());
            }

            this->label_values.insert(key, null_label_column);
        }
    }

    // find the index of the time, and add it as needed.
    // Times are sorted, so start from the end and work towards
    // the beginning (as we will normally be adding energies onto
    // the end)
    int idx = time_values.count();
    bool must_create = true;

    while (idx > 0)
    {
        if (t_round > time_values[idx - 1])
            break;

        else if (t_round == time_values[idx - 1])
        {
            SireBase::Console::warning(QObject::tr(
                "EnergyTrajectory::set: time %1 already exists in the trajectory. "
                "Overwriting existing values.").arg(t.toString()));
            must_create = false;
            idx = idx - 1;
            break;
        }

        idx -= 1;
    }

    if (idx >= time_values.count())
    {
        time_values.insert(idx, t_round);

        for (auto it = this->energy_values.begin();
             it != this->energy_values.end(); ++it)
        {
            it.value().append(NAN);
        }

        for (auto it = this->label_values.begin();
             it != this->label_values.end(); ++it)
        {
            it.value().append(QString());
        }
    }
    else if (must_create)
    {
        time_values.insert(idx, t_round);

        for (auto it = this->energy_values.begin();
             it != this->energy_values.end(); ++it)
        {
            it.value().insert(idx, NAN);
        }

        for (auto it = this->label_values.begin();
             it != this->label_values.end(); ++it)
        {
            it.value().insert(idx, QString());
        }
    }

    for (auto it = energies.constBegin();
         it != energies.constEnd(); ++it)
    {
        this->energy_values[it.key()][idx] = MolarEnergy(it.value()).value();
    }

    for (auto it = labels.constBegin();
         it != labels.constEnd(); ++it)
    {
        this->label_values[it.key()][idx] = it.value();
    }
}

/** Set the energies at time 'time' to the components contained
 *  in 'energies'
 */
void EnergyTrajectory::set(const GeneralUnit &time,
                           const QHash<QString, GeneralUnit> &energies)
{
    this->set(time, energies, QHash<QString, QString>());
}

/** Return whether or not this is empty */
bool EnergyTrajectory::isEmpty() const
{
    return this->time_values.isEmpty();
}

bool EnergyTrajectory::isNull() const
{
    return this->isEmpty();
}

/** Return the number of time values (number of rows) */
int EnergyTrajectory::count() const
{
    return this->time_values.count();
}

/** Return the number of time values (number of rows) */
int EnergyTrajectory::size() const
{
    return this->count();
}

/** Return the time and energy components at the ith row.
 *  Values are returned in internal units
 */
QHash<QString, double> EnergyTrajectory::operator[](int i) const
{
    i = SireID::Index(i).map(this->count());

    QHash<QString, double> vals;
    vals.reserve(1 + this->energy_values.count());

    vals.insert("time", this->time_values[i]);

    for (auto it = this->energy_values.constBegin();
         it != this->energy_values.constEnd();
         ++it)
    {
        vals.insert(it.key(), it.value()[i]);
    }

    return vals;
}

/** Return the time and energy components at the ith row.
 *  Values are returned in internal units
 */
QHash<QString, double> EnergyTrajectory::get(int i) const
{
    return this->operator[](i);
}

/** Return the time and energy components at the ith row.
 *  Values are returned in the specified units
 */
QHash<QString, double> EnergyTrajectory::get(int i,
                                             const GeneralUnit &time_unit,
                                             const GeneralUnit &energy_unit) const
{
    i = SireID::Index(i).map(this->count());

    QHash<QString, double> vals;
    vals.reserve(1 + this->energy_values.count());

    vals.insert("time", Time(this->time_values[i]).to(time_unit));

    for (auto it = this->energy_values.constBegin();
         it != this->energy_values.constEnd();
         ++it)
    {
        vals.insert(it.key(), MolarEnergy(it.value()[i]).to(energy_unit));
    }

    return vals;
}

/** Return the labels at the ith row */
QHash<QString, QString> EnergyTrajectory::getLabels(int i) const
{
    i = SireID::Index(i).map(this->count());

    QHash<QString, QString> vals;
    vals.reserve(this->label_values.count());

    for (auto it = this->label_values.constBegin();
         it != this->label_values.constEnd();
         ++it)
    {
        vals.insert(it.key(), it.value()[i]);
    }

    return vals;
}

/** Return the labels at the ith row as numbers */
QHash<QString, double> EnergyTrajectory::getLabelsAsNumbers(int i) const
{
    i = SireID::Index(i).map(this->count());

    QHash<QString, double> vals;
    vals.reserve(this->label_values.count());

    for (auto it = this->label_values.constBegin();
         it != this->label_values.constEnd();
         ++it)
    {
        if (it.value()[i].isEmpty())
        {
            vals.insert(it.key(), NAN);
            continue;
        }

        double val;
        bool ok = false;

        val = it.value()[i].toDouble(&ok);

        if (not ok)
            throw SireError::invalid_cast(QObject::tr(
                                              "Cannot convert label '%1' at index %2 to a number. Value is %3.")
                                              .arg(it.key())
                                              .arg(i)
                                              .arg(it.value()[i]),
                                          CODELOC);

        vals.insert(it.key(), val);
    }

    return vals;
}

/** Return all of the energy keys */
QStringList EnergyTrajectory::keys() const
{
    return this->energy_values.keys();
}

/** Return all of the label keys */
QStringList EnergyTrajectory::labelKeys() const
{
    return this->label_values.keys();
}

/** Return all of the time values (the time column). This is
 *  sorted from earliest to latest time, and is in the default internal
 *  unit
 */
QVector<double> EnergyTrajectory::times() const
{
    return this->time_values;
}

/** Return all of the label values for the passed key (the label-key column).
 *  This is in the same order as the times.
 */
QVector<QString> EnergyTrajectory::labels(const QString &key) const
{
    auto it = this->label_values.constFind(key);

    if (it == this->label_values.constEnd())
        throw SireError::invalid_key(QObject::tr(
                                         "There is no label component with key '%1'. Valid "
                                         "keys are [ %2 ].")
                                         .arg(key)
                                         .arg(this->label_values.keys().join(", ")),
                                     CODELOC);

    return it.value();
}

/** Return all of the label values for the passed key (the label-key column)
 *  as numbers. This is in the same order as the times.
 */
QVector<double> EnergyTrajectory::labelsAsNumbers(const QString &key) const
{
    auto it = this->label_values.constFind(key);

    if (it == this->label_values.constEnd())
        throw SireError::invalid_key(QObject::tr(
                                         "There is no label component with key '%1'. Valid "
                                         "keys are [ %2 ].")
                                         .arg(key)
                                         .arg(this->label_values.keys().join(", ")),
                                     CODELOC);

    const auto &l = it.value();

    QVector<double> vals;
    vals.reserve(l.count());

    for (const auto &label : l)
    {
        if (label.isEmpty())
        {
            vals.append(NAN);
            continue;
        }

        bool ok = false;
        double val = label.toDouble(&ok);

        if (not ok)
            throw SireError::invalid_cast(QObject::tr(
                                              "Cannot convert label '%1' to a number. Value is %2.")
                                              .arg(key)
                                              .arg(label),
                                          CODELOC);

        vals.append(val);
    }

    return vals;
}

/** Return all of the
 *  energy values for the passed key (the energy-key column).
 *  This is in the same order as the times, and is in the default internal
 *  unit
 */
QVector<double> EnergyTrajectory::energies(const QString &key) const
{
    auto it = this->energy_values.constFind(key);

    if (it == this->energy_values.constEnd())
        throw SireError::invalid_key(QObject::tr(
                                         "There is no energy component with key '%1'. Valid "
                                         "keys are [ %2 ].")
                                         .arg(key)
                                         .arg(this->energy_values.keys().join(", ")),
                                     CODELOC);

    return it.value();
}

/** Return all of the times converted to the passed unit */
QVector<double> EnergyTrajectory::times(const SireUnits::Dimension::GeneralUnit &time_unit) const
{
    const double factor = 1.0 / Time(time_unit).value();

    if (factor == 1)
        return this->times();

    auto t = this->times();

    for (auto &val : t)
    {
        val *= factor;
    }

    return t;
}

/** Return all of the energies fro the passed key converted to the
 *  passed unit
 */
QVector<double> EnergyTrajectory::energies(const QString &key,
                                           const SireUnits::Dimension::GeneralUnit &energy_unit) const
{
    const double factor = 1.0 / MolarEnergy(energy_unit).value();

    if (factor == 1)
        return this->energies(key);

    auto e = this->energies(key);

    for (auto &val : e)
    {
        val *= factor;
    }

    return e;
}

/** Set a property on the trajectory. This is used to store additional
 *  information about the trajectory, such as the simulation temperature
 */
void EnergyTrajectory::setProperty(const QString &key, const SireBase::Property &value)
{
    props.setProperty(key, value);
}

/** Get a property on the trajectory. This could be additional
 *  information about the trajectory, such as the simulation temperature
 */
const SireBase::Property &EnergyTrajectory::property(const SireBase::PropertyName &key) const
{
    return props.property(key);
}

/** Return whether or not the trajectory has a property with the passed key */
bool EnergyTrajectory::hasProperty(const SireBase::PropertyName &key)
{
    return props.hasProperty(key);
}

/** Return all of the properties on the trajectory */
const SireBase::Properties &EnergyTrajectory::properties() const
{
    return props;
}

/** Return all of the property keys on the trajectory */
QStringList EnergyTrajectory::propertyKeys() const
{
    return props.propertyKeys();
}

/** Remove a property from the trajectory */
void EnergyTrajectory::removeProperty(const QString &key)
{
    props.removeProperty(key);
}

/** Clear all of the properties on the trajectory */
void EnergyTrajectory::clearProperties()
{
    props = Properties();
}
