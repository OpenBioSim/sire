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
        << static_cast<const Property &>(etraj);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, EnergyTrajectory &etraj)
{
    VersionID v = readHeader(ds, r_etraj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> etraj.time_values >> etraj.energy_values >> static_cast<Property &>(etraj);
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
      time_values(other.time_values), energy_values(other.energy_values)
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
        Property::operator=(other);
    }

    return *this;
}

bool EnergyTrajectory::operator==(const EnergyTrajectory &other) const
{
    return time_values == other.time_values and
           energy_values == other.energy_values;
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
    keys.insert(0, "time");

    parts.append(keys.join("\t"));

    if (n <= 10)
    {
        for (int i = 0; i < n; ++i)
        {
            const auto vals = this->get(i, GeneralUnit(picosecond), GeneralUnit(kcal_per_mol));

            QStringList v;

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

            QStringList v;

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

            QStringList v;

            for (const auto &key : keys)
            {
                v.append(QString::number(vals[key]));
            }

            parts.append(v.join("\t"));
        }
    }

    return QObject::tr("SelectorMol( size=%1\n%2\n)").arg(n).arg(parts.join("\n"));
}

/** Set the energies at time 'time' to the components contained
 *  in 'energies'
 */
void EnergyTrajectory::set(const GeneralUnit &time,
                           const QHash<QString, GeneralUnit> &energies)
{
    if (energies.isEmpty())
        return;

    auto t = Time(time);

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

    // find the index of the time, and add it as needed.
    // Times are sorted, so start from the end and work towards
    // the beginning (as we will normally be adding energies onto
    // the end)
    int idx = time_values.count();
    bool must_create = true;

    while (idx > 0)
    {
        if (t > time_values[idx - 1])
            break;

        else if (t == time_values[idx - 1])
        {
            must_create = false;
            idx = idx - 1;
            break;
        }

        idx -= 1;
    }

    if (idx >= time_values.count())
    {
        time_values.append(t.value());

        for (auto it = this->energy_values.begin();
             it != this->energy_values.end(); ++it)
        {
            it.value().append(NAN);
        }
    }
    else if (must_create)
    {
        time_values.insert(idx, t.value());

        for (auto it = this->energy_values.begin();
             it != this->energy_values.end(); ++it)
        {
            it.value().insert(idx, NAN);
        }
    }

    for (auto it = energies.constBegin();
         it != energies.constEnd(); ++it)
    {
        this->energy_values[it.key()][idx] = MolarEnergy(it.value()).value();
    }
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

/** Return all of the energy keys */
QStringList EnergyTrajectory::keys() const
{
    return this->energy_values.keys();
}

/** Return all of the time values (the time column). This is
 *  sorted from earliest to latest time, and is in the default internal
 *  unit
 */
QVector<double> EnergyTrajectory::times() const
{
    return this->time_values;
}

/** Return all of the energy values for the passed key (the energy-key column).
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
