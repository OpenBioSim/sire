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

#ifndef SIREMATHS_ENERGYTRAJECTORY_H
#define SIREMATHS_ENERGYTRAJECTORY_H

#include "SireUnits/generalunit.h"
#include "SireUnits/dimensions.h"

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
    class EnergyTrajectory;
}

SIREMATHS_EXPORT QDataStream &operator<<(QDataStream &ds, const SireMaths::EnergyTrajectory &);
SIREMATHS_EXPORT QDataStream &operator>>(QDataStream &ds, SireMaths::EnergyTrajectory &);

namespace SireMaths
{
    /** This class holds the trajectory of energies, organised by
     *  timestep the energy was recorded and the types of energy
     *  (e.g. kinetic, potential, values at different lambda windows)
     */
    class SIREMATHS_EXPORT EnergyTrajectory
        : public SireBase::ConcreteProperty<EnergyTrajectory, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const EnergyTrajectory &);
        friend QDataStream & ::operator>>(QDataStream &, EnergyTrajectory &);

    public:
        EnergyTrajectory();
        EnergyTrajectory(const EnergyTrajectory &other);

        ~EnergyTrajectory();

        EnergyTrajectory &operator=(const EnergyTrajectory &other);

        bool operator==(const EnergyTrajectory &other) const;
        bool operator!=(const EnergyTrajectory &other) const;

        virtual EnergyTrajectory *clone() const;

        static const char *typeName();
        const char *what() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;

        QHash<QString, double> operator[](int i) const;

        QHash<QString, double> get(int i) const;
        QHash<QString, double> get(int i,
                                   const SireUnits::Dimension::GeneralUnit &time_unit,
                                   const SireUnits::Dimension::GeneralUnit &energy_unit) const;

        QHash<QString, QString> getLabels(int i) const;
        QHash<QString, double> getLabelsAsNumbers(int i) const;

        void set(const SireUnits::Dimension::GeneralUnit &time,
                 const QHash<QString, SireUnits::Dimension::GeneralUnit> &energies);

        void set(const SireUnits::Dimension::GeneralUnit &time,
                 const QHash<QString, SireUnits::Dimension::GeneralUnit> &energies,
                 const QHash<QString, QString> &labels);

        QStringList keys() const;
        QStringList labelKeys() const;

        QVector<double> times() const;
        QVector<double> energies(const QString &key) const;
        QVector<QString> labels(const QString &key) const;
        QVector<double> labelsAsNumbers(const QString &key) const;

        QVector<double> times(const SireUnits::Dimension::GeneralUnit &time_unit) const;
        QVector<double> energies(const QString &key,
                                 const SireUnits::Dimension::GeneralUnit &energy_unit) const;

    private:
        /** All of the time values */
        QVector<double> time_values;

        /** All of the energy values */
        QHash<QString, QVector<double>> energy_values;

        /** All of the label values */
        QHash<QString, QVector<QString>> label_values;
    };
}

Q_DECLARE_METATYPE(SireMaths::EnergyTrajectory)

SIRE_EXPOSE_CLASS(SireMaths::EnergyTrajectory)

SIRE_END_HEADER

#endif
