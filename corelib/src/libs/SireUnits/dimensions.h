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

#ifndef SIREUNITS_DIMENSIONS_H
#define SIREUNITS_DIMENSIONS_H

#include <QDataStream>
#include <QString>

#include "SireMaths/constants.h"

#include <boost/predef.h>
#include <boost/static_assert.hpp>

SIRE_BEGIN_HEADER

namespace SireUnits
{
    namespace Dimension
    {
        class Unit;
        class GeneralUnit;
    } // namespace Dimension
} // namespace SireUnits

QDataStream &operator<<(QDataStream &, const SireUnits::Dimension::Unit &);
QDataStream &operator>>(QDataStream &, SireUnits::Dimension::Unit &);

namespace SireUnits
{

    using SireMaths::pi;

    namespace Dimension
    {

        class TempBase;

        /** This is the base class of all units - at its heart, this is
            just a scale factor - how many times the base unit is the
            current value.

            @author Christopher Woods
        */
        class SIREUNITS_EXPORT Unit
        {

            friend SIREUNITS_EXPORT QDataStream & ::operator<<(QDataStream &, const Unit &);
            friend SIREUNITS_EXPORT QDataStream & ::operator>>(QDataStream &, Unit &);

        public:
            ~Unit()
            {
            }

            operator double() const
            {
                return sclfac;
            }

            double value() const
            {
                return sclfac;
            }

            double scaleFactor() const
            {
                return sclfac;
            }

            double convertToInternal(double value) const
            {
                return value * sclfac;
            }

            double convertFromInternal(double value) const
            {
                return value / sclfac;
            }

        protected:
            Unit(double scale_factor) : sclfac(scale_factor)
            {
            }

            Unit(const TempBase &temperature);

            void setScale(double scale_factor)
            {
                sclfac = scale_factor;
            }

        private:
            double sclfac;
        };

        SIREUNITS_EXPORT QPair<double, QString> getUnitString(int M, int L, int T, int C, int t, int Q, int A);

        namespace detail
        {
            template <class T>
            struct _simple_typename
            {
                static const char *typeName()
                {
                    static QString typenam = QString("SireUnits::Dimension::PhysUnit<%1-%2-%3-%4-%5-%6-%7>")
                                                 .arg(T::MASS())
                                                 .arg(T::LENGTH())
                                                 .arg(T::TIME())
                                                 .arg(T::CHARGE())
                                                 .arg(T::TEMPERATURE())
                                                 .arg(T::QUANTITY())
                                                 .arg(T::ANGLE());

                    static auto s = typenam.toStdString();

                    return s.c_str();
                }
            };

            template <class T>
            struct _registered_typename
            {
                static const char *typeName()
                {
                    return QMetaType::typeName(qMetaTypeId<T>());
                }
            };
        } // namespace detail

        /** Construct a physical unit with the specified
            Mass, Length, Time, Charge, temperature,
            Quantity and Angle dimensions

            @author Christopher Woods
        */
        template <int M, int L, int T, int C, int t, int Q, int A>
        class PhysUnit : public Unit
        {
        public:
            PhysUnit() : Unit(0)
            {
            }

            explicit PhysUnit(double scale_factor) : Unit(scale_factor)
            {
            }

            PhysUnit(const QString &value);

            explicit PhysUnit(const TempBase &temperature) : Unit(temperature)
            {
                // this must be a Temperature!
                // BOOST_STATIC_ASSERT( t == 1 and M == 0 and L == 0 and
                //                      T == 0 and C == 0 and Q == 0 and A == 0);
            }

            PhysUnit(const GeneralUnit &other);

            PhysUnit(const PhysUnit<M, L, T, C, t, Q, A> &other) : Unit(other)
            {
            }

            ~PhysUnit()
            {
            }

            static const char *typeName()
            {
                typedef PhysUnit<M, L, T, C, t, Q, A> PhysUnitT;

                return std::conditional<QMetaTypeId2<PhysUnitT>::Defined,
                                        detail::_registered_typename<PhysUnitT>,
                                        detail::_simple_typename<PhysUnitT>>::type::typeName();
            }

            const char *what() const
            {
                return PhysUnit<M, L, T, C, t, Q, A>::typeName();
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator=(const PhysUnit<M, L, T, C, t, Q, A> &other)
            {
                Unit::setScale(other.scaleFactor());
                return *this;
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator=(const GeneralUnit &other);

            bool operator==(const PhysUnit<M, L, T, C, t, Q, A> &other) const
            {
                return scaleFactor() == other.scaleFactor();
            }

            bool operator!=(const PhysUnit<M, L, T, C, t, Q, A> &other) const
            {
                return scaleFactor() != other.scaleFactor();
            }

            PhysUnit<M, L, T, C, t, Q, A> operator-() const
            {
                return PhysUnit<M, L, T, C, t, Q, A>(-scaleFactor());
            }

            PhysUnit<M, L, T, C, t, Q, A> operator+(const PhysUnit<M, L, T, C, t, Q, A> &other) const
            {
                return PhysUnit<M, L, T, C, t, Q, A>(scaleFactor() + other.scaleFactor());
            }

            PhysUnit<M, L, T, C, t, Q, A> operator-(const PhysUnit<M, L, T, C, t, Q, A> &other) const
            {
                return PhysUnit<M, L, T, C, t, Q, A>(scaleFactor() - other.scaleFactor());
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator+=(const PhysUnit<M, L, T, C, t, Q, A> &other)
            {
                Unit::setScale(scaleFactor() + other.scaleFactor());
                return *this;
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator-=(const PhysUnit<M, L, T, C, t, Q, A> &other)
            {
                Unit::setScale(scaleFactor() - other.scaleFactor());
                return *this;
            }

            PhysUnit<M, L, T, C, t, Q, A> operator*(double val) const
            {
                return PhysUnit<M, L, T, C, t, Q, A>(scaleFactor() * val);
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator*=(double val)
            {
                Unit::setScale(scaleFactor() * val);
                return *this;
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator/=(double val)
            {
                Unit::setScale(scaleFactor() / val);
                return *this;
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator*=(int val)
            {
                Unit::setScale(scaleFactor() * val);
                return *this;
            }

            PhysUnit<M, L, T, C, t, Q, A> &operator/=(int val)
            {
                Unit::setScale(scaleFactor() / val);
                return *this;
            }

            PhysUnit<M, L, T, C, t, Q, A> operator/(double val) const
            {
                return PhysUnit<M, L, T, C, t, Q, A>(scaleFactor() / val);
            }

            PhysUnit<M, L, T, C, t, Q, A> operator*(int val) const
            {
                return PhysUnit<M, L, T, C, t, Q, A>(scaleFactor() * val);
            }

            PhysUnit<M, L, T, C, t, Q, A> operator/(int val) const
            {
                return PhysUnit<M, L, T, C, t, Q, A>(scaleFactor() / val);
            }

            template <int _M, int _L, int _T, int _C, int _t, int _Q, int _A>
            PhysUnit<M + _M, L + _L, T + _T, C + _C, t + _t, Q + _Q, A + _A> operator*(
                const PhysUnit<_M, _L, _T, _C, _t, _Q, _A> &other) const
            {
                return PhysUnit<M + _M, L + _L, T + _T, C + _C, t + _t, Q + _Q, A + _A>(scaleFactor() * other.scaleFactor());
            }

            template <int _M, int _L, int _T, int _C, int _t, int _Q, int _A>
            PhysUnit<M - _M, L - _L, T - _T, C - _C, t - _t, Q - _Q, A - _A> operator/(
                const PhysUnit<_M, _L, _T, _C, _t, _Q, _A> &other) const
            {
                return PhysUnit<M - _M, L - _L, T - _T, C - _C, t - _t, Q - _Q, A - _A>(scaleFactor() / other.scaleFactor());
            }

            PhysUnit<-M, -L, -T, -C, -t, -Q, -A> invert() const
            {
                return PhysUnit<-M, -L, -T, -C, -t, -Q, -A>(1.0 / scaleFactor());
            }

            double in(const PhysUnit<M, L, T, C, t, Q, A> &units) const
            {
                return units.convertFromInternal(*this);
            }

            double to(const PhysUnit<M, L, T, C, t, Q, A> &units) const
            {
                return this->in(units);
            }

            QString unitString() const
            {
                return SireUnits::Dimension::getUnitString(M, L, T, C, t, Q, A).second;
            }

            QString toString() const
            {
                auto u = SireUnits::Dimension::getUnitString(M, L, T, C, t, Q, A);

                double v = value() / u.first;

                if (u.second.startsWith("°"))
                    return QString("%1%2").arg(v).arg(u.second);
                else
                    return QString("%1 %2").arg(v).arg(u.second);
            }

            static int MASS()
            {
                return M;
            }

            static int LENGTH()
            {
                return L;
            }

            static int TIME()
            {
                return T;
            }

            static int CHARGE()
            {
                return C;
            }

            static int TEMPERATURE()
            {
                return t;
            }

            static int QUANTITY()
            {
                return Q;
            }

            static int ANGLE()
            {
                return A;
            }
        };

        template <int M, int L, int T, int C, int t, int Q, int A>
        PhysUnit<M, L, T, C, t, Q, A> operator*(double val, const PhysUnit<M, L, T, C, t, Q, A> &unit)
        {
            return PhysUnit<M, L, T, C, t, Q, A>(val * unit.scaleFactor());
        }

        template <int M, int L, int T, int C, int t, int Q, int A>
        PhysUnit<M, L, T, C, t, Q, A> operator*(int val, const PhysUnit<M, L, T, C, t, Q, A> &unit)
        {
            return PhysUnit<M, L, T, C, t, Q, A>(val * unit.scaleFactor());
        }

        template <int M, int L, int T, int C, int t, int Q, int A>
        PhysUnit<-M, -L, -T, -C, -t, -Q, -A> operator/(double val, const PhysUnit<M, L, T, C, t, Q, A> &unit)
        {
            return PhysUnit<-M, -L, -T, -C, -t, -Q, -A>(val / unit.scaleFactor());
        }

        template <int M, int L, int T, int C, int t, int Q, int A>
        PhysUnit<-M, -L, -T, -C, -t, -Q, -A> operator/(int val, const PhysUnit<M, L, T, C, t, Q, A> &unit)
        {
            return PhysUnit<-M, -L, -T, -C, -t, -Q, -A>(val / unit.scaleFactor());
        }

/** Typedef the various unit dimensions (including derived units) */
#ifndef SKIP_BROKEN_GCCXML_PARTS
        typedef PhysUnit<0, 0, 0, 0, 0, 0, 0> Dimensionless;

        typedef PhysUnit<1, 0, 0, 0, 0, 0, 0> Mass;

        typedef PhysUnit<1, 0, 0, 0, 0, -1, 0> MolarMass;

        typedef PhysUnit<0, 1, 0, 0, 0, 0, 0> Length;

        typedef PhysUnit<0, 0, 1, 0, 0, 0, 0> Time;

        typedef PhysUnit<0, 0, 0, 1, 0, 0, 0> Charge;

        typedef PhysUnit<0, 0, 0, 1, 0, -1, 0> MolarCharge;

        typedef PhysUnit<0, 0, 0, 0, 1, 0, 0> Temperature;

        typedef PhysUnit<0, 0, 0, 0, 0, 1, 0> Quantity;

        typedef PhysUnit<0, 0, 0, 0, 0, 0, 1> Angle;

        typedef PhysUnit<0, 2, 0, 0, 0, 0, 0> Area;

        typedef PhysUnit<0, 3, 0, 0, 0, 0, 0> Volume;

        typedef PhysUnit<0, 3, 0, 0, 0, -1, 0> MolarVolume;

        typedef PhysUnit<0, 1, -1, 0, 0, 0, 0> Velocity;

        typedef PhysUnit<0, 0, -1, 0, 0, 0, 1> AngularVelocity;

        typedef PhysUnit<0, 1, -2, 0, 0, 0, 0> Acceleration;

        typedef PhysUnit<0, 0, -2, 0, 0, 0, 1> AngularAcceleration;

        typedef PhysUnit<1, 2, -2, 0, 0, 0, 0> Energy;

        typedef PhysUnit<1, 2, -2, 0, 0, -1, 0> MolarEnergy;

        typedef PhysUnit<1, 0, -2, 0, 0, -1, 0> HarmonicBondConstant;

        typedef PhysUnit<1, 2, -2, 0, 0, -1, -2> HarmonicAngleConstant;

        typedef PhysUnit<1, 2, -3, 0, 0, 0, 0> Power;

        typedef PhysUnit<1, 2, -3, 0, 0, -1, 0> MolarPower;

        typedef PhysUnit<1, -3, 0, 0, 0, 0, 0> Density;

        typedef PhysUnit<1, -3, 0, 0, 0, -1, 0> MolarDensity;

        typedef PhysUnit<1, 1, -2, 0, 0, 0, 0> Force;

        typedef PhysUnit<1, -1, -2, 0, 0, 0, 0> Pressure;

        typedef PhysUnit<0, 0, -1, 1, 0, 0, 0> Current;
        typedef PhysUnit<-1, -2, 2, 2, 0, 0, 0> Capacitance;
        typedef PhysUnit<1, 2, -2, -1, 0, 0, 0> Potential;

        typedef PhysUnit<-1, -3, 2, 2, 0, 0, 0> Constant1;
        typedef PhysUnit<1, 2, -1, 0, 0, 0, 0> Constant2;
        typedef PhysUnit<1, 1, 0, -2, 0, 0, 0> Constant3;
        typedef PhysUnit<-1, 3, -2, 0, 0, 0, 0> Constant4;
        typedef PhysUnit<1, 2, -2, 0, -1, -1, 0> Constant5;

#else // else with 'ifndef SKIP_BROKEN_GCCXML_PARTS'

        class Dimensionless;
        class Mass;
        class MolarMass;
        class Length;
        class Time;
        class Charge;
        class MolarCharge;
        class Temperature;
        class Quantity;
        class Angle;
        class Area;
        class Volume;
        class MolarVolume;
        class Velocity;
        class Acceleration;
        class Energy;
        class MolarEnergy;
        class Power;
        class MolarPower;
        class Density;
        class MolarDensity;
        class Force;
        class Pressure;
        class Capacitance;
        class Current;
        class Potential;
        class HarmonicBondConstant;
        class HarmonicAngleConstant;

        class Constant1;
        class Constant2;
        class Constant3;
        class Constant4;
        class Constant5;

#endif // end of 'ifndef SKIP_BROKEN_GCCXML_PARTS'

    } // end of namespace Dimension

} // namespace SireUnits

/** Serialise a unit to a binary datastream (this does not check
    the type of unit!) */
SIRE_ALWAYS_INLINE QDataStream &operator<<(QDataStream &ds, const SireUnits::Dimension::Unit &unit)
{
    ds << unit.sclfac;
    return ds;
}

/** Extract from a binary datastream (this does not check the type
    of unit!) */
SIRE_ALWAYS_INLINE QDataStream &operator>>(QDataStream &ds, SireUnits::Dimension::Unit &unit)
{
    ds >> unit.sclfac;
    return ds;
}

Q_DECLARE_METATYPE(SireUnits::Dimension::Dimensionless);
Q_DECLARE_METATYPE(SireUnits::Dimension::Mass);
Q_DECLARE_METATYPE(SireUnits::Dimension::MolarMass);
Q_DECLARE_METATYPE(SireUnits::Dimension::Length);
Q_DECLARE_METATYPE(SireUnits::Dimension::Time);
Q_DECLARE_METATYPE(SireUnits::Dimension::Charge);
Q_DECLARE_METATYPE(SireUnits::Dimension::MolarCharge);
Q_DECLARE_METATYPE(SireUnits::Dimension::Temperature);
Q_DECLARE_METATYPE(SireUnits::Dimension::Quantity);
Q_DECLARE_METATYPE(SireUnits::Dimension::Angle);
Q_DECLARE_METATYPE(SireUnits::Dimension::Area);
Q_DECLARE_METATYPE(SireUnits::Dimension::Volume);
Q_DECLARE_METATYPE(SireUnits::Dimension::MolarVolume);
Q_DECLARE_METATYPE(SireUnits::Dimension::Velocity);
Q_DECLARE_METATYPE(SireUnits::Dimension::AngularVelocity);
Q_DECLARE_METATYPE(SireUnits::Dimension::Acceleration);
Q_DECLARE_METATYPE(SireUnits::Dimension::AngularAcceleration);
Q_DECLARE_METATYPE(SireUnits::Dimension::Energy);
Q_DECLARE_METATYPE(SireUnits::Dimension::MolarEnergy);
Q_DECLARE_METATYPE(SireUnits::Dimension::Power);
Q_DECLARE_METATYPE(SireUnits::Dimension::MolarPower);
Q_DECLARE_METATYPE(SireUnits::Dimension::Density);
Q_DECLARE_METATYPE(SireUnits::Dimension::MolarDensity);
Q_DECLARE_METATYPE(SireUnits::Dimension::Force);
Q_DECLARE_METATYPE(SireUnits::Dimension::Pressure);
Q_DECLARE_METATYPE(SireUnits::Dimension::Capacitance);
Q_DECLARE_METATYPE(SireUnits::Dimension::Current);
Q_DECLARE_METATYPE(SireUnits::Dimension::Potential);
Q_DECLARE_METATYPE(SireUnits::Dimension::HarmonicBondConstant);
Q_DECLARE_METATYPE(SireUnits::Dimension::HarmonicAngleConstant);

Q_DECLARE_METATYPE(SireUnits::Dimension::Constant1);
Q_DECLARE_METATYPE(SireUnits::Dimension::Constant2);
Q_DECLARE_METATYPE(SireUnits::Dimension::Constant3);
Q_DECLARE_METATYPE(SireUnits::Dimension::Constant4);
Q_DECLARE_METATYPE(SireUnits::Dimension::Constant5);

SIRE_EXPOSE_CLASS(SireUnits::Dimension::Unit)

SIRE_END_HEADER

#endif
