
#include <QHash>
#include <QMutex>
#include <QStringList>

#include "generalunit.h"

#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"

#include "tostring.h"

#include <QDebug>

using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<GeneralUnit> r_genunit(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const SireUnits::Dimension::GeneralUnit &u)
{
    writeHeader(ds, r_genunit, 1);

    SharedDataStream sds(ds);

    qint8 a = u.ANGLE();
    qint8 c = u.CHARGE();
    qint8 l = u.LENGTH();
    qint8 m = u.MASS();
    qint8 t1 = u.TEMPERATURE();
    qint8 t2 = u.TIME();
    qint8 q = u.QUANTITY();

    sds << u.comps << a << c << l << m << t1 << t2 << q << static_cast<const Unit &>(u);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SireUnits::Dimension::GeneralUnit &u)
{
    VersionID v = readHeader(ds, r_genunit);

    if (v == 1)
    {
        qint8 a, c, l, m, t1, t2, q;

        SharedDataStream sds(ds);

        sds >> u.comps >> a >> c >> l >> m >> t1 >> t2 >> q >> static_cast<Unit &>(u);

        u.Angle = a;
        u.Charge = c;
        u.Length = l;
        u.Mass = m;
        u.temperature = t1;
        u.Time = t2;
        u.Quantity = q;
        u.log_value = 0;
    }
    else
        throw version_error(v, "1", r_genunit, CODELOC);

    return ds;
}

namespace SireUnits
{
    namespace Dimension
    {
        namespace detail
        {

            typedef QHash<QString, QString> TypenameRegistry;

            Q_GLOBAL_STATIC(TypenameRegistry, typename_registry);
            Q_GLOBAL_STATIC(QMutex, registry_mutex);

            static QString getKey(const GeneralUnit &unit)
            {
                return QString("%1-%2-%3-%4-%5-%6-%7")
                    .arg(unit.MASS())
                    .arg(unit.LENGTH())
                    .arg(unit.TIME())
                    .arg(unit.CHARGE())
                    .arg(unit.TEMPERATURE())
                    .arg(unit.QUANTITY())
                    .arg(unit.ANGLE());
            }

            void SIREUNITS_EXPORT registerTypeName(const GeneralUnit &unit, const char *typnam)
            {
                QString key = getKey(unit);

                QMutexLocker lkr(registry_mutex());
                if (not typename_registry()->contains(key))
                {
                    typename_registry()->insert(key, QString(typnam));
                }
            }

            static QString getTypeName(const GeneralUnit &unit)
            {
                QString key = getKey(unit);

                QMutexLocker lkr(registry_mutex());

                if (typename_registry()->contains(key))
                    return typename_registry()->value(key);
                else
                    return QString("SireUnits::Dimension::PhysUnit<%1,%2,%3,%4,%5,%6,%7>")
                        .arg(unit.MASS())
                        .arg(unit.LENGTH())
                        .arg(unit.TIME())
                        .arg(unit.CHARGE())
                        .arg(unit.TEMPERATURE())
                        .arg(unit.QUANTITY())
                        .arg(unit.ANGLE());
            }

        } // end of namespace detail
    }     // end of namespace Dimension
} // end of namespace SireUnits

GeneralUnit::GeneralUnit() : Unit(0)
{
    Mass = 0;
    Length = 0;
    Time = 0;
    Charge = 0;
    temperature = 0;
    Quantity = 0;
    Angle = 0;
    log_value = 0;
}

GeneralUnit::GeneralUnit(double value) : Unit(value)
{
    Mass = 0;
    Length = 0;
    Time = 0;
    Charge = 0;
    temperature = 0;
    Quantity = 0;
    Angle = 0;
    log_value = 0;
}

GeneralUnit::GeneralUnit(const TempBase &t) : Unit(t)
{
    Mass = 0;
    Length = 0;
    Time = 0;
    Charge = 0;
    temperature = 1;
    Quantity = 0;
    Angle = 0;
    log_value = 0;
}

GeneralUnit::GeneralUnit(double value, const QList<qint32> &dimensions)
    : Unit(value)
{
    if (dimensions.count() != 7)
    {
        throw SireError::invalid_arg(QObject::tr(
                                         "You can only create a GeneralUnit using 7 dimensions "
                                         "(M,L,T,C,t,Q,A). You have only supplied %1.")
                                         .arg(Sire::toString(dimensions)),
                                     CODELOC);
    }

    Mass = dimensions[0];
    Length = dimensions[1];
    Time = dimensions[2];
    Charge = dimensions[3];
    temperature = dimensions[4];
    Quantity = dimensions[5];
    Angle = dimensions[6];
    log_value = 0;
}

GeneralUnit::GeneralUnit(const GeneralUnit &other) : Unit(other), comps(other.comps)
{
    Mass = other.Mass;
    Length = other.Length;
    Time = other.Time;
    Charge = other.Charge;
    temperature = other.temperature;
    Quantity = other.Quantity;
    Angle = other.Angle;
    log_value = other.log_value;
}

GeneralUnit::~GeneralUnit()
{
}

const char *GeneralUnit::typeName()
{
    return QMetaType::typeName(qMetaTypeId<GeneralUnit>());
}

const char *GeneralUnit::what() const
{
    return GeneralUnit::typeName();
}

/** Return the C++ type that this particular GeneralUnit corresponds to */
QString GeneralUnit::_to_cpp_type() const
{
    return detail::getTypeName(*this);
}

void GeneralUnit::assertCompatible(const GeneralUnit &other) const
{
    if (this->isZero() or other.isZero())
        return;

    if (Mass != other.Mass or Length != other.Length or Time != other.Time or Charge != other.Charge or
        temperature != other.temperature or Quantity != other.Quantity or Angle != other.Angle)
    {
        throw SireError::incompatible_error(
            QObject::tr("Units for values %1 and %2 are incompatible.").arg(this->toString()).arg(other.toString()));
    }
}

QString GeneralUnit::unitString() const
{
    return SireUnits::Dimension::getUnitString(Mass, Length, Time, Charge, temperature, Quantity, Angle).second;
}

QString GeneralUnit::toString() const
{
    auto u = SireUnits::Dimension::getUnitString(Mass, Length, Time, Charge, temperature, Quantity, Angle);

    double v = value();

    if (u.first != 0)
        v /= u.first;

    if (u.second.startsWith("°"))
        return QString("%1%2").arg(v).arg(u.second);
    else
        return QString("%1 %2").arg(v).arg(u.second);
}

/** Return the physical dimensions of this unit, in the order
 *  (M,L,T,C,t,Q,A)
 */
QList<qint32> GeneralUnit::dimensions() const
{
    return QList<qint32>({Mass, Length, Time, Charge, temperature, Quantity, Angle});
}

double GeneralUnit::to(const GeneralUnit &units) const
{
    assertCompatible(units);
    return units.convertFromInternal(value());
}

double GeneralUnit::to(const TempBase &other) const
{
    // this must be a temperature!
    GeneralUnit general_temp;
    general_temp.temperature = 1;
    general_temp.setScale(other);

    this->assertCompatible(general_temp);

    return other.convertFromInternal(value()) / other.convertFromInternal();
}

bool GeneralUnit::isDimensionless() const
{
    return Mass == 0 and Length == 0 and Time == 0 and Charge == 0 and temperature == 0 and Quantity == 0 and
           Angle == 0;
}

// test for zero when multiplying
bool _isZero(const double v)
{
    // smallest double is ~1e-308
    return std::abs(v) < 1e-307;
}

// test for zero when doing addition
bool _isWithinEpsilonZero(const double v)
{
    return std::abs(v) < std::numeric_limits<double>::epsilon();
}

bool GeneralUnit::isWithinEpsilonZero() const
{
    return _isWithinEpsilonZero(this->value());
}

bool GeneralUnit::isZero() const
{
    return _isZero(this->value());
}

int GeneralUnit::MASS() const
{
    return Mass;
}

int GeneralUnit::LENGTH() const
{
    return Length;
}

int GeneralUnit::TIME() const
{
    return Time;
}

int GeneralUnit::CHARGE() const
{
    return Charge;
}

int GeneralUnit::TEMPERATURE() const
{
    return temperature;
}

int GeneralUnit::QUANTITY() const
{
    return Quantity;
}

int GeneralUnit::ANGLE() const
{
    return Angle;
}

/** Return the units of this value */
GeneralUnit GeneralUnit::units() const
{
    GeneralUnit ret(*this);
    ret.setScale(1.0);
    ret.log_value = 0;
    return ret;
}

bool GeneralUnit::hasSameUnits(const GeneralUnit &other) const
{
    return Mass == other.Mass and Length == other.Length and Time == other.Time and Charge == other.Charge and
           temperature == other.temperature and Quantity == other.Quantity and Angle == other.Angle;
}

GeneralUnit &GeneralUnit::operator=(const GeneralUnit &other)
{
    setScale(other.value());

    Mass = other.MASS();
    Length = other.LENGTH();
    Time = other.TIME();
    Charge = other.CHARGE();
    temperature = other.TEMPERATURE();
    Quantity = other.QUANTITY();
    Angle = other.ANGLE();
    log_value = other.log_value;
    comps = other.comps;

    return *this;
}

bool GeneralUnit::operator==(const GeneralUnit &other) const
{
    // note that we don't compare the componetns because
    // this leads to surprising behaviour where
    // two totals that are equal are not the same

    assertCompatible(other);
    return value() == other.value();
}

bool GeneralUnit::operator!=(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() != other.value();
}

bool GeneralUnit::operator>(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() > other.value();
}

bool GeneralUnit::operator>=(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() >= other.value();
}

bool GeneralUnit::operator<(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() < other.value();
}

bool GeneralUnit::operator<=(const GeneralUnit &other) const
{
    assertCompatible(other);
    return value() <= other.value();
}

GeneralUnit GeneralUnit::operator-() const
{
    GeneralUnit ret = *this;
    ret.setScale(-value());

    for (const auto &key : this->comps.keys())
    {
        ret.comps[key] = -(this->comps[key]);
    }

    return ret;
}

GeneralUnit &GeneralUnit::operator+=(const TempBase &other)
{
    return this->operator+=(GeneralUnit(other));
}

GeneralUnit &GeneralUnit::operator-=(const TempBase &other)
{
    return this->operator-=(GeneralUnit(other));
}

GeneralUnit &GeneralUnit::operator*=(const TempBase &other)
{
    return this->operator*=(GeneralUnit(other));
}

GeneralUnit &GeneralUnit::operator/=(const TempBase &other)
{
    return this->operator/=(GeneralUnit(other));
}

GeneralUnit GeneralUnit::operator+(const TempBase &other) const
{
    return this->operator+(GeneralUnit(other));
}

GeneralUnit GeneralUnit::operator-(const TempBase &other) const
{
    return this->operator-(GeneralUnit(other));
}

GeneralUnit GeneralUnit::operator*(const TempBase &other) const
{
    return this->operator*(GeneralUnit(other));
}

GeneralUnit GeneralUnit::operator/(const TempBase &other) const
{
    return this->operator/(GeneralUnit(other));
}

GeneralUnit &GeneralUnit::operator+=(const GeneralUnit &other)
{
    if (this->isWithinEpsilonZero())
    {
        this->operator=(other);
        return *this;
    }
    else if (other.isWithinEpsilonZero())
    {
        return *this;
    }

    assertCompatible(other);
    setScale(value() + other.value());
    log_value = 0;

    for (const auto &key : other.comps.keys())
    {
        this->comps[key] = this->comps.value(key) + other.comps[key];

        if (_isWithinEpsilonZero(this->comps[key]))
        {
            this->comps.remove(key);
        }
    }

    if (this->isWithinEpsilonZero() and this->comps.isEmpty())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit &GeneralUnit::operator-=(const GeneralUnit &other)
{
    if (this->isWithinEpsilonZero())
    {
        this->operator=(-other);
        return *this;
    }
    else if (other.isWithinEpsilonZero())
    {
        return *this;
    }

    assertCompatible(other);
    setScale(value() - other.value());
    log_value = 0;

    for (const auto &key : other.comps.keys())
    {
        this->comps[key] = this->comps.value(key) - other.comps[key];

        if (_isWithinEpsilonZero(this->comps[key]))
        {
            this->comps.remove(key);
        }
    }

    if (this->isWithinEpsilonZero() and this->comps.isEmpty())
        this->operator=(GeneralUnit());

    return *this;
}

GeneralUnit GeneralUnit::operator+(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret += other;
    return ret;
}

GeneralUnit GeneralUnit::operator-(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret -= other;
    return ret;
}

GeneralUnit &GeneralUnit::operator+=(double val)
{
    return this->operator+=(GeneralUnit(val));
}

GeneralUnit &GeneralUnit::operator-=(double val)
{
    return this->operator-=(GeneralUnit(val));
}

GeneralUnit GeneralUnit::operator+(double val) const
{
    GeneralUnit ret = *this;
    ret += val;
    return ret;
}

GeneralUnit GeneralUnit::operator-(double val) const
{
    GeneralUnit ret = *this;
    ret -= val;
    return ret;
}

bool GeneralUnit::useLog() const
{
    if (log_value != 0)
        return true;

    double v = std::abs(this->value());

    return (v > 1e6 or v < 1e-6);
}

double GeneralUnit::getLog() const
{
    if (log_value != 0)
        return log_value;
    else
    {
        double v = std::abs(this->value());

        if (v == 0)
            return 0;
        else
            return std::log10(v);
    }
}

double GeneralUnit::getLog()
{
    if (log_value == 0)
    {
        double v = std::abs(this->value());

        if (v != 0)
            log_value = std::log10(v);
    }

    return log_value;
}

GeneralUnit GeneralUnit::pow(int n) const
{
    if (n == 1)
    {
        return *this;
    }
    else if (n == 0)
    {
        return GeneralUnit(1.0);
    }
    else if (this->isZero())
    {
        return GeneralUnit(0);
    }

    GeneralUnit ret(*this);

    ret.comps.clear();

    double log_pow_n = ret.getLog() * n;

    ret.log_value = log_pow_n;

    if (this->value() > 0)
        ret.setScale(std::pow(10, log_pow_n));
    else
        ret.setScale(-std::pow(10, log_pow_n));

    ret.Mass *= n;
    ret.Length *= n;
    ret.Time *= n;
    ret.Charge *= n;
    ret.temperature *= n;
    ret.Quantity *= n;
    ret.Angle *= n;

    return ret;
}

GeneralUnit &GeneralUnit::operator*=(const GeneralUnit &other)
{
    if (this->isZero())
    {
        return *this;
    }
    else if (other.isZero())
    {
        this->operator=(GeneralUnit(0));
        return *this;
    }

    if (this->useLog() or other.useLog())
    {
        double log_prod = this->getLog() + other.getLog();
        log_value = log_prod;

        if ((this->value() > 0 and other.value() > 0) or (this->value() < 0 and other.value() < 0))
        {
            this->setScale(std::pow(10.0, log_value));
        }
        else
        {
            this->setScale(-std::pow(10.0, log_value));
        }
    }
    else
    {
        setScale(value() * other.value());
        log_value = 0;
    }

    Mass += other.Mass;
    Length += other.Length;
    Time += other.Time;
    Charge += other.Charge;
    temperature += other.temperature;
    Quantity += other.Quantity;
    Angle += other.Angle;

    if (other.comps.isEmpty())
    {
        // we can just multiply every component by the value
        for (const auto &key : this->comps.keys())
        {
            this->comps[key] *= other.value();

            if (_isWithinEpsilonZero(this->comps[key]))
                this->comps.remove(key);
        }
    }
    else if (this->comps.isEmpty())
    {
        // we can multiply by our constant
        const double v = this->value();

        this->comps = other.comps;

        for (const auto &key : this->comps.keys())
        {
            this->comps[key] *= v;

            if (_isWithinEpsilonZero(this->comps[key]))
                this->comps.remove(key);
        }
    }
    else
    {
        // we can't retain the components as they won't multiply
        this->comps.clear();
    }

    return *this;
}

GeneralUnit &GeneralUnit::operator/=(const GeneralUnit &other)
{
    if (this->isZero())
    {
        return *this;
    }
    else if (other.isZero())
    {
        // this will become infinite
        throw SireError::invalid_operation(QObject::tr(
                                               "You cannot divide %1 by 0 (%2)")
                                               .arg(this->toString())
                                               .arg(other.toString()),
                                           CODELOC);
    }

    if (this->useLog() or other.useLog())
    {
        double log_div = this->getLog() - other.getLog();
        log_value = log_div;

        if ((this->value() > 0 and other.value() > 0) or (this->value() < 0 and other.value() < 0))
        {
            this->setScale(std::pow(10.0, log_value));
        }
        else
        {
            this->setScale(-std::pow(10.0, log_value));
        }
    }
    else
    {
        setScale(value() / other.value());
        log_value = 0;
    }

    Mass -= other.Mass;
    Length -= other.Length;
    Time -= other.Time;
    Charge -= other.Charge;
    temperature -= other.temperature;
    Quantity -= other.Quantity;
    Angle -= other.Angle;

    if (other.comps.isEmpty())
    {
        // we can just multiply every component by the value
        for (const auto &key : this->comps.keys())
        {
            this->comps[key] /= other.value();

            if (_isWithinEpsilonZero(this->comps[key]))
                this->comps.remove(key);
        }
    }
    else if (this->comps.isEmpty())
    {
        // we can multiply by our constant
        const double v = this->value();

        this->comps = other.comps;

        for (const auto &key : this->comps.keys())
        {
            this->comps[key] /= v;

            if (_isWithinEpsilonZero(this->comps[key]))
                this->comps.remove(key);
        }
    }
    else
    {
        // we can't retain the components as they won't multiply
        this->comps.clear();
    }

    return *this;
}

SIREUNITS_EXPORT GeneralUnit Celsius::operator+(const GeneralUnit &other) const
{
    return other.operator+(*this);
}

SIREUNITS_EXPORT GeneralUnit Celsius::operator-(const GeneralUnit &other) const
{
    return GeneralUnit(*this) - other;
}

SIREUNITS_EXPORT GeneralUnit Celsius::operator*(const GeneralUnit &other) const
{
    return other.operator*(*this);
}

SIREUNITS_EXPORT GeneralUnit Celsius::operator/(const GeneralUnit &other) const
{
    return GeneralUnit(*this) / other;
}

SIREUNITS_EXPORT GeneralUnit Fahrenheit::operator+(const GeneralUnit &other) const
{
    return other.operator+(*this);
}

SIREUNITS_EXPORT GeneralUnit Fahrenheit::operator-(const GeneralUnit &other) const
{
    return GeneralUnit(*this) - other;
}

SIREUNITS_EXPORT GeneralUnit Fahrenheit::operator*(const GeneralUnit &other) const
{
    return other.operator*(*this);
}

SIREUNITS_EXPORT GeneralUnit Fahrenheit::operator/(const GeneralUnit &other) const
{
    return GeneralUnit(*this) / other;
}

GeneralUnit GeneralUnit::operator*(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret *= other;

    return ret;
}

GeneralUnit GeneralUnit::operator/(const GeneralUnit &other) const
{
    GeneralUnit ret = *this;
    ret /= other;
    return ret;
}

GeneralUnit &GeneralUnit::operator*=(double val)
{
    if (val == 0 or this->isZero())
    {
        this->operator=(GeneralUnit(0));
        return *this;
    }

    if (this->useLog())
    {
        double log_prod = this->getLog();
        log_prod += std::log10(std::abs(val));

        if ((this->value() < 0 and val < 0) or (this->value() > 0 and val > 0))
        {
            setScale(std::pow(10.0, log_prod));
        }
        else
        {
            setScale(-std::pow(10.0, log_prod));
        }

        log_value = log_prod;
    }
    else
    {
        setScale(value() * val);
        log_value = 0;
    }

    for (const auto &key : this->comps.keys())
    {
        this->comps[key] *= val;

        if (_isZero(this->comps[key]))
            this->comps.remove(key);
    }

    return *this;
}

GeneralUnit &GeneralUnit::operator/=(double val)
{
    if (this->isZero())
    {
        this->operator=(GeneralUnit(0));
        return *this;
    }
    else if (val == 0)
    {
        // this will become infinite
        throw SireError::invalid_operation(QObject::tr(
                                               "You cannot divide %1 by 0")
                                               .arg(this->toString()),
                                           CODELOC);
    }

    if (this->useLog())
    {
        double log_div = this->getLog() - std::log10(std::abs(val));

        if ((this->value() < 0 and val < 0) or (this->value() > 0 and val > 0))
        {
            setScale(std::pow(10.0, log_div));
        }
        else
        {
            setScale(-std::pow(10.0, log_div));
        }

        log_value = log_div;
    }
    else
    {
        setScale(value() / val);
        log_value = 0;
    }

    for (const auto &key : this->comps.keys())
    {
        this->comps[key] /= val;

        if (_isZero(this->comps[key]))
            this->comps.remove(key);
    }

    return *this;
}

GeneralUnit &GeneralUnit::operator*=(int val)
{
    return this->operator*=(double(val));
}

GeneralUnit &GeneralUnit::operator/=(int val)
{
    return this->operator/=(double(val));
}

GeneralUnit GeneralUnit::operator*(double val) const
{
    GeneralUnit ret = *this;
    ret *= val;
    return ret;
}

GeneralUnit GeneralUnit::operator/(double val) const
{
    GeneralUnit ret = *this;
    ret /= val;
    return ret;
}

GeneralUnit GeneralUnit::operator*(int val) const
{
    GeneralUnit ret = *this;
    ret *= val;
    return ret;
}

GeneralUnit GeneralUnit::operator/(int val) const
{
    GeneralUnit ret = *this;
    ret /= val;
    return ret;
}

GeneralUnit GeneralUnit::invert() const
{
    return this->pow(-1);
}

QHash<QString, GeneralUnit> GeneralUnit::components() const
{
    QHash<QString, GeneralUnit> c;
    c.reserve(this->comps.count());

    for (const auto &key : this->comps.keys())
    {
        GeneralUnit v(*this);
        v.comps = QHash<QString, double>();
        v.setScale(this->comps[key]);
        c.insert(key, v);
    }

    return c;
}

GeneralUnit GeneralUnit::getComponent(const QString &component) const
{
    GeneralUnit v(*this);
    v.comps = QHash<QString, double>();
    v.setScale(this->comps.value(component));
    return v;
}

void GeneralUnit::setComponent(const QString &component, const GeneralUnit &value)
{
    this->assertCompatible(value);

    if (this->isWithinEpsilonZero() and this->comps.isEmpty())
    {
        this->operator=(value);
        this->comps.clear();
        this->comps.insert(component, value.value());
        log_value = 0;
        return;
    }

    double delta = value.value() - this->comps.value(component);

    this->comps[component] = value.value();

    if (_isWithinEpsilonZero(this->comps[component]))
    {
        this->comps.remove(component);
    }

    this->setScale(this->value() + delta);

    log_value = 0;

    if (this->isWithinEpsilonZero() and this->comps.isEmpty())
        this->operator=(GeneralUnit());
}

void GeneralUnit::addComponent(const QString &component, const GeneralUnit &value)
{
    if (value.isWithinEpsilonZero())
        return;

    if (this->isWithinEpsilonZero() and this->comps.isEmpty())
    {
        this->operator=(value);
        this->comps.insert(component, value.value());
        log_value = 0;
        return;
    }

    this->assertCompatible(value);

    this->comps[component] += value.value();

    if (_isWithinEpsilonZero(this->comps[component]))
    {
        this->comps.remove(component);
    }

    this->setScale(this->value() + value.value());

    log_value = 0;

    if (this->isWithinEpsilonZero() and this->comps.isEmpty())
        this->operator=(GeneralUnit());
}

void GeneralUnit::subtractComponent(const QString &component, const GeneralUnit &value)
{
    this->addComponent(component, -value);
}
