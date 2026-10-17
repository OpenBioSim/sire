/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
  *
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published
  *  by the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  *
  *  For full details of the license please see the LICENSE file
  *  that should have come with this distribution.
  *
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include <QDataStream>

#include "complex.h"
#include "rational.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireMaths;

static const RegisterMetaType<Complex> r_complex(NO_ROOT);

/** Serialise a Complex to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const SireMaths::Complex &z)
{
    writeHeader(ds, r_complex, 1) << z.real() << z.imag();

    return ds;
}

/** Deserialise a Complex from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SireMaths::Complex &z)
{
    VersionID v = readHeader(ds, r_complex);

    if (v == 1)
    {
        double r, i;
        ds >> r >> i;
        z.setRectangular(r, i);
    }
    else
        throw version_error(v, "1", r_complex, CODELOC);

    return ds;
}

/** Construct the complex number  real + imag i */
Complex::Complex(double r, double i) : val(r, i)
{
}

/** Copy constructor */
Complex::Complex(const Complex &other) : val(other.val)
{
}

/** Destructor */
Complex::~Complex()
{
}

/** Return the real part of this number */
double Complex::real() const
{
    return val.real();
}

/** Return the imaginary part of this number */
double Complex::imag() const
{
    return val.imag();
}

/** Is this a real number (imag == 0) ? */
bool Complex::isReal() const
{
    return SireMaths::isZero(imag());
}

/** Is this a pure complex number (real == 0) */
bool Complex::isPurelyComplex() const
{
    return SireMaths::isZero(real());
}

/** Is this zero? */
bool Complex::isZero() const
{
    return isReal() and SireMaths::isZero(real());
}

/** Return a string representation of this Complex number */
QString Complex::toString() const
{
    if (isReal())
        return QString::number(real());
    else
    {
        double i = imag();
        double r = real();

        if (SireMaths::isZero(r))
        {
            if (SireMaths::areEqual(i, 1.0))
                return "i";
            else if (SireMaths::areEqual(i, -1.0))
                return "-i";
            else
                return QString("%1i").arg(i);
        }
        else if (i > 0)
            return QString("%1 + %2i").arg(r).arg(i);
        else
            return QString("%1 - %2i").arg(r).arg(-i);
    }
}

/** Return the complex number x + iy, given its rectangular components */
Complex Complex::rect(double x, double y)
{
    return Complex(x, y);
}

/** Return the complex number with modulus r and argument theta */
Complex Complex::polar(double r, double theta)
{
    return Complex(std::polar(r, theta));
}

/** Set this number to x + iy */
void Complex::setRectangular(double x, double y)
{
    val = std::complex<double>(x, y);
}

/** Set this number from its modulus r and argument theta */
void Complex::setPolar(double r, double theta)
{
    val = std::polar(r, theta);
}

/** Set the real part, keeping the imaginary part */
void Complex::setReal(double x)
{
    val.real(x);
}

/** Set the imaginary part, keeping the real part */
void Complex::setImag(double y)
{
    val.imag(y);
}

/** Comparison operator */
bool Complex::operator==(const Complex &other) const
{
    return real() == other.real() and imag() == other.imag();
}

/** Comparison operator */
bool Complex::operator!=(const Complex &other) const
{
    return not operator==(other);
}

/** Assignment operator */
Complex &Complex::operator=(const Complex &other)
{
    val = other.val;
    return *this;
}

/** self-addition */
Complex &Complex::operator+=(const Complex &other)
{
    val += other.val;
    return *this;
}

/** self-subtraction */
Complex &Complex::operator-=(const Complex &other)
{
    val -= other.val;
    return *this;
}

/** self-multiplication */
Complex &Complex::operator*=(const Complex &other)
{
    val *= other.val;
    return *this;
}

/** self-division */
Complex &Complex::operator/=(const Complex &other)
{
    val /= other.val;
    return *this;
}

/** Return the negation of this number */
Complex Complex::operator-() const
{
    return Complex(-val);
}

/** Comparison operator */
bool Complex::operator==(double r) const
{
    return isReal() and SireMaths::areEqual(r, real());
}

/** Comparison operator */
bool Complex::operator!=(double r) const
{
    return not operator==(r);
}

/** Assignment operator. Note that, like setReal(), this sets only the real
    part, leaving the imaginary part unchanged */
Complex &Complex::operator=(double r)
{
    val.real(r);
    return *this;
}

/** self-addition */
Complex &Complex::operator+=(double r)
{
    val += r;
    return *this;
}

/** self-subtraction */
Complex &Complex::operator-=(double r)
{
    val -= r;
    return *this;
}

/** self-multiplication */
Complex &Complex::operator*=(double r)
{
    val *= r;
    return *this;
}

/** self-division */
Complex &Complex::operator/=(double r)
{
    val /= r;
    return *this;
}

/** Return the argument (phase angle) of this number, in (-pi, pi] */
double Complex::arg() const
{
    return std::arg(val);
}

/** Return the modulus |z| of this number */
double Complex::abs() const
{
    return std::abs(val);
}

/** Return the squared modulus |z|^2 of this number */
double Complex::abs2() const
{
    return std::norm(val);
}

/** Return log|z|. Factoring out the larger component keeps this accurate
    for |z| close to one, where log(abs(z)) would lose precision */
double Complex::logAbs() const
{
    double xabs = std::abs(val.real());
    double yabs = std::abs(val.imag());
    double max, u;

    if (xabs >= yabs)
    {
        max = xabs;
        u = yabs / xabs;
    }
    else
    {
        max = yabs;
        u = xabs / yabs;
    }

    return std::log(max) + 0.5 * std::log1p(u * u);
}

/** Return the complex conjugate, x - iy */
Complex Complex::conjugate() const
{
    return Complex(std::conj(val));
}

/** Return the reciprocal, 1/z */
Complex Complex::inverse() const
{
    return Complex(1.0 / val);
}

/** Return the negation, -z */
Complex Complex::negative() const
{
    return Complex(-val);
}

const char *Complex::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Complex>());
}

namespace SireMaths
{
    /** Addition */
    Complex operator+(const Complex &z0, const Complex &z1)
    {
        Complex ret(z0);
        ret += z1;
        return ret;
    }

    /** Subtraction */
    Complex operator-(const Complex &z0, const Complex &z1)
    {
        Complex ret(z0);
        ret -= z1;
        return ret;
    }

    /** Multiplication */
    Complex operator*(const Complex &z0, const Complex &z1)
    {
        Complex ret(z0);
        ret *= z1;
        return ret;
    }

    /** Division */
    Complex operator/(const Complex &z0, const Complex &z1)
    {
        Complex ret(z0);
        ret /= z1;
        return ret;
    }

    /** Addition */
    Complex operator+(const Complex &z, double x)
    {
        Complex ret(z);
        ret += x;
        return ret;
    }

    /** Subtraction */
    Complex operator-(const Complex &z, double x)
    {
        Complex ret(z);
        ret -= x;
        return ret;
    }

    /** Multiplication */
    Complex operator*(const Complex &z, double x)
    {
        Complex ret(z);
        ret *= x;
        return ret;
    }

    /** Division */
    Complex operator/(const Complex &z, double x)
    {
        Complex ret(z);
        ret /= x;
        return ret;
    }

    /** Addition */
    Complex operator+(double x, const Complex &z)
    {
        return z + x;
    }

    /** Subtraction. Note this is Complex(x) - z, *not* z - x - unlike + and *,
        subtraction does not commute. (This was previously a bug, returning
        z - x). */
    Complex operator-(double x, const Complex &z)
    {
        return Complex(x) - z;
    }

    /** Multiplication */
    Complex operator*(double x, const Complex &z)
    {
        return z * x;
    }

    /** Division */
    Complex operator/(double x, const Complex &z)
    {
        return Complex(x) / z;
    }

    /** Return the principal square root of z, which has a non-negative
        real part. The branch cut lies along the negative real axis */
    Complex sqrt(const Complex &z)
    {
        return Complex(std::sqrt(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the square root of the real number x, which is purely
        imaginary for negative x */
    Complex sqrt_real(double x)
    {
        if (x >= 0.0)
            return Complex(std::sqrt(x), 0.0);
        else
            return Complex(0.0, std::sqrt(-x));
    }

    /** Return z raised to the complex power a */
    Complex pow(const Complex &z, const Complex &a)
    {
        return Complex(std::pow(std::complex<double>(z.real(), z.imag()), std::complex<double>(a.real(), a.imag())));
    }

    /** Return z raised to the real power x */
    Complex pow(const Complex &z, double x)
    {
        return Complex(std::pow(std::complex<double>(z.real(), z.imag()), x));
    }

    /** Return z raised to the integer power n */
    Complex pow(const Complex &z, int n)
    {
        if (n >= 0)
        {
            switch (n)
            {
            case 0:
                return Complex(1);
            case 1:
                return z;
            case 2:
                return z * z;
            default:
                return pow(z, double(n));
            }
        }
        else
            return Complex(1) / pow(z, -n);
    }

    /** Return z raised to the rational power r */
    Complex pow(const Complex &z, const Rational &r)
    {
        if (r.denominator() == 1)
            return pow(z, r.numerator());
        else
            return pow(z, SireMaths::toDouble(r));
    }

    /** Return the real number x raised to the complex power z */
    Complex pow(double x, const Complex &z)
    {
        return pow(Complex(x), z);
    }

    /** Return e raised to the power z */
    Complex exp(const Complex &z)
    {
        return Complex(std::exp(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the principal natural logarithm of z. The branch cut lies
        along the negative real axis */
    Complex log(const Complex &z)
    {
        return Complex(std::log(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the principal base-10 logarithm of z */
    Complex log10(const Complex &z)
    {
        return Complex(std::log10(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the logarithm of z to the complex base b, log(z)/log(b) */
    Complex log_b(const Complex &z, const Complex &b)
    {
        return log(z) / log(b);
    }

    /** Return the sine of z */
    Complex sin(const Complex &z)
    {
        if (z.isReal())
            return std::sin(z.real());
        else
            return Complex(std::sin(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the cosine of z */
    Complex cos(const Complex &z)
    {
        if (z.isReal())
            return std::cos(z.real());
        else
            return Complex(std::cos(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the tangent of z */
    Complex tan(const Complex &z)
    {
        if (z.isReal())
            return std::tan(z.real());
        else
            return Complex(std::tan(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the secant of z, 1/cos(z) */
    Complex sec(const Complex &z)
    {
        return Complex(1.0) / cos(z);
    }

    /** Return the cosecant of z, 1/sin(z) */
    Complex csc(const Complex &z)
    {
        return Complex(1.0) / sin(z);
    }

    /** Return the cotangent of z, 1/tan(z) */
    Complex cot(const Complex &z)
    {
        return Complex(1.0) / tan(z);
    }

    /** Return the principal inverse sine of z. The branch cuts lie on the
        real axis outside [-1, 1] */
    Complex arcsin(const Complex &z)
    {
        return Complex(std::asin(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the inverse sine of the real number z, which is complex
        for |z| > 1 */
    Complex arcsin_real(double z)
    {
        return Complex(std::asin(std::complex<double>(z, 0.0)));
    }

    /** Return the principal inverse cosine of z. The branch cuts lie on the
        real axis outside [-1, 1] */
    Complex arccos(const Complex &z)
    {
        return Complex(std::acos(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the inverse cosine of the real number z, which is complex
        for |z| > 1 */
    Complex arccos_real(double z)
    {
        return Complex(std::acos(std::complex<double>(z, 0.0)));
    }

    /** Return the principal inverse tangent of z. The branch cuts lie on the
        imaginary axis outside [-i, i] */
    Complex arctan(const Complex &z)
    {
        return Complex(std::atan(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the inverse secant of z, arccos(1/z) */
    Complex arcsec(const Complex &z)
    {
        return arccos(Complex(1.0) / z);
    }

    /** Return the inverse secant of the real number z, arccos(1/z) */
    Complex arcsec_real(double z)
    {
        return arccos(Complex(1.0 / z));
    }

    /** Return the inverse cosecant of z, arcsin(1/z) */
    Complex arccsc(const Complex &z)
    {
        return arcsin(Complex(1.0) / z);
    }

    /** Return the inverse cosecant of the real number z, arcsin(1/z) */
    Complex arccsc_real(double z)
    {
        return arcsin(Complex(1.0 / z));
    }

    /** Return the inverse cotangent of z, arctan(1/z) */
    Complex arccot(const Complex &z)
    {
        return arctan(Complex(1.0) / z);
    }

    /** Return the hyperbolic sine of z */
    Complex sinh(const Complex &z)
    {
        return Complex(std::sinh(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the hyperbolic cosine of z */
    Complex cosh(const Complex &z)
    {
        return Complex(std::cosh(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the hyperbolic tangent of z */
    Complex tanh(const Complex &z)
    {
        return Complex(std::tanh(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the hyperbolic secant of z, 1/cosh(z) */
    Complex sech(const Complex &z)
    {
        return Complex(1.0) / cosh(z);
    }

    /** Return the hyperbolic cosecant of z, 1/sinh(z) */
    Complex csch(const Complex &z)
    {
        return Complex(1.0) / sinh(z);
    }

    /** Return the hyperbolic cotangent of z, 1/tanh(z) */
    Complex coth(const Complex &z)
    {
        return Complex(1.0) / tanh(z);
    }

    /** Return the principal inverse hyperbolic sine of z. The branch cuts lie
        on the imaginary axis outside [-i, i] */
    Complex arcsinh(const Complex &z)
    {
        return Complex(std::asinh(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the principal inverse hyperbolic cosine of z. The branch cut
        lies on the real axis below 1 */
    Complex arccosh(const Complex &z)
    {
        return Complex(std::acosh(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the inverse hyperbolic cosine of the real number z, which is
        complex for z < 1 */
    Complex arccosh_real(double z)
    {
        return Complex(std::acosh(std::complex<double>(z, 0.0)));
    }

    /** Return the principal inverse hyperbolic tangent of z. The branch cuts
        lie on the real axis outside [-1, 1] */
    Complex arctanh(const Complex &z)
    {
        return Complex(std::atanh(std::complex<double>(z.real(), z.imag())));
    }

    /** Return the inverse hyperbolic tangent of the real number z, which is
        complex for |z| > 1 */
    Complex arctanh_real(double z)
    {
        return Complex(std::atanh(std::complex<double>(z, 0.0)));
    }

    /** Return the inverse hyperbolic secant of z, arccosh(1/z) */
    Complex arcsech(const Complex &z)
    {
        return arccosh(Complex(1.0) / z);
    }

    /** Return the inverse hyperbolic cosecant of z, arcsinh(1/z) */
    Complex arccsch(const Complex &z)
    {
        return arcsinh(Complex(1.0) / z);
    }

    /** Return the inverse hyperbolic cotangent of z, arctanh(1/z) */
    Complex arccoth(const Complex &z)
    {
        return arctanh(Complex(1.0) / z);
    }
} // namespace SireMaths
