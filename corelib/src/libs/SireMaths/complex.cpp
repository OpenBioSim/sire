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

/** This function uses the rectangular cartesian components (x,y) to
    return the complex number z = x + i y. */
Complex Complex::rect(double x, double y)
{
    return Complex(x, y);
}

/** This function returns the complex number
    z = r \exp(i \theta) = r (\cos(\theta) + i \sin(\theta))
    from the polar representation (r,theta). */
Complex Complex::polar(double r, double theta)
{
    return Complex(std::polar(r, theta));
}

/** This function sets the rectangular cartesian components (x,y) to
    the complex number z = x + i y. */
void Complex::setRectangular(double x, double y)
{
    val = std::complex<double>(x, y);
}

/** This function sets the complex number to
    z = r \exp(i \theta) = r (\cos(\theta) + i \sin(\theta))
    from the polar representation (r,theta). */
void Complex::setPolar(double r, double theta)
{
    val = std::polar(r, theta);
}

/** This function sets the real part of the complex number, leaving the
    imaginary part unchanged */
void Complex::setReal(double x)
{
    val.real(x);
}

/** This function sets the imaginary part of the complex number, leaving
    the real part unchanged */
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

/** This function returns the negative of the complex
    number z, -z = (-x) + i(-y). */
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

/** This function returns the argument of the complex number z,
    \arg(z), where -\pi < \arg(z) <= \pi. */
double Complex::arg() const
{
    return std::arg(val);
}

/** This function returns the magnitude of the complex number z, |z|. */
double Complex::abs() const
{
    return std::abs(val);
}

/** This function returns the squared magnitude of the complex number z, |z|^2. */
double Complex::abs2() const
{
    return std::norm(val);
}

/** This function returns the natural logarithm of the magnitude of the
    complex number z, \log|z|. It allows an accurate evaluation of \log|z|
    when |z| is close to one. The direct evaluation of log(abs(z))
    would lead to a loss of precision in this case, so this uses the same
    factor-out-the-largest-component technique as the GSL routine this used
    to call. */
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

/** This function returns the complex conjugate of the complex
    number z, z^* = x - i y. */
Complex Complex::conjugate() const
{
    return Complex(std::conj(val));
}

/** This function returns the inverse, or reciprocal, of the
    complex number z, 1/z = (x - i y)/(x^2 + y^2). */
Complex Complex::inverse() const
{
    return Complex(1.0 / val);
}

/** This function returns the negative of the complex
    number z, -z = (-x) + i(-y). */
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

    /** This function returns the square root of the complex number z, \sqrt z.
        The branch cut is the negative real axis. The result always lies in the
        right half of the complex plane. */
    Complex sqrt(const Complex &z)
    {
        return Complex(std::sqrt(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex square root of the real number x,
        where x may be negative. */
    Complex sqrt_real(double x)
    {
        if (x >= 0.0)
            return Complex(std::sqrt(x), 0.0);
        else
            return Complex(0.0, std::sqrt(-x));
    }

    /** The function returns the complex number z raised to the complex power a,
        z^a. This is computed as \exp(\log(z)*a) using complex logarithms and complex
        exponentials. */
    Complex pow(const Complex &z, const Complex &a)
    {
        return Complex(std::pow(std::complex<double>(z.real(), z.imag()), std::complex<double>(a.real(), a.imag())));
    }

    /** This function returns the complex number z raised to the real power x, z^x. */
    Complex pow(const Complex &z, double x)
    {
        return Complex(std::pow(std::complex<double>(z.real(), z.imag()), x));
    }

    /** This function returns the complex number z raised to an integer power n, z^n */
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

    /** This function returns the complex number z raised to a rational power r, z^r */
    Complex pow(const Complex &z, const Rational &r)
    {
        if (r.denominator() == 1)
            return pow(z, r.numerator());
        else
            return pow(z, SireMaths::toDouble(r));
    }

    /** This function returns the real number x raised to the complex power z, x^z. */
    Complex pow(double x, const Complex &z)
    {
        return pow(Complex(x), z);
    }

    /** This function returns the complex exponential of the complex number z,
        \exp(z). */
    Complex exp(const Complex &z)
    {
        return Complex(std::exp(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex natural logarithm (base e) of the complex
        number z, \log(z). The branch cut is the negative real axis. */
    Complex log(const Complex &z)
    {
        return Complex(std::log(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex base-10 logarithm of the complex number z,
        \log_10 (z). */
    Complex log10(const Complex &z)
    {
        return Complex(std::log10(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex base-b logarithm of the complex number z,
       \log_b(z). This quantity is computed as the ratio \log(z)/\log(b). */
    Complex log_b(const Complex &z, const Complex &b)
    {
        return log(z) / log(b);
    }

    /** This function returns the complex sine of the complex number
        z, \sin(z) = (\exp(iz) - \exp(-iz))/(2i). */
    Complex sin(const Complex &z)
    {
        if (z.isReal())
            return std::sin(z.real());
        else
            return Complex(std::sin(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex cosine of the complex number
        z, \cos(z) = (\exp(iz) + \exp(-iz))/2. */
    Complex cos(const Complex &z)
    {
        if (z.isReal())
            return std::cos(z.real());
        else
            return Complex(std::cos(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex tangent of the complex number z,
        \tan(z) = \sin(z)/\cos(z). */
    Complex tan(const Complex &z)
    {
        if (z.isReal())
            return std::tan(z.real());
        else
            return Complex(std::tan(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex secant of the complex number z,
        \sec(z) = 1/\cos(z). */
    Complex sec(const Complex &z)
    {
        return Complex(1.0) / cos(z);
    }

    /** This function returns the complex cosecant of the complex number z,
        \csc(z) = 1/\sin(z). */
    Complex csc(const Complex &z)
    {
        return Complex(1.0) / sin(z);
    }

    /** This function returns the complex cotangent of the complex number z,
        \cot(z) = 1/\tan(z). */
    Complex cot(const Complex &z)
    {
        return Complex(1.0) / tan(z);
    }

    /** This function returns the complex arcsine of the complex number z,
        \arcsin(z). The branch cuts are on the real axis, less than -1 and greater
        than 1. */
    Complex arcsin(const Complex &z)
    {
        return Complex(std::asin(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex arcsine of the real number z,
        \arcsin(z). For z between -1 and 1, the function returns a real
        value in the range [-\pi/2,\pi/2]. For z less than -1 the result
        has a real part of -\pi/2 and a positive imaginary part. For z greater
        than 1 the result has a real part of \pi/2 and a negative imaginary part. */
    Complex arcsin_real(double z)
    {
        return Complex(std::asin(std::complex<double>(z, 0.0)));
    }

    /** This function returns the complex arccosine of the complex number z,
        \arccos(z). The branch cuts are on the real axis, less than -1 and
         greater than 1. */
    Complex arccos(const Complex &z)
    {
        return Complex(std::acos(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex arccosine of the real number z,
        \arccos(z). For z between -1 and 1, the function returns a real value in
        the range [0,\pi]. For z less than -1 the result has a real part of \pi
        and a negative imaginary part. For z greater than 1 the result is purely
        imaginary and positive. */
    Complex arccos_real(double z)
    {
        return Complex(std::acos(std::complex<double>(z, 0.0)));
    }

    /** This function returns the complex arctangent of the complex number z,
        \arctan(z). The branch cuts are on the imaginary axis, below -i and above i. */
    Complex arctan(const Complex &z)
    {
        return Complex(std::atan(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex arcsecant of the complex number z,
        \arcsec(z) = \arccos(1/z). */
    Complex arcsec(const Complex &z)
    {
        return arccos(Complex(1.0) / z);
    }

    /** This function returns the complex arcsecant of the real number z,
        \arcsec(z) = \arccos(1/z). */
    Complex arcsec_real(double z)
    {
        return arccos(Complex(1.0 / z));
    }

    /** This function returns the complex arccosecant of the complex number z,
        \arccsc(z) = \arcsin(1/z). */
    Complex arccsc(const Complex &z)
    {
        return arcsin(Complex(1.0) / z);
    }

    /** This function returns the complex arccosecant of the real number z,
        \arccsc(z) = \arcsin(1/z). */
    Complex arccsc_real(double z)
    {
        return arcsin(Complex(1.0 / z));
    }

    /** This function returns the complex arccotangent of the complex number z,
        \arccot(z) = \arctan(1/z). */
    Complex arccot(const Complex &z)
    {
        return arctan(Complex(1.0) / z);
    }

    /** This function returns the complex hyperbolic sine of the complex number z,
        \sinh(z) = (\exp(z) - \exp(-z))/2. */
    Complex sinh(const Complex &z)
    {
        return Complex(std::sinh(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex hyperbolic cosine of the complex number z,
        \cosh(z) = (\exp(z) + \exp(-z))/2.  */
    Complex cosh(const Complex &z)
    {
        return Complex(std::cosh(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex hyperbolic tangent of the complex number z,
        \tanh(z) = \sinh(z)/\cosh(z). */
    Complex tanh(const Complex &z)
    {
        return Complex(std::tanh(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex hyperbolic secant of the complex number z,
        \sech(z) = 1/\cosh(z). */
    Complex sech(const Complex &z)
    {
        return Complex(1.0) / cosh(z);
    }

    /** This function returns the complex hyperbolic cosecant of the complex number z,
        \csch(z) = 1/\sinh(z). */
    Complex csch(const Complex &z)
    {
        return Complex(1.0) / sinh(z);
    }

    /** This function returns the complex hyperbolic cotangent of the complex number z,
        \coth(z) = 1/\tanh(z). */
    Complex coth(const Complex &z)
    {
        return Complex(1.0) / tanh(z);
    }

    /** This function returns the complex hyperbolic arcsine of the complex number z,
        \arcsinh(z). The branch cuts are on the imaginary axis, below -i and above i. */
    Complex arcsinh(const Complex &z)
    {
        return Complex(std::asinh(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex hyperbolic arccosine of the complex number z,
        \arccosh(z). The branch cut is on the real axis, less than 1. */
    Complex arccosh(const Complex &z)
    {
        return Complex(std::acosh(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex hyperbolic arccosine of the real number z,
        \arccosh(z). */
    Complex arccosh_real(double z)
    {
        return Complex(std::acosh(std::complex<double>(z, 0.0)));
    }

    /** This function returns the complex hyperbolic arctangent of the complex number z,
        \arctanh(z). The branch cuts are on the real axis, less than -1 and greater than 1. */
    Complex arctanh(const Complex &z)
    {
        return Complex(std::atanh(std::complex<double>(z.real(), z.imag())));
    }

    /** This function returns the complex hyperbolic arctangent of the real number z,
        \arctanh(z). */
    Complex arctanh_real(double z)
    {
        return Complex(std::atanh(std::complex<double>(z, 0.0)));
    }

    /** This function returns the complex hyperbolic arcsecant of the complex number z,
       \arcsech(z) = \arccosh(1/z). */
    Complex arcsech(const Complex &z)
    {
        return arccosh(Complex(1.0) / z);
    }

    /** This function returns the complex hyperbolic arccosecant of the complex number z,
        \arccsch(z) = \arcsinh(1/z). */
    Complex arccsch(const Complex &z)
    {
        return arcsinh(Complex(1.0) / z);
    }

    /** This function returns the complex hyperbolic arccotangent of the complex number z,
        \arccoth(z) = \arctanh(1/z). */
    Complex arccoth(const Complex &z)
    {
        return arctanh(Complex(1.0) / z);
    }
} // namespace SireMaths
