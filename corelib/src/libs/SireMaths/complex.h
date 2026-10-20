/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIRECAS_COMPLEX_H
#define SIRECAS_COMPLEX_H

#include <QString>

#include <complex>

#include "SireMaths/maths.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
    class Complex;
}

class QDataStream;
SIREMATHS_EXPORT QDataStream &operator<<(QDataStream &, const SireMaths::Complex &);
SIREMATHS_EXPORT QDataStream &operator>>(QDataStream &, SireMaths::Complex &);

namespace SireMaths
{

    class Rational;

    /**
    This class represents a complex number to the same precision as 'double'.
    Internally, it is backed by a std::complex<double>.

    @author Christopher Woods
    */
    class SIREMATHS_EXPORT Complex
    {
    public:
        Complex(double r = 0.0, double i = 0.0);

        template <typename T>
        Complex(const std::complex<T> &stdcomplex);

        Complex(const Complex &other);

        ~Complex();

        static const char *typeName();

        const char *what() const
        {
            return Complex::typeName();
        }

        double real() const;
        double imag() const;

        template <typename T>
        operator std::complex<T>() const;

        bool isReal() const;
        bool isPurelyComplex() const;

        bool isZero() const;

        QString toString() const;

        static Complex rect(double x, double y);
        static Complex polar(double r, double theta);

        void setRectangular(double x, double y);
        void setPolar(double r, double theta);

        void setReal(double x);
        void setImag(double y);

        bool operator==(const Complex &other) const;
        bool operator!=(const Complex &other) const;

        template <typename T>
        bool operator==(const std::complex<T> &stdcomplex) const;

        template <typename T>
        bool operator!=(const std::complex<T> &stdcomplex) const;

        Complex &operator=(const Complex &other);

        template <typename T>
        Complex &operator=(const std::complex<T> &stdcomplex);

        Complex &operator+=(const Complex &other);
        Complex &operator-=(const Complex &other);
        Complex &operator*=(const Complex &other);
        Complex &operator/=(const Complex &other);

        Complex operator-() const;

        bool operator==(double r) const;
        bool operator!=(double r) const;

        Complex &operator=(double r);

        Complex &operator+=(double r);
        Complex &operator-=(double r);
        Complex &operator*=(double r);
        Complex &operator/=(double r);

        double arg() const;
        double abs() const;
        double abs2() const;

        double logAbs() const;

        Complex conjugate() const;

        Complex inverse() const;

        Complex negative() const;

    private:
        std::complex<double> val;
    };

    SIREMATHS_EXPORT Complex operator+(const Complex &z0, const Complex &z1);
    SIREMATHS_EXPORT Complex operator-(const Complex &z0, const Complex &z1);
    SIREMATHS_EXPORT Complex operator*(const Complex &z0, const Complex &z1);
    SIREMATHS_EXPORT Complex operator/(const Complex &z0, const Complex &z1);
    SIREMATHS_EXPORT Complex operator+(const Complex &z, double x);
    SIREMATHS_EXPORT Complex operator-(const Complex &z, double x);
    SIREMATHS_EXPORT Complex operator*(const Complex &z, double x);
    SIREMATHS_EXPORT Complex operator/(const Complex &z, double x);
    SIREMATHS_EXPORT Complex operator+(double x, const Complex &z);
    SIREMATHS_EXPORT Complex operator-(double x, const Complex &z);
    SIREMATHS_EXPORT Complex operator*(double x, const Complex &z);
    SIREMATHS_EXPORT Complex operator/(double x, const Complex &z);
    SIREMATHS_EXPORT Complex sqrt(const Complex &z);
    SIREMATHS_EXPORT Complex sqrt_real(double x);
    SIREMATHS_EXPORT Complex pow(const Complex &z, const Complex &a);
    SIREMATHS_EXPORT Complex pow(const Complex &z, double x);
    SIREMATHS_EXPORT Complex pow(const Complex &z, int n);
    SIREMATHS_EXPORT Complex pow(const Complex &z, const Rational &r);
    SIREMATHS_EXPORT Complex pow(double x, const Complex &z);
    SIREMATHS_EXPORT Complex exp(const Complex &z);
    SIREMATHS_EXPORT Complex log(const Complex &z);
    SIREMATHS_EXPORT Complex log10(const Complex &z);
    SIREMATHS_EXPORT Complex log_b(const Complex &z, const Complex &b);
    SIREMATHS_EXPORT Complex sin(const Complex &z);
    SIREMATHS_EXPORT Complex cos(const Complex &z);
    SIREMATHS_EXPORT Complex tan(const Complex &z);
    SIREMATHS_EXPORT Complex sec(const Complex &z);
    SIREMATHS_EXPORT Complex csc(const Complex &z);
    SIREMATHS_EXPORT Complex cot(const Complex &z);
    SIREMATHS_EXPORT Complex arcsin(const Complex &z);
    SIREMATHS_EXPORT Complex arcsin_real(double z);
    SIREMATHS_EXPORT Complex arccos(const Complex &z);
    SIREMATHS_EXPORT Complex arccos_real(double z);
    SIREMATHS_EXPORT Complex arctan(const Complex &z);
    SIREMATHS_EXPORT Complex arcsec(const Complex &z);
    SIREMATHS_EXPORT Complex arcsec_real(double z);
    SIREMATHS_EXPORT Complex arccsc(const Complex &z);
    SIREMATHS_EXPORT Complex arccsc_real(double z);
    SIREMATHS_EXPORT Complex arccot(const Complex &z);
    SIREMATHS_EXPORT Complex sinh(const Complex &z);
    SIREMATHS_EXPORT Complex cosh(const Complex &z);
    SIREMATHS_EXPORT Complex tanh(const Complex &z);
    SIREMATHS_EXPORT Complex sech(const Complex &z);
    SIREMATHS_EXPORT Complex csch(const Complex &z);
    SIREMATHS_EXPORT Complex coth(const Complex &z);
    SIREMATHS_EXPORT Complex arcsinh(const Complex &z);
    SIREMATHS_EXPORT Complex arccosh(const Complex &z);
    SIREMATHS_EXPORT Complex arccosh_real(double z);
    SIREMATHS_EXPORT Complex arctanh(const Complex &z);
    SIREMATHS_EXPORT Complex arctanh_real(double z);
    SIREMATHS_EXPORT Complex arcsech(const Complex &z);
    SIREMATHS_EXPORT Complex arccsch(const Complex &z);
    SIREMATHS_EXPORT Complex arccoth(const Complex &z);

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

    /** Construct from a std::complex */
    template <typename T>
    Complex::Complex(const std::complex<T> &stdcomplex)
        : val(stdcomplex.real(), stdcomplex.imag())
    {
    }

    /** Assignment from a std::complex */
    template <typename T>
    Complex &Complex::operator=(const std::complex<T> &stdcomplex)
    {
        setReal(stdcomplex.real());
        setImag(stdcomplex.imag());
        return *this;
    }

    /** Implicit conversion to a std::complex */
    template <typename T>
    Complex::operator std::complex<T>() const
    {
        return std::complex<T>(real(), imag());
    }

    /** Comparison to std::complex */
    template <typename T>
    bool Complex::operator==(const std::complex<T> &stdcomplex) const
    {
        return SireMaths::areEqual(real(), stdcomplex.real()) and SireMaths::areEqual(imag(), stdcomplex.imag());
    }

    /** Comparison to std::complex */
    template <typename T>
    bool Complex::operator!=(const std::complex<T> &stdcomplex) const
    {
        return not operator==(stdcomplex);
    }

    /** Comparison with std::complex */
    template <typename T>
    SIRE_INLINE_TEMPLATE bool operator==(const std::complex<T> &stdcomplex, const Complex &complex)
    {
        return complex == stdcomplex;
    }

    /** Comparison with std::complex */
    template <typename T>
    SIRE_INLINE_TEMPLATE bool operator!=(const std::complex<T> &stdcomplex, const Complex &complex)
    {
        return complex != stdcomplex;
    }

    /** This is the std sqrt function. This helps resolve namespace issues... */
    SIRE_ALWAYS_INLINE double sqrt(double x)
    {
        return std::sqrt(x);
    }

#endif // SIRE_SKIP_INLINE_FUNCTIONS

} // namespace SireMaths

Q_DECLARE_METATYPE(SireMaths::Complex)
Q_DECLARE_TYPEINFO(SireMaths::Complex, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS(SireMaths::Complex)

SIRE_END_HEADER

#endif
