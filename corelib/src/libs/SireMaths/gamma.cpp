/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include <cmath>

#include "gamma.h"

#include <boost/math/special_functions/gamma.hpp>

namespace SireMaths
{

    /** Return the value of the Gamma function
        \Gamma(\alpha) = \int_0^{\infty} t^{\alpha-1} e^{-t} dt
    */
    double Gamma(double alpha)
    {
        return boost::math::tgamma(alpha);
    }

    /** Synonym for SireMaths::Gamma */
    double gamma(double alpha)
    {
        return Gamma(alpha);
    }

    /** Return the incomplete gamma function
        \Gamma(\alpha) = \int_x^{\infty} t^{\alpha-1} e^{-t} dt */
    double Gamma(double alpha, double x)
    {
        // boost::math::gamma_q is the regularized upper incomplete gamma
        // function Q(alpha,x) = Gamma(alpha,x)/Gamma(alpha), matching GSL's
        // gsl_sf_gamma_inc_Q
        return boost::math::gamma_q(alpha, x) * boost::math::tgamma(alpha);
    }

    /** Return the incomplete gamma function
        \Gamma(\alpha) = \int_0^{x} t^{\alpha-1} e^{-t} dt */
    double gamma(double alpha, double x)
    {
        // boost::math::gamma_p is the regularized lower incomplete gamma
        // function P(alpha,x) = gamma(alpha,x)/Gamma(alpha), matching GSL's
        // gsl_sf_gamma_inc_P
        return boost::math::gamma_p(alpha, x) * boost::math::tgamma(alpha);
    }

    /** Synonym for gamma(alpha,x) */
    double incomplete_gamma_lower(double alpha, double x)
    {
        return gamma(alpha, x);
    }

    /** Synonym for Gamma(alpha,x) */
    double incomplete_gamma_higher(double alpha, double x)
    {
        return Gamma(alpha, x);
    }

} // end of namespace SireMaths
