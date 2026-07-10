/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMATHS_CONSTANTS_H
#define SIREMATHS_CONSTANTS_H

#include <limits>

#include "sireglobal.h"

SIRE_BEGIN_HEADER

/** These used to be exposed via GSL's gsl_math.h #define constants; now
    hardcoded directly (to the same precision GSL/glibc's math.h use) so
    this doesn't depend on GSL, or on non-standard math.h extensions like
    M_PI_2/M_SQRT2 that aren't guaranteed portable (e.g. not defined by
    MSVC without _USE_MATH_DEFINES) */
namespace SireMaths
{

    const double e = 2.71828182845904523536; /* e */

    const double log2_e = 1.44269504088896340736; /* log_2 (e) */

    const double log10_e = 0.43429448190325182765; /* log_10 (e) */

    const double sqrt_two = 1.41421356237309504880; /* sqrt(2) */

    const double sqrt_half = 0.70710678118654752440; /* sqrt(1/2) */

    const double sqrt_three = 1.73205080756887729353; /* sqrt(3) */

    const double pi = 3.141592653589793238462643383279; /* pi */

    const double two_pi = double(2) * pi; /* 2 * pi */

    const double pi_over_two = 1.57079632679489661923; /* pi/2 */

    const double pi_4 = 0.78539816339744830962; /* pi/4 */

    const double sqrtpi = 1.77245385090551602730; /* sqrt(pi) */

    const double two_sqrtpi = 1.12837916709551257390; /* 2/sqrt(pi) */

    const double one_over_pi = 0.31830988618379067154; /* 1/pi */

    const double two_over_pi = 0.63661977236758134308; /* 2/pi */

    const double ln_ten = 2.30258509299404568402; /* ln(10) */

    const double ln_two = 0.69314718055994530942; /* ln(2) */

    const double ln_pi = 1.14472988584940017414; /* ln(pi) */

    const double euler = 0.57721566490153286061; /* Euler constant */

    // now expose some sizes...
    /** A small number */
    const double small = 0.001;

    /** A tiny number */
    const double tiny = 1.0e-8;

    /** The smallest possible number */
    const double smallest = std::numeric_limits<double>::min();

    /** A large number */
    const double large = 1000.0;

    /** A huge number */
    const double huge = 1.0e8;

    /** The largest possible number */
    const double largest = std::numeric_limits<double>::max();

    /** The smallest possible int */
    const int smallest_int = std::numeric_limits<int>::min();

    /** The largest possible int */
    const int largest_int = std::numeric_limits<int>::max();

} // namespace SireMaths

SIRE_END_HEADER

#endif
