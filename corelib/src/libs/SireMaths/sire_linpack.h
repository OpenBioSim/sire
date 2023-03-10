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

#ifndef SIREMATHS_SIRE_LINPACK_H
#define SIREMATHS_SIRE_LINPACK_H

#include "sireglobal.h"

#include <QVector>

#include <utility>

SIRE_BEGIN_HEADER

namespace SireMaths
{

    class NMatrix;
    class NVector;

    std::pair<NMatrix, QVector<int>> dgeco(const NMatrix &A);

    NMatrix dgedi_inverse(const NMatrix &A, const QVector<int> &IPVT);

    double dgedi_determinant(const NMatrix &A, const QVector<int> &IPVT);

    std::pair<double, NMatrix> dgedi(const NMatrix &A, const QVector<int> &IPVT);

} // namespace SireMaths

SIRE_END_HEADER

#endif
