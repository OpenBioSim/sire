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

#ifndef SIREOPENMM_EMLECALLBACK_H
#define SIREOPENMM_EMLECALLBACK_H

#include "boost/python.hpp"
#include <boost/tuple/tuple.hpp>

#include <QVector>

#include "sireglobal.h"

namespace bp = boost::python;

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    // A callback wrapper class to allow use of electrostatic embedding of
    // machine learning potentials via emle-engine.
    class EMLECallback
    {
    public:
        //! Constructor
        /*! \param py_object
                A Python object that contains the callback function.

            \param callback
                The name of a callback method that take the following arguments:
                    - numbers_qm: A list of atomic numbers for the atoms in the ML region.
                    - charges_mm: A list of the MM charges.
                    - xyz_qm: A vector of positions for the atoms in the ML region.
                    - xyz_mm: A vector of positions for the atoms in the MM region.
         */
        EMLECallback(bp::object object, QString callback);

        //! Constructor
        /*! \param numbers_qm
                A vector of atomic numbers for the atoms in the ML region.

            \param charges_mm
                A vector of the charges on the MM atoms.

            \param xyz_qm
                A vector of positions for the atoms in the ML region.

            \param xyz_mm
                A vector of positions for the atoms in the MM region.

            \returns
                A flattened vector of forces for the QM and MM atoms.
         */
        boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>> call(
                QVector<int> numbers_qm,
                QVector<double> charges_mm,
                QVector<QVector<double>> xyz_qm,
                QVector<QVector<double>> xyz_mm
        );

    private:
        bp::object py_object;
        QString callback;
    };
}

SIRE_END_HEADER

#endif
