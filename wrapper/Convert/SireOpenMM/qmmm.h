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

#ifndef SIREOPENMM_QMMM_H
#define SIREOPENMM_QMMM_H

#include "OpenMM.h"
#include "openmm/Force.h"

#include "sireglobal.h"

#include "SireBase/property.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    class QMMMEngine : public SireBase::Property, public OpenMM::Force
    {
    public:
        virtual ~QMMMEngine();

        //! Get the QM cutoff distance.
        virtual SireUnits::Dimension::Length getCutoff() const = 0;

        //! Set the QM cutoff distance.
        virtual void setCutoff(SireUnits::Dimension::Length cutoff) = 0;

        //! Get the indices of the atoms in the QM region.
        virtual QVector<int> getAtoms() const = 0;

        //! Set the list of atom indices for the QM region.
        virtual void setAtoms(QVector<int>) = 0;

        //! Get the atomic numbers of the atoms in the QM region.
        virtual QVector<int> getNumbers() const = 0;

        //! Set the list of atomic numbers for the QM region.
        virtual void setNumbers(QVector<int>) = 0;

        //! Get the atomic charges of all atoms in the system.
        virtual QVector<double> getCharges() const = 0;

        //! Set the atomic charges of all atoms in the system.
        virtual void setCharges(QVector<double>) = 0;

    protected:
        virtual OpenMM::ForceImpl *createImpl() const = 0;
    };

    typedef SireBase::PropPtr<SireOpenMM::QMMMEngine> QMMEnginePtr;
}

SIRE_END_HEADER

#endif
