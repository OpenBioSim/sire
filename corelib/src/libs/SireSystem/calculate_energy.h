/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIRESYSTEM_CALCULATE_ENERGY_H
#define SIRESYSTEM_CALCULATE_ENERGY_H

#include "SireBase/propertymap.h"
#include "SireFF/forcefields.h"
#include "SireMol/molecules.h"
#include "SireMol/moleculeview.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireSystem
{
  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(SireFF::ForceFields &ffields);

  SIREMM_EXPORT SireFF::ForceFields create_forcefield(const SireMol::MoleculeView &mol, const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireFF::ForceFields create_forcefield(const SireMol::Molecules &mols, const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireFF::ForceFields create_forcefield(const SireMol::MoleculeView &mol0,
                                                      const SireMol::MoleculeView &mol1,
                                                      const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireFF::ForceFields create_forcefield(const SireMol::MoleculeView &mol0, const SireMol::Molecules &mols1,
                                                      const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireFF::ForceFields create_forcefield(const SireMol::Molecules &mols0, const SireMol::Molecules &mols1,
                                                      const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::MoleculeView &mol);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::MoleculeView &mol,
                                                                   const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::Molecules &mols);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::Molecules &mols,
                                                                   const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                                                   const SireMol::MoleculeView &mol1);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                                                   const SireMol::Molecules &mols1);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::Molecules &mols0,
                                                                   const SireMol::Molecules &mols1);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                                                   const SireMol::MoleculeView &mol1,
                                                                   const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0,
                                                                   const SireMol::Molecules &mols1,
                                                                   const SireBase::PropertyMap &map);

  SIREMM_EXPORT SireUnits::Dimension::GeneralUnit calculate_energy(const SireMol::Molecules &mols0,
                                                                   const SireMol::Molecules &mols1,
                                                                   const SireBase::PropertyMap &map);

  SIREMM_EXPORT QVector<SireUnits::Dimension::GeneralUnit> calculate_trajectory_energy(const SireFF::ForceFields &ff,
                                                                                       const QList<qint64> &frames,
                                                                                       const SireBase::PropertyMap &map);

  SIREMM_EXPORT QVector<QVector<SireUnits::Dimension::GeneralUnit>> calculate_trajectory_energies(
      const QVector<SireFF::ForceFields> &ff, const QList<qint64> &frames, const SireBase::PropertyMap &map);

} // namespace SireMM

SIRE_EXPOSE_FUNCTION(SireSystem::create_forcefield)
SIRE_EXPOSE_FUNCTION(SireSystem::calculate_energy)
SIRE_EXPOSE_FUNCTION(SireSystem::calculate_trajectory_energy)
SIRE_EXPOSE_FUNCTION(SireSystem::calculate_trajectory_energies)

SIRE_END_HEADER

#endif
