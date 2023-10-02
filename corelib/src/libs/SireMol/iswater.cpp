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

#include "iswater.h"

#include "atomelements.h"

#include "SireBase/parallel.h"

using namespace SireMol;
using namespace SireBase;

namespace SireMol
{
    bool _is_water(const AtomElements &elements)
    {
        // Counters for the number of hydrogens, oxygens, and protons in the molecule.
        int num_hydrogen = 0;
        int num_oxygen = 0;
        int num_protons = 0;

        // Loop over all cut-groups associated with the elements.
        for (int i = 0; i < elements.nCutGroups(); ++i)
        {
            // Create the cut-group index.
            CGIdx cg(i);

            // Extract the data for this cut-group.
            auto data = elements.constData(cg);

            // Loop over all atoms in this cut-group.
            for (int j = 0; j < elements.nAtoms(cg); ++j)
            {
                // Get the element.
                const auto element = data[j];

                // Update the number of protons.
                num_protons += element.nProtons();

                // Hydrogen.
                if (element.nProtons() == 1)
                    num_hydrogen++;
                // Oxygen.
                else if (element.nProtons() == 8)
                    num_oxygen++;

                // Not a water molecule, abort!
                if (num_oxygen > 1 or num_hydrogen > 2 or num_protons > 10)
                {
                    return false;
                }
            }
        }

        return (num_oxygen == 1 and num_hydrogen == 2 and num_protons == 10);
    }

    SIREMOL_EXPORT
    bool is_water(const MoleculeView &molecule, const PropertyMap &map)
    {
        // Convert to a molecule.
        const auto &moldata = molecule.data();

        if (moldata.info().nAtoms() > 6)
            // we won't check molecules that are full of dummy atoms
            return false;

        const auto element_property = map["element"];

        // Skip if there is no element property.
        if (not moldata.hasProperty(element_property))
            return false;

        return _is_water(moldata.property(element_property).asA<AtomElements>());
    }

    SIREMOL_EXPORT
    QVector<bool> is_water(const SelectorMol &molecules, const PropertyMap &map)
    {
        QVector<bool> result;

        if (molecules.isEmpty())
            return result;

        const int nmols = molecules.count();

        result = QVector<bool>(nmols, false);
        auto result_data = result.data();

        const auto element_property = map["element"];

        if (SireBase::should_run_in_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    const auto &moldata = molecules[i].data();

                    if (moldata.info().nAtoms() <= 6 and moldata.hasProperty(element_property))
                    {
                        result_data[i] = _is_water(moldata.property(element_property).asA<AtomElements>());
                    }
            } });
        }
        else
        {
            for (int i = 0; i < molecules.count(); ++i)
            {
                const auto &moldata = molecules[i].data();

                if (moldata.info().nAtoms() <= 6 and moldata.hasProperty(element_property))
                {
                    result_data[i] = _is_water(moldata.property(element_property).asA<AtomElements>());
                }
            }
        }

        return result;
    }

    SIREMOL_EXPORT
    QVector<bool> is_water(const SelectResult &molecules, const PropertyMap &map)
    {
        QVector<bool> result;

        if (molecules.isEmpty())
            return result;

        const auto &mols = molecules.toList();

        const int nmols = mols.count();

        result = QVector<bool>(nmols, false);
        auto result_data = result.data();

        const auto element_property = map["element"];

        if (SireBase::should_run_in_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    const auto &moldata = mols.at(i).read().data();

                    if (moldata.info().nAtoms() <= 6 and moldata.hasProperty(element_property))
                    {
                        result_data[i] = _is_water(moldata.property(element_property).asA<AtomElements>());
                    }
            } });
        }
        else
        {
            for (int i = 0; i < nmols; ++i)
            {
                const auto &moldata = mols.at(i).read().data();

                if (moldata.info().nAtoms() <= 6 and moldata.hasProperty(element_property))
                {
                    result_data[i] = _is_water(moldata.property(element_property).asA<AtomElements>());
                }
            }
        }

        return result;
    }

}
