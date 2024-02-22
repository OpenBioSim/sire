/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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

#include "SireSystem/merge.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"

using namespace SireMol;
using namespace SireBase;

namespace SireSystem
{
    /**
     * @brief Merge function that combines multiple molecules into a single molecule.
     *
     * @param mols The AtomMapping object that contains the molecules to be merged.
     * @param as_new_molecule Flag indicating whether the merged molecule should be created as a new molecule.
     * @param allow_ring_breaking Whether to allow the opening/closing of rings during a merge.
     * @param allow_ring_size_change Whether to allow changes in ring size.
     * @param force       Whether to try to force the merge, even when the molecular
     *                   connectivity changes not as the result of a ring transformation.
     *                  This will likely lead to an unstable perturbation. This option
     *                 takes precedence over 'allow_ring_breaking' and
     *                'allow_ring_size_change'.
     * @param map The PropertyMap object that contains additional properties for the merged molecule.
     * @return The merged molecule.
     */
    Molecule merge(const AtomMapping &mols, bool as_new_molecule,
                   bool allow_ring_breaking, bool allow_ring_size_change,
                   bool force, const PropertyMap &map)
    {
        if (not mols.isSingleMolecule())
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "You can only create a merged molecule from a mapping that "
                                                    "refers to a single molecule. You cannot use:\n%1")
                                                    .arg(mols.toString()));
        }

        if (map.specified("as_new_molecule"))
        {
            as_new_molecule = map["as_new_molecule"].value().asABoolean();
        }

        if (map.specified("allow_ring_breaking"))
        {
            allow_ring_breaking = map["allow_ring_breaking"].value().asABoolean();
        }

        if (map.specified("allow_ring_size_change"))
        {
            allow_ring_size_change = map["allow_ring_size_change"].value().asABoolean();
        }

        if (map.specified("as_new_molecule"))
        {
            as_new_molecule = map["as_new_molecule"].value().asABoolean();
        }

        if (map.specified("force"))
        {
            force = map["force"].value().asABoolean();
        }

        if (force)
        {
            allow_ring_breaking = true;
            allow_ring_size_change = true;
        }

        // get the forwards and backwards map
        auto forwards_map = mols;
        auto backwards_map = mols.swap();

        // the list of mapped atoms
        const auto mapped_atoms0 = mols.mappedAtoms0();
        const auto mapped_atoms1 = mols.mappedAtoms1();

        if (mapped_atoms0.count() != mapped_atoms1.count())
        {
            throw SireError::program_bug(QObject::tr(
                                             "The number of atoms in the forward and backward mappings "
                                             "are not the same. This is a bug!."),
                                         CODELOC);
        }

        const int nmapped = mapped_atoms0.count();

        // get the merged maps for the reference and perturbed states
        auto map0 = map.merge(mols.propertyMap0());
        auto map1 = map.merge(mols.propertyMap1());

        // get an editable copy of the molecule to be changed
        MolStructureEditor mol(mols.atoms0().toSingleMolecule().molecule());

        // and a handle on the whole reference and perturbed molecule
        const auto mol0 = mols.atoms0().toSingleMolecule().molecule();
        const auto mol1 = mols.atoms1().toSingleMolecule().molecule();

        // copy the properties from the reference state to both states
        QStringList merged_properties = {"charge", "LJ", "atomtype", "intrascale",
                                         "coordinates", "mass", "element",
                                         "bond", "angle", "dihedral",
                                         "improper"};

        /*
                for (int i = 0; i < mol.nAtoms(); ++i)
                {
                    auto atom = mol.atom(AtomIdx(i));

                    for (const auto &prop : merged_properties)
                    {
                        if (mol0.hasProperty(map0[prop]))
                        {
                            const auto &val = mol0.property(map0[prop]);

                            atom.setProperty(map[prop + "0"].source(), val);
                            atom.setProperty(map[prop + "1"].source(), val);
                        }
                    }
                }

                QVector<AtomStructureEditor> matched_atoms;
                matched_atoms.reserve(nmapped);

                // now go through and update the values of properties for
                // the atoms that mutate
                ResStructureEditor last_res;
                bool changed_res = false;

                for (int i = 0; i < nmapped; ++i)
                {
                    const auto atom0 = mapped_atoms0[i];
                    const auto atom1 = mapped_atoms1[i];

                    matched_atoms.append(mol.atom(atom0.index()));
                    auto &atom = matched_atoms.last();

                    // save the perturbed state atom and residue names into new properties, so
                    // that we can use these when extracting the end states
                    atom.setAlternateName(atom1.name());

                    try
                    {
                        auto next_res = atom.residue();

                        if (next_res != last_res)
                        {
                            changed_res = true;
                            last_res = next_res;
                        }
                    }
                    catch (...)
                    {
                    }

                    if (changed_res)
                    {
                        last_res.setAlternateName(atom1.residue().name());
                    }

                    // check if we need to change residue
                }
        */

        if (as_new_molecule)
        {
            mol.renumber();
        }

        auto editmol = mol.commit().edit();

        // add the reference and perturbed molecules as 'molecule0' and 'molecule1'
        editmol.setProperty(map["molecule0"].source(), mol0);
        editmol.setProperty(map["molecule1"].source(), mol1);

        // set the flag that this is a perturbable molecule
        editmol.setProperty(map["is_perturbable"].source(), BooleanProperty(true));

        return editmol.commit();
    }
}
