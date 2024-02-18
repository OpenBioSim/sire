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

#include "SireSystem/mutate.h"

#include "SireMol/moleditor.h"
#include "SireMol/core.h"

#include "SireError/errors.h"

using namespace SireMol;
using namespace SireBase;

namespace SireSystem
{
    /**
     * @brief Mutate function that transforms a molecule into a new molecule.
     *
     * @param mols The AtomMapping object that contains the molecule to be mutated.
     * @param as_new_molecule Flag indicating whether the mutated molecule should be created as a new molecule.
     * @param map The PropertyMap object that contains additional properties for the mutated molecule.
     * @return The mutated molecule.
     */
    Molecule mutate(const AtomMapping &mols, bool as_new_molecule, const PropertyMap &map)
    {
        if (mols.atoms0().molecules().count() != 1 or mols.atoms1().molecules().count() != 1)
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "The passed mapping must contain just a single molecule to be mutated."),
                                                CODELOC);
        }

        auto mol = MolStructureEditor(mols.atoms0().molecules()[0]);

        // we will look up from mol1 to mol0
        const auto swapped_mols = mols.swap();

        // we will match residue by residue
        auto residues1 = mols.atoms1().molecules()[0].residues();

        Residue residue0;
        ResStructureEditor res;
        CGStructureEditor cg;

        for (int i = 0; i < residues1.count(); ++i)
        {
            const auto &residue1 = residues1(i);
            const auto atoms1 = residue1.atoms();

            QSet<AtomIdx> found_atoms;
            found_atoms.reserve(atoms1.count());

            for (int j = 0; j < atoms1.count(); ++j)
            {
                const auto &atom1 = atoms1(j);

                AtomStructureEditor atom;

                if (swapped_mols.contains(atom1))
                {
                    // update the existing atom
                    auto atom0 = swapped_mols[atom1];
                    found_atoms.insert(atom0.index());

                    if (residue0.isNull())
                    {
                        residue0 = atom0.residue();
                        qDebug() << "MUTATING" << residue0.toString() << "TO" << residue1.toString();

                        // rename the residue in the new molecule (we keep the same residue number)
                        res = mol.residue(residue0.index());
                        res.rename(residue1.name());
                        cg = res.atom(0).cutGroup();
                    }
                    else if (residue0 != atom0.residue())
                    {
                        throw SireError::incompatible_error(QObject::tr(
                                                                "The atoms in the mapping must belong to the same residue. "
                                                                "Currently mapping the atoms in %1 to %2, but atom %3 belongs to %4.")
                                                                .arg(residue0.toString())
                                                                .arg(residue1.toString())
                                                                .arg(atom0.toString())
                                                                .arg(atom0.residue().toString()),
                                                            CODELOC);
                    }

                    qDebug() << "UPDATING" << atom0.toString() << "TO" << atom1.toString();

                    // rename the atom in the new molecule (we keep the same atom number)
                    atom = mol.atom(atom0.index());
                    atom.rename(atom1.name());
                }
                else
                {
                    // add the new atom
                    qDebug() << "ADDING" << atom1.toString() << "TO" << residue0.toString();
                    atom = cg.add(atom1.name());
                    atom.reparent(residue0.index());

                    // calculate the coordinates of the new atom
                }
            } // for each atom in residue

            // find any atoms in the residue that weren't matched, and so need to be removed
            for (int j = 0; j < residue0.atoms().count(); ++j)
            {
                const auto &atom0 = residue0.atoms()(j);

                if (not found_atoms.contains(atom0.index()))
                {
                    qDebug() << "REMOVING" << atom0.toString();
                    mol.atom(atom0.index()).remove();
                }
            }

        } // for each residue

        if (as_new_molecule)
        {
            qDebug() << "RENUMBER";
            mol.renumber();
        }

        qDebug() << "COMMIT";
        return mol.commit();
    }

} // namespace SireSystem
