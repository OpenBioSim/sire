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
#include "SireMol/atomidxmapping.h"

#include "SireMM/mmdetail.h"

using namespace SireMol;
using namespace SireBase;

namespace SireSystem
{
    /**
     * @brief Merge function that combines multiple molecules into a single molecule.
     *
     * @param mols The AtomMapping object that contains the molecules to be merged.
     * @param properties The list of properties to be merged. If this is empty then
     *                   a default set of properties will be merged.
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
    Molecule merge(const AtomMapping &mols,
                   const QStringList &properties,
                   bool as_new_molecule,
                   bool allow_ring_breaking, bool allow_ring_size_change,
                   bool force, const PropertyMap &input_map)
    {
        if (not mols.isSingleMolecule())
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "You can only create a merged molecule from a mapping that "
                                                    "refers to a single molecule. You cannot use:\n%1")
                                                    .arg(mols.toString()));
        }

        PropertyMap map = input_map;

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

        map.set("as_new_molecule", BooleanProperty(as_new_molecule));
        map.set("allow_ring_breaking", BooleanProperty(allow_ring_breaking));

        // see which properties to merge
        QStringList props = properties;

        if (props.isEmpty())
        {
            props = QStringList({
                "angle",
                "ambertype",
                "atomtype",
                "bond",
                "charge",
                "connectivity",
                "coordinates",
                "dihedral",
                "element",
                "improper",
                "intrascale",
                "LJ",
                "mass",
            });
        }

        // get the forwards and backwards map
        auto forwards_map = mols;
        auto backwards_map = mols.swap();

        // the list of mapped atoms
        const auto mapped_atoms0 = mols.mappedAtoms0().toSingleMolecule();
        const auto mapped_atoms1 = mols.mappedAtoms1().toSingleMolecule();

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

        // get the MolEditor that can be used to set properties
        MolEditor editmol = mols.atoms0().toSingleMolecule().molecule().edit();

        // check and set the forcefields
        SireMM::MMDetail ffield0;
        SireMM::MMDetail ffield1;

        bool have_ffield0 = false;
        bool have_ffield1 = false;

        try
        {
            ffield0 = mapped_atoms0.data().property(map0["forcefield"]).asA<SireMM::MMDetail>();
            have_ffield0 = true;
        }
        catch (...)
        {
        }

        try
        {
            ffield1 = mapped_atoms1.data().property(map1["forcefield"]).asA<SireMM::MMDetail>();
            have_ffield1 = true;
        }
        catch (...)
        {
        }

        if (not(have_ffield0 and have_ffield1))
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "You must specify a forcefield for both the reference and "
                                                    "perturbed states in order to merge the molecules."),
                                                CODELOC);
        }
        else if (not have_ffield1)
        {
            ffield1 = ffield0;
        }
        else
        {
            ffield0 = ffield1;
        }

        if (not ffield0.isCompatibleWith(ffield1))
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "The forcefields for the reference and perturbed states are "
                                                    "incompatible. You cannot merge the molecules. The forcefields "
                                                    "are %1 and %2")
                                                    .arg(ffield0.toString())
                                                    .arg(ffield1.toString()),
                                                CODELOC);
        }

        // and a handle on the whole reference and perturbed molecule
        const auto mol0 = mols.atoms0().toSingleMolecule().molecule();
        const auto mol1 = mols.atoms1().toSingleMolecule().molecule();

        if (mol0.isEmpty())
            return mol0;

        // find the largest AtomNum in mol0
        AtomNum largest_atomnum;

        const auto &molinfo = mol0.info();

        for (int i = 0; i < molinfo.nAtoms(); ++i)
        {
            const auto num = molinfo.number(AtomIdx(i));

            if (largest_atomnum.isNull())
            {
                largest_atomnum = num;
            }
            else if (not num.isNull())
            {
                if (num.value() > largest_atomnum.value())
                    largest_atomnum = num;
            }
        }

        // use a property to track which atoms have been mapped -
        // a value of -1 means that this atom is not mapped
        editmol.setProperty("_mol0_index", AtomIntProperty(mol0.info(), -1));
        editmol.setProperty("_mol1_index", AtomIntProperty(mol1.info(), -1));

        // get an editable copy of the molecule to be changed
        MolStructureEditor mol(editmol);

        // a map from the residue index in the mapped state back to the
        // residue index in the merged molecule
        QHash<ResIdx, ResIdx> pert_to_merge_residx;

        // all of the residue indicies that we have seen
        QHash<ResIdx, CGIdx> residx_to_cgidx;

        // go through all of the common atoms and save their indicies
        // and set the atom and residue names
        for (int i = 0; i < nmapped; ++i)
        {
            const auto atom0 = mapped_atoms0(i);
            const auto atom1 = mapped_atoms1(i);

            auto atom = mol.atom(atom0.index());

            // save the index of this atom in both mol0 and mol1
            atom.setProperty<qint64>("_mol0_index", atom0.index().value());
            atom.setProperty<qint64>("_mol1_index", atom1.index().value());

            // save the perturbed state atom and residue names into new properties, so
            // that we can use these when extracting the end states
            atom.setAlternateName(atom1.name());

            ResIdx residx;

            try
            {
                residx = atom0.residue().index();
            }
            catch (...)
            {
            }

            if (not(residx.isNull() or residx_to_cgidx.contains(residx)))
            {
                // we haven't seen this residue before - assume that
                // all residues that are mapped from this residue
                // exist in the same equivalent residue in the
                // perturbed molecule (using the cutgroup of the
                // first atom in this residue)
                residx_to_cgidx.insert(residx, atom0.cutGroup().index());

                // first, get an editor for this residue
                // and save the alternate residue name for the mapped state
                auto res = mol.residue(residx);
                res.setAlternateName(atom1.residue().name());

                // now save the mapping from perturbed residue index
                // to merged residue index
                pert_to_merge_residx[atom1.residue().index()] = residx;
            }
        }

        // now go through the unmapped atoms of the reference molecule and
        // save their indicies
        const auto unmapped_atoms0 = mols.unmappedAtoms0().toSingleMolecule();

        for (int i = 0; i < unmapped_atoms0.count(); ++i)
        {
            const auto atom0 = unmapped_atoms0(i);

            auto atom = mol.atom(atom0.index());

            // unmapped atoms are called "Xxx"
            atom.setAlternateName("Xxx");
            atom.setProperty<qint64>("_mol0_index", atom0.index().value());
            atom.setProperty<qint64>("_mol1_index", -1);
        }

        // now go through the unmapped atoms of the perturbed molecule and
        // add them to the merged molecule, saving their indicies
        const auto unmapped_atoms1 = mols.unmappedAtoms1().toSingleMolecule();

        for (int i = 0; i < unmapped_atoms1.count(); ++i)
        {
            const auto atom1 = unmapped_atoms1(i);

            // we should have seen this residue before...
            auto residx = pert_to_merge_residx.value(atom1.residue().index());

            if (residx.isNull())
            {
                // we haven't seen this residue before, so we don't know
                // really where to add the atoms. The best thing to do
                // is add this to the last residue that we saw in the
                // molecule (the one with the highest index)
                if (residx_to_cgidx.isEmpty())
                {
                    throw SireError::program_bug(QObject::tr(
                                                     "We have not seen any residues before, so we don't know "
                                                     "where to add the atoms. This is a bug!"),
                                                 CODELOC);
                }

                auto residxs = residx_to_cgidx.keys();
                std::sort(residxs.begin(), residxs.end());

                residx = residxs.last();
            }

            auto cgidx = residx_to_cgidx.value(residx);

            if (cgidx.isNull())
            {
                throw SireError::program_bug(QObject::tr(
                                                 "We don't know the CutGroup for the residue, so we don't know "
                                                 "where to add the atoms. This is a bug!"),
                                             CODELOC);
            }

            auto res = mol.residue(residx);

            // add the atom - it has the name "Xxx" as it doesn't exist
            // in the reference state
            auto atom = res.add(AtomName("Xxx"));
            largest_atomnum = AtomNum(largest_atomnum.value() + 1);
            atom.renumber(largest_atomnum);

            // reparent this atom to the CutGroup for this residue
            atom.reparent(cgidx);

            // save the name in the perturbed state
            atom.setAlternateName(atom1.name());
            atom.setProperty<qint64>("_mol0_index", -1);
            atom.setProperty<qint64>("_mol1_index", atom1.index().value());
        }

        if (as_new_molecule)
        {
            mol.renumber();
        }

        editmol = mol.commit().edit();

        // now we have the merged molecule, we need to work out the mapping
        // of atoms from the reference to the perturbed state in this
        // merged molecule
        QList<AtomIdxMappingEntry> idx_entries;

        const auto &molinfo0 = editmol.info();
        const auto &molinfo1 = mol1.info();

        for (int i = 0; i < molinfo0.nAtoms(); ++i)
        {
            const auto atom = editmol.atom(AtomIdx(i));

            const auto index0 = atom.property<qint64>("_mol0_index");
            const auto index1 = atom.property<qint64>("_mol1_index");

            if (index0 == -1 and index1 == -1)
            {
                // this atom is not involved in the mapping
                continue;
            }

            const bool is_unmapped_in_reference = (index0 == -1);

            AtomIdx atomidx1(index1);

            if (index1 == -1)
            {
                // this atom is not in the perturbed state, so this
                // index should be set to null
                atomidx1 = AtomIdx();
            }

            idx_entries.append(AtomIdxMappingEntry(AtomIdx(i), atomidx1,
                                                   molinfo0, molinfo1,
                                                   is_unmapped_in_reference));
        }

        AtomIdxMapping entries(idx_entries);
        idx_entries.clear();

        // now go through all of the properties that we want to merge
        // and merge them using the AtomIdxMapping object - remove
        // the common property of these
        for (const auto &prop : props)
        {
            if (editmol.hasProperty(map0[prop]) and mol1.hasProperty(map1[prop]))
            {
                // we both have the property, so it should be mergeable
                const auto &prop0 = editmol.property(map0[prop]);
                const auto &prop1 = mol1.property(map1[prop]);

                if (prop0.what() != prop1.what())
                {
                    // the properties are not the same type
                    throw SireError::incompatible_error(QObject::tr(
                                                            "Cannot merge the molecule because the property %1 is not the "
                                                            "same type in both molecules. It is %2=%3 in the reference "
                                                            "and %4=%5 in the perturbed molecule.")
                                                            .arg(prop)
                                                            .arg(map0[prop].source())
                                                            .arg(prop0.what())
                                                            .arg(map1[prop].source())
                                                            .arg(prop1.what()),
                                                        CODELOC);
                }

                if (prop0.isA<MolViewProperty>())
                {
                    // they can be properly merged - get the value of the ghost parameter for this property
                    QString ghost_param;

                    if (map.specified("ghost_" + prop))
                    {
                        ghost_param = map["ghost_" + prop].source();
                    }
                    else if (prop == "atomtype" or prop == "ambertype")
                    {
                        ghost_param = "Xx";
                    }

                    auto merged = prop0.asA<MolViewProperty>().merge(prop1.asA<MolViewProperty>(),
                                                                     entries, ghost_param, map);

                    if (merged.count() != 2)
                        throw SireError::program_bug(QObject::tr(
                                                         "The merge of the property %1 did not return two properties. "
                                                         "This is a bug!")
                                                         .arg(prop),
                                                     CODELOC);

                    editmol.removeProperty(map0[prop]);
                    editmol.setProperty(map[prop + "0"].source(), merged[0]);
                    editmol.setProperty(map[prop + "1"].source(), merged[1]);
                }
                else
                {
                    // they are normal properties - they cannot be merged
                    // so just add them as the two end states
                    editmol.removeProperty(map0[prop]);
                    editmol.setProperty(map[prop + "0"].source(), prop0);
                    editmol.setProperty(map[prop + "1"].source(), prop1);
                }
            }
        }

        // add the reference and perturbed molecules as 'molecule0' and 'molecule1'
        editmol.setProperty(map["molecule0"].source(), mol0);
        editmol.setProperty(map["molecule1"].source(), mol1);

        // add the forcefields for the two molecules
        editmol.setProperty(map["forcefield0"].source(), ffield0);
        editmol.setProperty(map["forcefield1"].source(), ffield1);

        // remove any property called "parameters"
        if (editmol.hasProperty(map["parameters"].source()))
        {
            editmol.removeProperty(map["parameters"].source());
        }

        // set the flag that this is a perturbable molecule
        editmol.setProperty(map["is_perturbable"].source(), BooleanProperty(true));

        return editmol.commit();
    }
}
