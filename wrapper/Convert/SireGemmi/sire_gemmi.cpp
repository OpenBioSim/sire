
#include "sire_gemmi.h"

namespace cif = gemmi::cif;

namespace SireGemi
{
    SireSystem::System gemmi_to_sire(const gemmi::Structure &structure,
                                     const SireBase::PropertyMap &map)
    {
        // TODO
        return SireSystem::System();
    }

    gemmi::Structure sire_to_gemmi(const SireSystem::System &system,
                                   const SireBase::PropertyMap &map)
    {
        // TODO
        return gemmi::Structure();
    }

    void register_pdbx_loader()
    {
        // TODO
    }
}

/*
void parseBlock(cif::Block &block, MoleculeGroup &mols)
{
    auto structure = gemmi::make_structure_from_block(block);

    gemmi::standardize_crystal_frame(structure);

    for (auto model : structure.models)
    {
        qDebug() << "model" << QString::fromStdString(model.name.c_str());

        for (auto chain : model.chains)
        {
            qDebug() << "chain" << QString::fromStdString(chain.name.c_str());

            qDebug() << "ligands" << chain.get_ligands().size();
            qDebug() << "polymers" << chain.get_polymer().size();
            qDebug() << "water" << chain.get_waters().size();

            for (auto residue : chain.residues)
            {
                qDebug() << "residue" << QString::fromStdString(residue.name.c_str())
                         << "is_water" << residue.is_water();

                for (auto atom : residue.atoms)
                {
                    qDebug() << "atom" << QString::fromStdString(atom.name.c_str());
                }
            }
        }
    }
}

void PDBx::parseLines(const PropertyMap &map)
{
    // assemble all of the line into a single string
    auto input_string = QStringList(lines().toList()).join("\n").toStdString();

    QMutexLocker lkr(&mutex);

    doc.reset(new cif::Document(cif::read_string(input_string)));

    input_string.clear();

    // construct each molecule
    MoleculeGroup mols("all");

    for (auto &block : doc->blocks)
    {
        parseBlock(block, mols);
    }

    this->setScore(nAtoms());
}
*/
