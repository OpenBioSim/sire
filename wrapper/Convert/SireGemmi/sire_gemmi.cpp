
#include "sire_gemmi.h"

#include "gemmi/cif.hpp"
#include "gemmi/modify.hpp"
#include "gemmi/polyheur.hpp"

#include "SireIO/pdbx.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"

#include "SireError/errors.h"

namespace cif = gemmi::cif;

namespace SireGemmi
{
    void parse_waters(const gemmi::Entity &entity, SireMol::MoleculeGroup &mols,
                      const QHash<QString, int> &subchains,
                      const gemmi::Structure &structure)
    {
        for (const auto &subchain : entity.subchains)
        {
            auto model_id = subchains.value(QString::fromStdString(subchain), 0);

            for (const auto &residue : structure.models[model_id].get_subchain(subchain))
            {
                auto cg = SireMol::Molecule().edit().rename("WAT").add(SireMol::CGName("0"));
                auto res = cg.molecule().add(SireMol::ResNum(1));
                res.rename(SireMol::ResName("WAT"));

                for (const auto &atom : residue.atoms)
                {
                    auto atm = cg.add(SireMol::AtomName(QString::fromStdString(atom.name.c_str())));
                    atm.reparent(res.index());
                }

                mols.add(cg.molecule().commit());
            }
        }
    }

    SireSystem::System gemmi_to_sire(const gemmi::Structure &orig_structure,
                                     const SireBase::PropertyMap &map)
    {
        auto mols = SireMol::MoleculeGroup("all");

        gemmi::Structure structure = orig_structure;

        structure.merge_chain_parts();

        // this does not seem to exist in conda version?
        // gemmi::ensure_entities(structure);

        // create a dictionary of to locate which model contains which subchain
        QHash<QString, int> subchains;

        for (int i = 0; i < structure.models.size(); ++i)
        {
            for (auto subchain : structure.models[i].subchains())
            {
                auto id = QString::fromStdString(subchain.subchain_id());

                if (not subchains.contains(id))
                    subchains.insert(id, i);
            }
        }

        // this function doesn't appear to be in the gemmi API?
        // structure.standardize_crystal_frame();

        for (const auto &entity : structure.entities)
        {
            qDebug() << "entity" << QString::fromStdString(entity.name.c_str());

            switch (entity.entity_type)
            {
            case gemmi::EntityType::Polymer:
                qDebug() << "POLYMER";
                break;
            case gemmi::EntityType::NonPolymer:
                qDebug() << "NONPOLYMER";
                break;
            case gemmi::EntityType::Branched:
                qDebug() << "LIGAND";
                break;
            case gemmi::EntityType::Water:
                qDebug() << "WATER";
                parse_waters(entity, mols, subchains, structure);
                break;
            default:
                qDebug() << "UNKNOWN";
            }
        }

        SireSystem::System system;
        system.setName(QString::fromStdString(structure.name.c_str()));

        system.add(mols);

        return system;
    }

    gemmi::Structure sire_to_gemmi(const SireSystem::System &system,
                                   const SireBase::PropertyMap &map)
    {
        // TODO
        return gemmi::Structure();
    }

    SireSystem::System pdbx_reader_function(const QStringList &lines,
                                            const SireBase::PropertyMap &map)
    {
        // assemble all of the line into a single string
        auto input_string = lines.join("\n").toStdString();

        auto doc = cif::Document(cif::read_string(input_string));

        input_string.clear();

        // mmCIF files for deposition may have more than one block:
        // coordinates in the first block and restraints in the others.
        for (size_t i = 1; i < doc.blocks.size(); ++i)
            if (doc.blocks[i].has_tag("_atom_site.id"))
                throw SireError::unsupported(
                    QObject::tr("2+ blocks are ok if only the first one has coordinates. "
                                "_atom_site in block #%1 : %2")
                        .arg(i + 1)
                        .arg(QString::fromStdString(doc.source)),
                    CODELOC);

        auto structure = gemmi::make_structure_from_block(doc.blocks.at(0));

        return gemmi_to_sire(structure, map);
    }

    QStringList pdbx_writer_function(const SireSystem::System &system,
                                     const SireBase::PropertyMap &map)
    {
        // TODO
        return QStringList();
    }

    void register_pdbx_loader()
    {
        SireIO::PDBxReaderFunction reader_function(&pdbx_reader_function);
        SireIO::PDBxWriterFunction writer_function(&pdbx_writer_function);

        SireIO::register_pdbx_loader_functions(writer_function, reader_function);
    }
}
