
#include "sire_gemmi.h"

#include "gemmi/cif.hpp"
#include "gemmi/modify.hpp"
#include "gemmi/polyheur.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/to_mmcif.hpp"

#include "SireIO/pdbx.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"
#include "SireMol/element.h"

#include "SireMol/atomproperty.hpp"
#include "SireMol/atomelements.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "SireMol/connectivity.h"
#include "SireMol/bondhunter.h"

#include "SireError/errors.h"

#include <string>
#include <strstream>

namespace cif = gemmi::cif;

namespace SireGemmi
{
    template <class T, class V>
    void set_prop(T &atom, const QString &key, const V &value,
                  const SireBase::PropertyMap &map)
    {
        const auto name = map[key];

        if (name.hasSource())
            atom.setProperty(name.source(), value);
    }

    SireMaths::Vector from_vec(const gemmi::Position &pos)
    {
        return SireMaths::Vector(pos.x, pos.y, pos.z);
    }

    void populate_atom(SireMol::AtomStructureEditor &atm,
                       const gemmi::Atom &atom,
                       const QString &is_hetatm,
                       const SireBase::PropertyMap &map)
    {
        // the charge is a signed-char
        double chg = int(atom.charge);

        set_prop(atm, "element", SireMol::Element(atom.element.atomic_number()), map);
        set_prop(atm, "formal_charge", SireUnits::Dimension::Charge(chg), map);
        set_prop(atm, "occupancy", double(atom.occ), map);
        set_prop(atm, "beta_factor", double(atom.b_iso), map);
        set_prop(atm, "is_het", is_hetatm, map);
        set_prop(atm, "alt_loc", QString(atom.altloc), map);
        set_prop(atm, "coordinates", from_vec(atom.pos), map);
    }

    void parse_polymer(const gemmi::Entity &entity, SireMol::MoleculeGroup &mols,
                       const QHash<QString, int> &subchains,
                       const gemmi::Structure &structure,
                       QHash<int, SireMol::MolNum> &serial_to_molnum,
                       const QString &molname,
                       const SireBase::PropertyMap &map)
    {
        // parse this into a single polymer comprised of multiple chains
        auto mol = SireMol::MolStructureEditor(SireMol::Molecule().edit());
        mol.renumber();
        mol.rename(SireMol::MolName(molname));

        auto cg0 = mol.add(SireMol::CGName("0"));

        QHash<QString, SireMol::SegStructureEditor> segments;

        const auto molnum = cg0.molecule().number();

        int cg_num = 0;

        for (const auto &subchain : entity.subchains)
        {
            auto model_id = subchains.value(QString::fromStdString(subchain), 0);

            auto chain = mol.add(SireMol::ChainName(QString::fromStdString(subchain.c_str())));

            for (const auto &residue : structure.models[model_id].get_subchain(subchain))
            {
                auto cg = cg0;

                if (cg_num > 0)
                {
                    cg = mol.add(SireMol::CGName(QString::number(cg_num)));
                }

                cg_num += 1;

                int resnum = cg_num;

                if (residue.seqid.num.has_value())
                    resnum = residue.seqid.num.value;

                auto res = mol.add(SireMol::ResName(QString::fromStdString(residue.name)));
                res.reparent(chain.index());
                res.renumber(SireMol::ResNum(resnum));

                if (residue.seqid.has_icode())
                    set_prop(res, "insert_code", QString(residue.seqid.icode), map);

                SireMol::SegStructureEditor seg;
                bool residue_in_segment = false;

                if (not residue.segment.empty())
                {
                    // this residue belongs in a segment
                    auto segname = QString::fromStdString(residue.segment);

                    if (not segments.contains(segname))
                    {
                        // create a new segment
                        auto seg = cg.molecule().add(SireMol::SegName(segname));
                        segments.insert(segname, seg);
                    }

                    seg = segments[segname];
                    residue_in_segment = true;
                }

                QString is_hetatm("False");

                if (residue.het_flag == 'H')
                    is_hetatm = "True";

                for (const auto &atom : residue.atoms)
                {
                    auto atm = cg.add(SireMol::AtomNum(atom.serial));
                    atm.reparent(res.index());
                    atm.rename(SireMol::AtomName(QString::fromStdString(atom.name)));
                    populate_atom(atm, atom, is_hetatm, map);

                    serial_to_molnum.insert(atom.serial, molnum);

                    if (residue_in_segment)
                        atm.reparent(seg.index());
                }
            }
        }

        mols.add(mol.commit());
    }

    void parse_molecules(const gemmi::Entity &entity, SireMol::MoleculeGroup &mols,
                         const QHash<QString, int> &subchains,
                         const gemmi::Structure &structure,
                         QHash<int, SireMol::MolNum> &serial_to_molnum,
                         const SireBase::PropertyMap &map)
    {
        // parse each chain into a separate molecule

        for (const auto &subchain : entity.subchains)
        {
            auto model_id = subchains.value(QString::fromStdString(subchain), 0);

            auto mol = SireMol::MolStructureEditor(SireMol::Molecule().edit());
            mol.renumber();

            auto cg0 = mol.add(SireMol::CGName("0"));

            const auto molnum = mol.number();

            QHash<QString, SireMol::SegStructureEditor> segments;

            int cg_num = 0;

            for (const auto &residue : structure.models[model_id].get_subchain(subchain))
            {
                auto cg = cg0;

                if (cg_num == 0)
                {
                    // first residue - use the residue name as the atom name
                    mol.rename(SireMol::MolName(QString::fromStdString(residue.name)));
                }
                else
                {
                    cg = mol.add(SireMol::CGName(QString::number(cg_num)));
                }

                cg_num += 1;

                int resnum = cg_num;

                if (residue.seqid.num.has_value())
                    resnum = residue.seqid.num.value;

                auto res = mol.add(SireMol::ResNum(resnum));
                res.rename(SireMol::ResName(QString::fromStdString(residue.name)));

                if (residue.seqid.has_icode())
                    set_prop(res, "insert_code", QString(residue.seqid.icode), map);

                SireMol::SegStructureEditor seg;
                bool residue_in_segment = false;

                if (not residue.segment.empty())
                {
                    // this residue belongs in a segment
                    auto segname = QString::fromStdString(residue.segment);

                    if (not segments.contains(segname))
                    {
                        // create a new segment
                        auto seg = mol.add(SireMol::SegName(segname));
                        segments.insert(segname, seg);
                    }

                    seg = segments[segname];
                    residue_in_segment = true;
                }

                QString is_hetatm("False");

                if (residue.het_flag == 'H')
                    is_hetatm = "True";

                for (const auto &atom : residue.atoms)
                {
                    auto atm = cg.add(SireMol::AtomNum(atom.serial));
                    atm.reparent(res.index());
                    atm.rename(SireMol::AtomName(QString::fromStdString(atom.name)));
                    populate_atom(atm, atom, is_hetatm, map);

                    serial_to_molnum.insert(atom.serial, molnum);

                    if (residue_in_segment)
                        atm.reparent(seg.index());
                }
            }

            auto m = mol.commit();
            mols.add(m);
        }
    }

    void parse_waters(const gemmi::Entity &entity, SireMol::MoleculeGroup &mols,
                      const QHash<QString, int> &subchains,
                      const gemmi::Structure &structure,
                      QHash<int, SireMol::MolNum> &serial_to_molnum,
                      const SireBase::PropertyMap &map)
    {
        // parse each residue in each chain into a separate water molecule

        for (const auto &subchain : entity.subchains)
        {
            auto model_id = subchains.value(QString::fromStdString(subchain), 0);

            int num_waters = 0;

            for (const auto &residue : structure.models[model_id].get_subchain(subchain))
            {
                num_waters += 1;

                auto mol = SireMol::MolStructureEditor(SireMol::Molecule().edit());
                mol.renumber();
                mol.rename(SireMol::MolName("WAT"));

                auto cg = mol.add(SireMol::CGName("0"));

                const auto molnum = cg.molecule().number();

                int resnum = num_waters;

                if (residue.seqid.num.has_value())
                    resnum = residue.seqid.num.value;

                auto res = mol.add(SireMol::ResNum(resnum));
                res.rename(SireMol::ResName(QString::fromStdString(residue.name)));

                if (residue.seqid.has_icode())
                    set_prop(res, "insert_code", QString(residue.seqid.icode), map);

                QString is_hetatm("False");

                if (residue.het_flag == 'H')
                    is_hetatm = "True";

                for (const auto &atom : residue.atoms)
                {
                    auto atm = cg.add(SireMol::AtomNum(atom.serial));
                    atm.reparent(res.index());
                    atm.rename(SireMol::AtomName(QString::fromStdString(atom.name)));
                    populate_atom(atm, atom, is_hetatm, map);

                    serial_to_molnum.insert(atom.serial, molnum);
                }

                auto m = mol.commit().edit();
                m.setProperty(map["is_water"], SireBase::BooleanProperty(true));

                mols.add(m.commit());
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

        QHash<int, SireMol::MolNum> serial_to_molnum;

        // this function doesn't appear to be in the gemmi API?
        // structure.standardize_crystal_frame();

        auto structure_name = QString::fromStdString(structure.name);

        QStringList polynames;

        for (const auto &entity : structure.entities)
        {
            if (entity.entity_type == gemmi::EntityType::Polymer)
            {
                polynames.append(structure_name);
            }
        }

        if (polynames.size() > 1)
        {
            int npolys = 0;

            for (const auto &entity : structure.entities)
            {
                if (entity.entity_type == gemmi::EntityType::Polymer)
                {
                    polynames[npolys] = QString("[%1]-%2")
                                            .arg(QString::fromStdString(entity.name))
                                            .arg(polynames[npolys]);
                    npolys += 1;
                }
            }
        }

        int npolys = 0;

        for (const auto &entity : structure.entities)
        {
            switch (entity.entity_type)
            {
            case gemmi::EntityType::Polymer:
                parse_polymer(entity, mols, subchains, structure,
                              serial_to_molnum, polynames[npolys], map);
                npolys += 1;
                break;
            case gemmi::EntityType::NonPolymer:
                parse_molecules(entity, mols, subchains, structure,
                                serial_to_molnum, map);
                break;
            case gemmi::EntityType::Branched:
                parse_molecules(entity, mols, subchains, structure,
                                serial_to_molnum, map);
                break;
            case gemmi::EntityType::Water:
                parse_waters(entity, mols, subchains, structure,
                             serial_to_molnum, map);
                break;
            default:
                parse_molecules(entity, mols, subchains, structure,
                                serial_to_molnum, map);
            }
        }

        // now need to handle any connections
        const auto connectivity_property = map["connectivity"];

        if (connectivity_property.hasSource())
        {
            bool auto_connect = true;

            if (map.specified("auto_connect"))
            {
                auto_connect = map["auto_connect"].value().asABoolean();
            }

            QHash<SireMol::MolNum, SireMol::ConnectivityEditor> connectivities;

            for (const auto &connection : structure.connections)
            {
                const auto model_num = subchains.value(QString::fromStdString(connection.partner1.chain_name), 0);

                // the connection must be between two atoms in the same model
                const auto &atom1 = structure.models[model_num].find_cra(connection.partner1);
                const auto &atom2 = structure.models[model_num].find_cra(connection.partner2);

                if (atom1.atom == 0 or atom2.atom == 0)
                {
                    // one of the atoms is missing
                    continue;
                }

                // get the atom numbers of both atoms
                const auto serial1 = atom1.atom->serial;
                const auto serial2 = atom2.atom->serial;

                // find the atoms - the serial is the AtomNum and should be unique
                const auto molnum1 = serial_to_molnum.value(serial1, SireMol::MolNum(0));
                const auto molnum2 = serial_to_molnum.value(serial2, SireMol::MolNum(0));

                if (molnum1 == SireMol::MolNum(0) or molnum1 != molnum2)
                {
                    // bonds between molecules? Or something else wrong?
                    continue;
                }

                // get the connectivity for this molecule
                if (not connectivities.contains(molnum1))
                {
                    const auto mol = mols.molecule(molnum1).molecule();

                    // create a new connectivity
                    if (auto_connect)
                    {
                        try
                        {
                            // if we are auto-connecting, then we need to auto-generate
                            // the connectivity now, so that the CONECT records are added
                            // to this. If we don't, then the later stage of code won't
                            // want to overwrite the connectivity
                            auto hunter = SireMol::CovalentBondHunter();
                            connectivities.insert(molnum1, hunter(mol).edit());
                        }
                        catch (...)
                        {
                            connectivities.insert(molnum1, SireMol::ConnectivityEditor(mol));
                        }
                    }
                    else
                    {
                        connectivities.insert(molnum1, SireMol::ConnectivityEditor(mol));
                    }
                }

                auto &connectivity = connectivities[molnum1];
                connectivity.connect(SireMol::AtomNum(serial1), SireMol::AtomNum(serial2));
            }

            for (auto it = connectivities.constBegin(); it != connectivities.constEnd(); ++it)
            {
                auto connectivity = it.value().commit();

                auto molnum = it.key();

                auto mol = mols.molecule(molnum).molecule().edit();

                mol.setProperty(connectivity_property, connectivity);

                mols.update(mol.commit());
            }
        }

        SireSystem::System system;
        system.setName(structure_name);

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
        auto structure = sire_to_gemmi(system, map);

        auto doc = gemmi::make_mmcif_document(structure);

        std::string s;
        std::stringstream stream(s);

        gemmi::cif::write_cif_to_stream(stream, doc);

        return QString::fromStdString(s).split("\n");
    }

    void register_pdbx_loader()
    {
        SireIO::PDBxReaderFunction reader_function(&pdbx_reader_function);
        SireIO::PDBxWriterFunction writer_function(&pdbx_writer_function);

        SireIO::register_pdbx_loader_functions(writer_function, reader_function);
    }
}
