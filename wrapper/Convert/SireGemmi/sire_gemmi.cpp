
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

#include "SireMol/iswater.h"

#include "SireBase/propertylist.h"
#include "SireBase/stringproperty.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include <string>
#include <sstream>
#include <iostream>

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

            auto subchain_residues = structure.models[model_id].get_subchain(subchain);

            if (subchain_residues.empty())
            {
                continue;
            }

            auto chain = mol.add(SireMol::ChainName(QString::fromStdString(subchain.c_str())));

            for (const auto &residue : subchain_residues)
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

                auto res = chain.add(SireMol::ResName(QString::fromStdString(residue.name)));
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

            auto chain = mol.add(SireMol::ChainName(QString::fromStdString(subchain.c_str())));

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
                res.reparent(chain.name());

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
                auto chain = mol.add(SireMol::ChainName(QString::fromStdString(subchain.c_str())));

                auto cg = mol.add(SireMol::CGName("0"));

                const auto molnum = cg.molecule().number();

                int resnum = num_waters;

                if (residue.seqid.num.has_value())
                    resnum = residue.seqid.num.value;

                auto res = mol.add(SireMol::ResNum(resnum));
                res.rename(SireMol::ResName(QString::fromStdString(residue.name)));
                res.reparent(chain.name());

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

        gemmi::ensure_entities(structure);

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

        if (map.specified("metadata"))
        {
            auto metadata = map["metadata"].value().asA<SireBase::Properties>();
            system.setProperty("metadata", metadata);
        }

        return system;
    }

    bool populate_atom(gemmi::Atom &gemmi_atom,
                       const SireMol::Atom &atom,
                       const SireBase::PropertyMap &map)
    {
        gemmi_atom.name = atom.name().value().toStdString();
        gemmi_atom.serial = atom.number().value();

        auto coords = atom.property<SireMaths::Vector>(map["coordinates"]);
        gemmi_atom.pos = gemmi::Position(coords.x(), coords.y(), coords.z());

        try
        {
            auto element = atom.property<SireMol::Element>(map["element"]);
            gemmi_atom.element = gemmi::Element(element.nProtons());
        }
        catch (...)
        {
        }

        try
        {
            auto chg = atom.property<SireUnits::Dimension::Charge>(map["formal_charge"]).to(SireUnits::mod_electron);
            gemmi_atom.charge = int(chg);
        }
        catch (...)
        {
        }

        try
        {
            auto occ = atom.property<double>(map["occupancy"]);
            gemmi_atom.occ = occ;
        }
        catch (...)
        {
        }

        try
        {
            auto b_iso = atom.property<double>(map["beta_factor"]);
            gemmi_atom.b_iso = b_iso;
        }
        catch (...)
        {
        }

        bool is_hetatm = false;

        try
        {
            auto is_hetatm = atom.property<QString>(map["is_het"]);

            if (is_hetatm == "True")
                is_hetatm = true;
        }
        catch (...)
        {
        }

        try
        {
            auto alt_loc = atom.property<QString>(map["alt_loc"]);

            if (not alt_loc.isEmpty())
                gemmi_atom.altloc = alt_loc.toStdString()[0];
        }
        catch (...)
        {
        }

        return is_hetatm;
    }

    void convert_polymer(int molid,
                         const SireMol::Molecule &mol, gemmi::Chain &chain,
                         QHash<QString, gemmi::Chain> &chains,
                         const SireBase::PropertyMap &map)
    {
        const auto residues = mol.residues();

        const auto entity_id = QString::number(molid).toStdString();

        for (int i = 0; i < residues.count(); ++i)
        {
            const auto residue = residues(i);

            gemmi::Residue gemmi_residue;
            gemmi_residue.entity_type = gemmi::EntityType::Polymer;
            gemmi_residue.entity_id = entity_id;

            gemmi_residue.name = residue.name().value().toStdString();
            gemmi_residue.seqid.num = residue.number().value();

            const auto atoms = residue.atoms();

            bool is_hetatm_residue = false;

            QString seg;

            for (int j = 0; j < atoms.count(); ++j)
            {
                const auto atom = atoms(j);

                gemmi::Atom gemmi_atom;
                auto is_hetatm = populate_atom(gemmi_atom, atom, map);

                is_hetatm_residue = is_hetatm_residue or is_hetatm;

                gemmi_residue.atoms.push_back(gemmi_atom);

                if (atom.isWithinSegment())
                    seg = atom.segment().name().value();
            }

            if (not seg.isEmpty())
                gemmi_residue.segment = seg.toStdString();

            if (is_hetatm_residue)
                gemmi_residue.het_flag = 'H';
            else
                gemmi_residue.het_flag = 'A';

            if (residue.isWithinChain())
            {
                gemmi_residue.subchain = residue.chain().name().value().toStdString();
                chains.find(residue.chain().name().value()).value().residues.push_back(gemmi_residue);
            }
            else
            {
                chain.residues.push_back(gemmi_residue);
            }
        }
    }

    void convert_molecule(int molid,
                          const SireMol::Molecule &mol, gemmi::Chain &chain,
                          QHash<QString, gemmi::Chain> &chains,
                          const SireBase::PropertyMap &map)
    {
        const auto residues = mol.residues();

        const auto entity_id = QString::number(molid).toStdString();

        for (int i = 0; i < residues.count(); ++i)
        {
            const auto residue = residues(i);

            gemmi::Residue gemmi_residue;
            gemmi_residue.entity_type = gemmi::EntityType::NonPolymer;
            gemmi_residue.entity_id = entity_id;

            gemmi_residue.name = residue.name().value().toStdString();
            gemmi_residue.seqid.num = residue.number().value();

            const auto atoms = residue.atoms();

            bool is_hetatm_residue = false;

            QString seg;

            for (int j = 0; j < atoms.count(); ++j)
            {
                const auto atom = atoms(j);

                gemmi::Atom gemmi_atom;
                auto is_hetatm = populate_atom(gemmi_atom, atom, map);

                is_hetatm_residue = is_hetatm_residue or is_hetatm;

                if (atom.isWithinSegment())
                    seg = atom.segment().name().value();

                gemmi_residue.atoms.push_back(gemmi_atom);
            }

            if (not seg.isEmpty())
                gemmi_residue.segment = seg.toStdString();

            if (is_hetatm_residue)
                gemmi_residue.het_flag = 'H';
            else
                gemmi_residue.het_flag = 'A';

            if (residue.isWithinChain())
            {
                gemmi_residue.subchain = residue.chain().name().value().toStdString();
                chains.find(residue.chain().name().value()).value().residues.push_back(gemmi_residue);
            }
            else
            {
                chain.residues.push_back(gemmi_residue);
            }
        }
    }

    void convert_water(int molid,
                       const SireMol::Molecule &mol, gemmi::Chain &chain,
                       QHash<QString, gemmi::Chain> &chains,
                       const SireBase::PropertyMap &map)
    {
        gemmi::Residue residue;
        residue.entity_type = gemmi::EntityType::Water;
        residue.entity_id = QString::number(molid).toStdString();

        const auto atoms = mol.atoms();

        if (atoms.isEmpty())
            return;

        auto first_res = atoms(0).residue();
        residue.name = first_res.name().value().toStdString();
        residue.seqid.num = first_res.number().value();

        bool is_hetatm_residue = false;

        for (int i = 0; i < atoms.count(); ++i)
        {
            const auto atom = atoms(i);

            gemmi::Atom gemmi_atom;
            auto is_hetatm = populate_atom(gemmi_atom, atom, map);

            is_hetatm_residue = is_hetatm_residue or is_hetatm;

            residue.atoms.push_back(gemmi_atom);
        }

        if (is_hetatm_residue)
            residue.het_flag = 'H';
        else
            residue.het_flag = 'A';

        if (first_res.isWithinChain())
        {
            residue.subchain = first_res.chain().name().value().toStdString();
            chains.find(first_res.chain().name().value()).value().residues.push_back(residue);
        }
        else
        {
            chain.residues.push_back(residue);
        }
    }

    gemmi::Structure sire_to_gemmi(const SireSystem::System &system,
                                   const SireBase::PropertyMap &map)
    {
        if (system.nAtoms() == 0)
            return gemmi::Structure();

        const auto mols = SireMol::SelectorMol(system);

        gemmi::Structure structure;
        gemmi::Model model(system.name().value().toStdString());

        auto name = system.name().value().toStdString();

        if (name.empty())
        {
            // try to find a name from the molecules
            for (const auto &mol : mols)
            {
                name = mol.name().value().toStdString();

                if (not name.empty())
                    break;
            }

            if (name.empty())
                name = "sire_system";
        }

        structure.name = name;

        // create a real chain for each named chain in the molecules
        QHash<QString, gemmi::Chain> chains;

        for (const auto &mol : mols)
        {
            if (mol.nChains() > 0)
            {
                const auto c = mol.chains();

                for (int i = 0; i < c.count(); ++i)
                {
                    auto name = c(i).name().value();

                    if (not chains.contains(name))
                        chains.insert(name, gemmi::Chain(name.toStdString()));
                }
            }
        }

        // also create a single chain for all other molecules
        int chain_num = 1;
        auto chain_name = QString::number(chain_num).toStdString();

        while (true)
        {
            if (not chains.contains(QString::number(chain_num)))
            {
                break;
            }

            chain_num += 1;
            chain_name = QString::number(chain_num).toStdString();
        }

        auto chain = gemmi::Chain(chain_name);

        QHash<SireMol::MolNum, int> molnum_to_chain;

        int molid = 0;

        for (const auto &mol : mols)
        {
            molid += 1;

            if (SireMol::is_water(mol, map))
            {
                convert_water(molid, mol, chain, chains, map);
            }
            else if (mol.nAtoms() == 1)
            {
                if (mol.atoms()(0).property<SireMol::Element>(map["element"]) == SireMol::Element(8))
                {
                    // single oxygen is a water without hydrogens
                    convert_water(molid, mol, chain, chains, map);
                }
                else
                {
                    // convert as a normal molecule
                    convert_molecule(molid, mol, chain, chains, map);
                }
            }
            else if (mol.nResidues() > 1)
            {
                convert_polymer(molid, mol, chain, chains, map);
            }
            else
            {
                // convert as a normal molecule
                convert_molecule(molid, mol, chain, chains, map);
            }
        }

        auto chain_names = chains.keys();
        chain_names.sort();

        for (const auto &name : chain_names)
        {
            auto &chain = chains.find(name).value();

            if (not chain.residues.empty())
                model.chains.push_back(chain);
        }

        if (not chain.residues.empty())
            model.chains.push_back(chain);

        structure.models.push_back(model);

        structure.renumber_models();
        gemmi::setup_entities(structure);
        gemmi::assign_serial_numbers(structure, true);

        // now we have done this, we need to add in all of the bonds
        const auto connectivity_property = map["connectivity"];

        for (const auto &mol : mols)
        {
            try
            {
                const auto atoms = mol.atoms();
                const int nats = atoms.count();

                std::vector<std::string> chain_names;

                const auto residues = mol.residues();

                for (int i = 0; i < residues.count(); ++i)
                {
                    const auto res = residues(i);

                    if (res.isWithinChain())
                    {
                        auto name = res.chain().name().value().toStdString();
                        chain_names.push_back(name);
                    }
                    else
                    {
                        chain_names.push_back(chain_name);
                    }
                }

                const auto &connectivity = mol.property(connectivity_property).asA<SireMol::Connectivity>();

                for (int i = 0; i < nats - 1; ++i)
                {
                    const auto idx0 = SireMol::AtomIdx(i);
                    const auto atom0 = atoms(i);

                    auto connections = connectivity.connectionsTo(idx0);

                    if (connections.isEmpty())
                        continue;

                    gemmi::AtomAddress addr0;

                    auto res0 = atom0.residue();

                    addr0.chain_name = chain_names[res0.index().value()];
                    addr0.res_id.name = res0.name().value().toStdString();
                    addr0.res_id.seqid.num = res0.number().value();
                    addr0.atom_name = atom0.name().value().toStdString();

                    auto cra0 = model.find_cra(addr0, true);

                    for (const auto &idx1 : connections)
                    {
                        if (idx1.value() <= i)
                            continue;

                        const auto atom1 = atoms(idx1.value());
                        const auto res1 = atom1.residue();

                        gemmi::AtomAddress addr1;

                        addr1.chain_name = chain_names[res1.index().value()];
                        addr1.res_id.name = res1.name().value().toStdString();
                        addr1.res_id.seqid.num = res1.number().value();
                        addr1.atom_name = atom1.name().value().toStdString();

                        auto cra1 = model.find_cra(addr1, true);

                        gemmi::Connection connection;
                        connection.partner1 = addr0;
                        connection.partner2 = addr1;

                        structure.connections.push_back(connection);
                    }
                }
            }
            catch (...)
            {
                // no connectivity for this molecule
                continue;
            }
        }

        return structure;
    }

    QString _string_to_property(const std::string &s)
    {
        QString q = QString::fromStdString(s);

        if (q.startsWith("\"") and q.endsWith("\""))
            q = q.mid(1, q.size() - 2);
        else if (q.startsWith("'") and q.endsWith("'"))
            q = q.mid(1, q.size() - 2);

        return q;
    }

    SireSystem::System pdbx_reader_function(const QStringList &lines,
                                            const SireBase::PropertyMap &map)
    {
        // assemble all of the line into a single string
        auto input_string = lines.join("\n").toStdString();

        auto doc = cif::Document(cif::read_string(input_string));

        input_string.clear();

        int structure_block = -1;

        // mmCIF files for deposition may have more than one block:
        // coordinates in the first block and restraints in the others.
        for (size_t i = 0; i < doc.blocks.size(); ++i)
        {
            if (doc.blocks[i].has_tag("_atom_site.id"))
            {
                if (structure_block != -1)
                    throw SireError::unsupported(
                        QObject::tr("2+ blocks are ok if only the first one has coordinates. "
                                    "_atom_site in block #%1 : %2")
                            .arg(i + 1)
                            .arg(QString::fromStdString(doc.source)),
                        CODELOC);

                structure_block = i;
            }
        }

        if (structure_block == -1)
            // just use the first block
            structure_block = 0;

        auto structure = gemmi::make_structure_from_block(doc.blocks.at(structure_block));

        // check for any metadata - if there is, then add it to the map
        auto *block = doc.find_block("sire");
        SireBase::Properties metadata;

        if (block != nullptr)
        {
            for (const auto &item : block->items)
            {
                switch (item.type)
                {
                case cif::ItemType::Pair:
                {
                    const auto &p = item.pair;

                    QString tag = QString::fromStdString(p[0]);
                    QString value = _string_to_property(p[1]);

                    if (tag.startsWith("_"))
                        tag = tag.mid(1);

                    metadata.setProperty(tag, SireBase::StringProperty(value));
                    break;
                }
                case cif::ItemType::Loop:
                {
                    const auto &l = item.loop;

                    QStringList tags;

                    for (const auto &tag : l.tags)
                    {
                        tags.append(QString::fromStdString(tag));
                    }

                    QStringList values;

                    for (const auto &value : l.values)
                    {
                        values.append(_string_to_property(value));
                    }

                    if (tags.size() == 1 and tags[0].endsWith(".value"))
                    {
                        // this is an array of values
                        auto array = SireBase::PropertyList();

                        for (const auto &value : values)
                        {
                            array.append(SireBase::StringProperty(value));
                        }

                        auto tag = tags[0];

                        if (tag.startsWith("_"))
                            tag = tag.mid(1);

                        // remove the .value on the end
                        tag = tag.left(tag.size() - 6);

                        metadata.setProperty(tag, array);
                    }
                    else
                    {
                        // this is a set of property values
                        auto props = SireBase::Properties();

                        if (tags.size() == values.size())
                        {
                            // one value per key
                            for (int i = 0; i < tags.size(); ++i)
                            {
                                auto subtag = tags[i].split(".").mid(1).join(".");
                                props.setProperty(subtag, SireBase::StringProperty(values[i]));
                            }
                        }
                        else
                        {
                            // multiple values per key
                            for (int i = 0; i < tags.size(); ++i)
                            {
                                auto subtag = tags[i].split(".").mid(1).join(".");

                                auto subvals = SireBase::StringArrayProperty();

                                for (int j = i; j < values.size(); j += tags.size())
                                {
                                    subvals.append(values[j]);
                                }

                                props.setProperty(subtag, subvals);
                            }
                        }

                        auto tag = tags[0].split(".").at(0);

                        if (tag.startsWith("_"))
                            tag = tag.mid(1);

                        metadata.setProperty(tag, props);
                    }

                    break;
                }
                }
            }
        }

        auto m = map;

        if (not metadata.isEmpty())
        {
            m.set("metadata", metadata);
        }

        return gemmi_to_sire(structure, m);
    }

    std::string _property_to_string(const SireBase::Property &p)
    {
        QString s;

        try
        {
            s = p.asAString();
        }
        catch (...)
        {
            s = p.toString();
        }

        // make sure we put any strings that contain spaces into quotes
        if (s.contains(" "))
            s = QString("\"%1\"").arg(s);

        return s.toStdString();
    }

    QStringList pdbx_writer_function(const SireSystem::System &system,
                                     const SireBase::PropertyMap &map)
    {
        auto structure = sire_to_gemmi(system, map);

        auto doc = gemmi::make_mmcif_document(structure);

        if (system.containsProperty(map["metadata"]))
        {
            auto &block = doc.add_new_block("sire");

            auto metadata = system.property(map["metadata"]).asA<SireBase::Properties>();

            auto keys = metadata.propertyKeys();
            keys.sort();

            for (const auto &key : keys)
            {
                const auto &value = metadata.property(key);

                if (value.isA<SireBase::Properties>())
                {
                    const auto &props2 = value.asA<SireBase::Properties>();
                    auto keys2 = props2.propertyKeys();
                    keys2.sort();

                    std::vector<std::string> tags;

                    int nrows = 0;

                    for (const auto &key2 : keys2)
                    {
                        tags.push_back(QString(".%1").arg(key2.simplified().replace(" ", "_")).toStdString());

                        const auto &value2 = props2.property(key2);

                        if (value2.isAnArray())
                        {
                            nrows = std::max(nrows, value2.asAnArray().count());
                        }
                    }

                    auto &loop = block.init_loop(QString("_%1").arg(key.simplified().replace(" ", "_")).toStdString(), tags);

                    for (int i = 0; i < nrows; ++i)
                    {
                        std::vector<std::string> values;

                        for (const auto &key2 : keys2)
                        {
                            const auto &value2 = props2.property(key2);

                            if (value2.isAnArray())
                            {
                                const auto &array2 = value2.asAnArray();

                                if (i < array2.count())
                                {
                                    values.push_back(_property_to_string(array2[i]));
                                }
                                else
                                {
                                    values.push_back("\"\"");
                                }
                            }
                            else
                            {
                                values.push_back(_property_to_string(value2));
                            }
                        }

                        loop.add_row(values);
                    }
                }
                else if (value.isAnArray())
                {
                    auto tag = QString("_%1").arg(key.simplified().replace(" ", "_"));

                    auto array = value.asAnArray();

                    if (array.count() == 1)
                    {
                        block.set_pair(tag.toStdString(), _property_to_string(array[0]));
                    }
                    else
                    {
                        auto &loop = block.init_loop(tag.toStdString(), {".value"});

                        for (int i = 0; i < array.size(); ++i)
                        {
                            loop.add_row({_property_to_string(array[i])});
                        }
                    }
                }
                else
                {
                    auto tag = QString("_%1").arg(key.simplified().replace(" ", "_"));

                    doc.blocks[0].set_pair(tag.toStdString(), _property_to_string(value));
                }
            }
        }

        std::stringstream stream;

        gemmi::cif::write_cif_to_stream(stream, doc);

        stream.flush();

        auto lines = QString::fromStdString(stream.str()).split("\n");

        return lines;
    }

    void register_pdbx_loader()
    {
        SireIO::PDBxReaderFunction reader_function(&pdbx_reader_function);
        SireIO::PDBxWriterFunction writer_function(&pdbx_writer_function);

        SireIO::register_pdbx_loader_functions(writer_function, reader_function);
    }
}
