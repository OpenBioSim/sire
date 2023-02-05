
#include "sire_rdkit.h"

#include "GraphMol/molops.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomproperty.hpp"
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"

#include "SireMM/selectorbond.h"

#include "SireMaths/vector.h"

#include "SireBase/propertylist.h"
#include "SireBase/stringproperty.h"

#include "SireUnits/units.h"

#include <map>

#include <QDebug>

using RDKit::ROMOL_SPTR;
using SireBase::PropertyMap;
using SireMol::Molecule;
using SireMol::SelectorMol;

RDKit::Atom::ChiralType string_to_chiral(const QString &typ)
{
    const static std::map<QString, RDKit::Atom::ChiralType> types = {
        {"CHI_UNSPECIFIED", RDKit::Atom::CHI_UNSPECIFIED},
        {"CHI_TETRAHEDRAL_CW", RDKit::Atom::CHI_TETRAHEDRAL_CW},
        {"CHI_TETRAHEDRAL_CCW", RDKit::Atom::CHI_TETRAHEDRAL_CCW},
        {"CHI_OTHER", RDKit::Atom::CHI_OTHER},
        {"CHI_TETRAHEDRAL", RDKit::Atom::CHI_OTHER},
        {"CHI_ALLENE", RDKit::Atom::CHI_ALLENE},
        {"CHI_SQUAREPLANAR", RDKit::Atom::CHI_SQUAREPLANAR},
        {"CHI_TRIGONALBIPYRAMIDAL", RDKit::Atom::CHI_TRIGONALBIPYRAMIDAL},
        {"CHI_OCTAHEDRAL", RDKit::Atom::CHI_OCTAHEDRAL},
        {"CHI_UNKNOWN", RDKit::Atom::CHI_UNSPECIFIED}};

    auto it = types.find(typ);

    if (it == types.end())
        return RDKit::Atom::CHI_UNSPECIFIED;
    else
        return it->second;
}

QString chiral_to_string(RDKit::Atom::ChiralType typ)
{
    switch (typ)
    {
    case RDKit::Atom::CHI_UNSPECIFIED:
        return "CHI_UNSPECIFIED";
    case RDKit::Atom::CHI_TETRAHEDRAL_CW:
        return "CHI_TETRAHEDRAL_CW";
    case RDKit::Atom::CHI_TETRAHEDRAL_CCW:
        return "CHI_TETRAHEDRAL_CCW";
    case RDKit::Atom::CHI_OTHER:
        return "CHI_OTHER";
    case RDKit::Atom::CHI_TETRAHEDRAL:
        return "CHI_TETRAHEDRAL";
    case RDKit::Atom::CHI_ALLENE:
        return "CHI_ALLENE";
    case RDKit::Atom::CHI_SQUAREPLANAR:
        return "CHI_SQUAREPLANAR";
    case RDKit::Atom::CHI_TRIGONALBIPYRAMIDAL:
        return "CHI_TRIGONALBIPYRAMIDAL";
    case RDKit::Atom::CHI_OCTAHEDRAL:
        return "CHI_OCTAHEDRAL";
    default:
        return "CHI_UNKNOWN";
    }
}

RDKit::Atom::HybridizationType string_to_hybridization(const QString &typ)
{
    const static std::map<QString, RDKit::Atom::HybridizationType> types = {
        {"UNSPECIFIED", RDKit::Atom::UNSPECIFIED},
        {"S", RDKit::Atom::S},
        {"SP", RDKit::Atom::SP},
        {"SP2", RDKit::Atom::SP2},
        {"SP2D", RDKit::Atom::SP2D},
        {"SP3D", RDKit::Atom::SP3D},
        {"SP3D2", RDKit::Atom::SP3D2},
        {"OTHER", RDKit::Atom::OTHER},
        {"UNKNOWN", RDKit::Atom::UNSPECIFIED}};

    auto it = types.find(typ);

    if (it == types.end())
        return RDKit::Atom::UNSPECIFIED;
    else
        return it->second;
}

QString hybridization_to_string(RDKit::Atom::HybridizationType typ)
{
    switch (typ)
    {
    case RDKit::Atom::UNSPECIFIED:
        return "UNSPECIFIED";
    case RDKit::Atom::S:
        return "S";
    case RDKit::Atom::SP:
        return "SP";
    case RDKit::Atom::SP2:
        return "SP2";
    case RDKit::Atom::SP2D:
        return "SP2D";
    case RDKit::Atom::SP3D:
        return "SP3D";
    case RDKit::Atom::SP3D2:
        return "SP3D2";
    case RDKit::Atom::OTHER:
        return "OTHER";
    default:
        return "UNKNOWN";
    }
}

QString bondtype_to_string(RDKit::Bond::BondType typ)
{
    switch (typ)
    {
    case RDKit::Bond::UNSPECIFIED:
        return "UNSPECIFIED";
    case RDKit::Bond::SINGLE:
        return "SINGLE";
    case RDKit::Bond::DOUBLE:
        return "DOUBLE";
    case RDKit::Bond::TRIPLE:
        return "TRIPLE";
    case RDKit::Bond::QUADRUPLE:
        return "QUADRUPLE";
    case RDKit::Bond::QUINTUPLE:
        return "QUINTUPLE";
    case RDKit::Bond::HEXTUPLE:
        return "HEXTUPLE";
    case RDKit::Bond::ONEANDAHALF:
        return "ONEANDAHALF";
    case RDKit::Bond::TWOANDAHALF:
        return "TWOANDAHALF";
    case RDKit::Bond::THREEANDAHALF:
        return "THREEANDAHALF";
    case RDKit::Bond::FOURANDAHALF:
        return "FOURANDAHALF";
    case RDKit::Bond::FIVEANDAHALF:
        return "FIVEANDAHALF";
    case RDKit::Bond::AROMATIC:
        return "AROMATIC";
    case RDKit::Bond::IONIC:
        return "IONIC";
    case RDKit::Bond::HYDROGEN:
        return "HYDROGEN";
    case RDKit::Bond::THREECENTER:
        return "THREECENTER";
    case RDKit::Bond::DATIVEONE:
        return "DATIVEONE";
    case RDKit::Bond::DATIVE:
        return "DATIVEONE";
    case RDKit::Bond::DATIVEL:
        return "DATIVEL";
    case RDKit::Bond::DATIVER:
        return "DATIVER";
    case RDKit::Bond::OTHER:
        return "OTHER";
    case RDKit::Bond::ZERO:
        return "ZERO";
    default:
        return "UNSPECIFIED";
    }
}

RDKit::Bond::BondType string_to_bondtype(const QString &typ)
{
    const static std::map<QString, RDKit::Bond::BondType> types = {
        {"UNSPECIFIED", RDKit::Bond::UNSPECIFIED},
        {"SINGLE", RDKit::Bond::SINGLE},
        {"DOUBLE", RDKit::Bond::DOUBLE},
        {"TRIPLE", RDKit::Bond::TRIPLE},
        {"QUADRUPLE", RDKit::Bond::QUADRUPLE},
        {"QUINTUPLE", RDKit::Bond::QUINTUPLE},
        {"HEXTUPLE", RDKit::Bond::HEXTUPLE},
        {"ONEANDAHALF", RDKit::Bond::ONEANDAHALF},
        {"TWOANDAHALF", RDKit::Bond::TWOANDAHALF},
        {"THREEANDAHALF", RDKit::Bond::THREEANDAHALF},
        {"FOURANDAHALF", RDKit::Bond::FOURANDAHALF},
        {"FIVEANDAHALF", RDKit::Bond::FIVEANDAHALF},
        {"AROMATIC", RDKit::Bond::AROMATIC},
        {"IONIC", RDKit::Bond::IONIC},
        {"HYDROGEN", RDKit::Bond::HYDROGEN},
        {"THREECENTER", RDKit::Bond::THREECENTER},
        {"DATIVEONE", RDKit::Bond::DATIVEONE},
        {"DATIVE", RDKit::Bond::DATIVE},
        {"DATIVEL", RDKit::Bond::DATIVEL},
        {"DATIVER", RDKit::Bond::DATIVER},
        {"OTHER", RDKit::Bond::OTHER},
        {"ZERO", RDKit::Bond::ZERO}};

    auto it = types.find(typ);

    if (it == types.end())
        return RDKit::Bond::UNSPECIFIED;
    else
        return it->second;
}

template <class T, class V>
void set_prop(T &atom, const QString &key, const V &value,
              const PropertyMap &map)
{
    const auto name = map[key];

    if (name.hasSource())
        atom.setProperty(name.source(), value);
}

ROMOL_SPTR sire_to_rdkit(const Molecule &mol, const PropertyMap &map)
{
    RDKit::RWMol molecule;
    molecule.beginBatchEdit();

    molecule.setProp<std::string>("_Name", mol.name().value().toStdString());

    const auto atoms = mol.atoms();

    for (int i = 0; i < atoms.count(); ++i)
    {
        const auto atom = atoms(i);

        molecule.addAtom(true);
        auto a = molecule.getActiveAtom();

        a->setProp<std::string>("_Name", atom.name().value().toStdString());
        a->setAtomicNum(atom.property<SireMol::Element>(map["element"]).nProtons());

        try
        {
            a->setFormalCharge(atom.property<SireUnits::Dimension::Charge>(map["formal_charge"]).to(SireUnits::mod_electron));
        }
        catch (...)
        {
        }

        try
        {
            a->setIsAromatic(atom.property<qint64>(map["is_aromatic"]) != 0);
        }
        catch (...)
        {
        }

        try
        {
            a->setIsotope(atom.property<qint64>(map["isotope"]));
        }
        catch (...)
        {
        }

        try
        {
            a->setChiralTag(string_to_chiral(atom.property<QString>(map["chiral_tag"])));
        }
        catch (...)
        {
        }

        try
        {
            a->setHybridization(string_to_hybridization(atom.property<QString>(map["hybridization"])));
        }
        catch (...)
        {
        }
    }

    const auto bonds = SireMM::SelectorBond(mol, map);

    for (int i = 0; i < bonds.count(); ++i)
    {
        const auto bond = bonds(i);

        RDKit::Bond::BondType bondtype = RDKit::Bond::UNSPECIFIED;

        try
        {
            bondtype = string_to_bondtype(bond.property(map["bondtype"]).asAString());
        }
        catch (...)
        {
        }

        molecule.addBond(bond.atom0().index().value(),
                         bond.atom1().index().value(),
                         bondtype);
    }

    molecule.commitBatchEdit();

    RDKit::MolOps::sanitizeMol(molecule);

    return ROMOL_SPTR(new RDKit::ROMol(molecule));
}

Molecule rdkit_to_sire(const ROMOL_SPTR &mol, const PropertyMap &map)
{
    if (mol.get() == 0)
        return Molecule();

    QString molname;

    if (mol->hasProp("_Name"))
    {
        molname = QString::fromStdString(mol->getProp<std::string>("_Name"));
    }

    auto cg = Molecule().edit().rename(molname).add(SireMol::CGName("0"));
    auto res = cg.molecule().add(SireMol::ResNum(1));
    res.rename(SireMol::ResName("LIG"));

    int n = 0;

    QHash<const RDKit::Atom *, SireMol::AtomIdx> atom_to_idx;

    QHash<SireMol::Element, int> element_count;

    for (const auto &atom : mol->atoms())
    {
        atom_to_idx[atom] = n;
        n += 1;
        auto a = cg.add(SireMol::AtomNum(n));
        a.reparent(res.number());

        if (atom->hasProp("_Name"))
        {
            a.rename(SireMol::AtomName(
                QString::fromStdString(
                    atom->getProp<std::string>("_Name"))));
        }
        else
        {
            a.rename(SireMol::AtomName(
                QString("%1%2")
                    .arg(QString::fromStdString(atom->getSymbol()))
                    .arg(n)));
        }

        set_prop(a, "element", SireMol::Element(atom->getAtomicNum()), map);
        set_prop(a, "formal_charge", atom->getFormalCharge() * SireUnits::mod_electron, map);
        set_prop(a, "is_aromatic", qint64(atom->getIsAromatic()), map);
        set_prop(a, "mass", atom->getMass() * SireUnits::g_per_mol, map);
        set_prop(a, "isotope", qint64(atom->getIsotope()), map);
        set_prop(a, "chiral_tag", chiral_to_string(atom->getChiralTag()), map);
        set_prop(a, "hybridization", hybridization_to_string(atom->getHybridization()), map);
    }

    // we've built the main structure now - just need to add other
    // properties
    auto molecule = cg.molecule().commit().edit();

    if (map["connectivity"].hasSource())
    {
        auto connectivity = SireMol::Connectivity(molecule.info()).edit();

        const auto bondtype_as_double = map["bondtype_as_double"];
        const auto bondtype = map["bondtype"];

        for (const auto &bond : mol->bonds())
        {
            const auto atom0 = atom_to_idx[bond->getBeginAtom()];
            const auto atom1 = atom_to_idx[bond->getEndAtom()];

            connectivity.connect(atom0, atom1);

            const auto b = SireMol::BondID(atom0, atom1);

            if (bondtype_as_double.hasSource())
                connectivity.setProperty(b, bondtype_as_double.source(),
                                         SireBase::wrap(bond->getBondTypeAsDouble()));

            if (bondtype.hasSource())
                connectivity.setProperty(b, bondtype.source(),
                                         SireBase::wrap(bondtype_to_string(bond->getBondType())));
        }

        molecule.setProperty(map["connectivity"], connectivity.commit());
    }

    if (map["coordinates"].hasSource())
    {
        const int nconformers = mol->getNumConformers();

        if (nconformers > 0)
        {
            const auto &conformer = mol->getConformer(0);

            SireMol::AtomCoords coords(molecule.info());

            for (int i = 0; i < molecule.nAtoms(); ++i)
            {
                const auto pos = conformer.getAtomPos(i);
                coords.set(molecule.info().cgAtomIdx(SireMol::AtomIdx(i)),
                           SireMaths::Vector(pos.x, pos.y, pos.z));
            }

            molecule.setProperty(map["coordinates"].source(), coords);
        }
    }

    return molecule.commit();
}

SelectorMol rdkit_to_sire(const QList<ROMOL_SPTR> &mols, const PropertyMap &map)
{
    QList<Molecule> sire_mols;

    for (const auto &mol : mols)
    {
        sire_mols.append(rdkit_to_sire(mol, map));
    }

    return sire_mols;
}

QList<ROMOL_SPTR> sire_to_rdkit(const SelectorMol &mols, const PropertyMap &map)
{
    QList<ROMOL_SPTR> rdkit_mols;

    for (const auto &mol : mols)
    {
        rdkit_mols.append(sire_to_rdkit(mol, map));
    }

    return rdkit_mols;
}
