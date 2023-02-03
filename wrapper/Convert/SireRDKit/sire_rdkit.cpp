
#include "sire_rdkit.h"

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

#include "SireMaths/vector.h"

#include "SireBase/propertylist.h"
#include "SireBase/stringproperty.h"

#include "SireUnits/units.h"

#include <QDebug>

using RDKit::ROMol;
using SireBase::PropertyMap;
using SireMol::Molecule;
using SireMol::SelectorMol;

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

template <class T, class V>
void set_prop(T &atom, const QString &key, const V &value,
              const PropertyMap &map)
{
    const auto name = map[key];

    if (name.hasSource())
        atom.setProperty(name.source(), value);
}

ROMol sire_to_rdkit(const Molecule &mol, const PropertyMap &map)
{
    RDKit::RWMol molecule;
    molecule.beginBatchEdit();

    const auto atoms = mol.atoms();

    for (int i = 0; i < atoms.count(); ++i)
    {
        const auto atom = atoms(i);

        molecule.addAtom(true);
        auto a = molecule.getActiveAtom();

        a->setAtomicNum(atom.property<SireMol::Element>(map["element"]).nProtons());
    }

    molecule.commitBatchEdit();

    return molecule;
}

Molecule rdkit_to_sire(const ROMol &mol, const PropertyMap &map)
{
    auto cg = Molecule().edit().add(SireMol::CGName("0"));
    auto res = cg.molecule().add(SireMol::ResNum(1));
    res.rename(SireMol::ResName("LIG"));

    int n = 0;

    QHash<const RDKit::Atom *, SireMol::AtomIdx> atom_to_idx;

    QHash<SireMol::Element, int> element_count;

    for (const auto &atom : mol.atoms())
    {
        atom_to_idx[atom] = n;
        n += 1;
        auto a = cg.add(SireMol::AtomNum(n));
        a.reparent(res.number());

        a.rename(SireMol::AtomName(
            QString("%1%2")
                .arg(QString::fromStdString(atom->getSymbol()))
                .arg(n)));

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

        for (const auto &bond : mol.bonds())
        {
            const auto atom0 = atom_to_idx[bond->getBeginAtom()];
            const auto atom1 = atom_to_idx[bond->getEndAtom()];

            connectivity.connect(atom0, atom1);

            const auto b = SireMol::BondID(atom0, atom1);

            if (bondtype_as_double.hasSource())
                connectivity.setProperty(b, bondtype_as_double.source(),
                                         SireBase::wrap(bond->getBondTypeAsDouble()));
        }

        molecule.setProperty(map["connectivity"], connectivity.commit());
    }

    if (map["coordinates"].hasSource())
    {
        const int nconformers = mol.getNumConformers();

        if (nconformers > 0)
        {
            const auto &conformer = mol.getConformer(0);

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

SelectorMol rdkit_to_sire(const QList<ROMol> &mols, const PropertyMap &map)
{
    QList<Molecule> sire_mols;

    for (const auto &mol : mols)
    {
        sire_mols.append(rdkit_to_sire(mol, map));
    }

    return sire_mols;
}

QList<ROMol> sire_to_rdkit(const SelectorMol &mols, const PropertyMap &map)
{
    QList<ROMol> rdkit_mols;

    for (const auto &mol : mols)
    {
        rdkit_mols.append(sire_to_rdkit(mol, map));
    }

    return rdkit_mols;
}
