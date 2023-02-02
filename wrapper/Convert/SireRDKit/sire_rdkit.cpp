
#include "sire_rdkit.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"

#include <QDebug>

using RDKit::ROMol;
using SireMol::Molecule;
using SireMol::SelectorMol;

Molecule rdkit_to_sire(const ROMol &mol)
{
    auto cg = Molecule().edit().add(SireMol::CGName("0"));
    auto res = cg.molecule().add(SireMol::ResNum(1));
    res.rename(SireMol::ResName("LIG"));

    int n = 0;

    for (const auto &atom : mol.atoms())
    {
        n += 1;
        auto a = cg.add(SireMol::AtomNum(n));
        a.reparent(res.number());
        a.rename(SireMol::AtomName(QString::fromStdString(atom->getSymbol())));
    }

    return cg.molecule().commit();
}

SelectorMol rdkit_to_sire(const QList<ROMol> &mols)
{
    QList<Molecule> sire_mols;

    for (const auto &mol : mols)
    {
        sire_mols.append(rdkit_to_sire(mol));
    }

    return sire_mols;
}

QList<ROMol> sire_to_rdkit(const SelectorMol &mols)
{
    return QList<ROMol>();
}
