#ifndef SIRE_RDKIT_H
#define SIRE_RDKIT_H

#include "GraphMol/GraphMol.h"

#include "SireMol/selectormol.h"

#include <QList>

SireMol::SelectorMol rdkit_to_sire(const QList<RDKit::ROMol> &mols);

QList<RDKit::ROMol> sire_to_rdkit(const SireMol::SelectorMol &mols);

#endif
