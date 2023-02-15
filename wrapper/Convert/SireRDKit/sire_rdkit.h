#ifndef SIRE_RDKIT_H
#define SIRE_RDKIT_H

#include "GraphMol/GraphMol.h"

#include "SireMol/selectormol.h"

#include "SireBase/propertymap.h"

#include <QList>

SireMol::SelectorMol rdkit_to_sire(const QList<RDKit::ROMOL_SPTR> &mols,
                                   const SireBase::PropertyMap &map);

QList<RDKit::ROMOL_SPTR> sire_to_rdkit(const SireMol::SelectorMol &mols,
                                       const SireBase::PropertyMap &map);

QStringList rdkit_to_smiles(const QList<RDKit::ROMOL_SPTR> &mols,
                            bool ignore_errors = true);

QList<RDKit::ROMOL_SPTR> rdkit_remove_hydrogens(const QList<RDKit::ROMOL_SPTR> &mols,
                                                bool ignore_errors = true);

#endif
