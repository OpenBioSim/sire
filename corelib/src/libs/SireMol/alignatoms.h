#ifndef SIREMOL_ALIGNATOMS_H
#define SIREMOL_ALIGNATOMS_H

#include "selector.hpp"
#include "selectorm.hpp"
#include "atom.h"

#include "SireMaths/align.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    SireMaths::Transform get_atom_alignment(const SelectorM<Atom> &atoms0,
                                            const SelectorM<Atom> &atoms1,
                                            bool fit = true);

    SireMaths::Transform get_atom_alignment(const SelectorM<Atom> &atoms0,
                                            const SelectorM<Atom> &atoms1,
                                            const PropertyMap &map,
                                            bool fit = true);

    SireMaths::Transform get_atom_alignment(const SelectorM<Atom> &atoms0,
                                            const SelectorM<Atom> &atoms1,
                                            const PropertyMap &map0,
                                            const PropertyMap &map1,
                                            bool fit = true);

}

SIRE_EXPOSE_FUNCTION(SireMol::align_atoms)

SIRE_END_HEADER

#endif
