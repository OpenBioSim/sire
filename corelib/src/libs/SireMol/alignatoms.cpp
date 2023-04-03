
#include "alignatoms.h"

using namespace SireMaths;
using namespace SireMol;

/** Return the transform needed to align atoms1 on top of atom0,
 *  assuming that the atoms are mapped by index in the two sets.
 *  You need to rearrange one or the other set if you want
 *  to use a different mapping.
 */
SIREMOL_EXPORT Transform SireMol::get_atom_alignment(
    const SelectorM<Atom> &atoms0,
    const SelectorM<Atom> &atoms1,
    const PropertyMap &map0,
    const PropertyMap &map1,
    bool fit)
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return Transform();

    // get the coordinates of these atoms
    QVector<Vector> coords0 = atoms0.property<Vector>(map0["coordinates"]);
    QVector<Vector> coords1 = atoms1.property<Vector>(map1["coordinates"]);

    // must make sure that the two sets of coordinates have the
    // same size - truncate to matching sizes
    if (coords0.count() != coords1.count())
    {
        const auto min_size = std::mid(coords0.count(), coords1.count());

        coords0.resize(min_size);
        coords1.resize(min_size);
    }

    return SireMaths::getAlignment(coords0, coords1, fit);
}
