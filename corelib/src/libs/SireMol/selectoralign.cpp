
#include "SireMol/core.h"
#include "SireMol/selectorm.hpp"

using namespace SireMaths;
using namespace SireMol;

Transform SireMol::detail::_getAlignment(const SelectorM<Atom> &atoms0,
                                         const SelectorM<Atom> &atoms1,
                                         const PropertyMap &map0,
                                         const PropertyMap &map1)
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

    return SireMaths::getAlignment(coords0, coords1, true);
}
