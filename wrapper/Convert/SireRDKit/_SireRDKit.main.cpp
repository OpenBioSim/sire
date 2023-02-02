
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "sire_rdkit.h"

#include "Helpers/convertlist.hpp"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_SireRDKit)
{
    bp::def("sire_to_rdkit", &sire_to_rdkit);
    bp::def("rdkit_to_sire", &rdkit_to_sire);

    register_list<QList<RDKit::ROMol>>();
}
