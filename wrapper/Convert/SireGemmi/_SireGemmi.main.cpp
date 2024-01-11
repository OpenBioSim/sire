
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "boost/python/converter/registry.hpp"

#include "sire_gemmi.h"

#include "Helpers/pyboost11.hpp"

namespace bp = boost::python;

using namespace SireGemmi;

#include <gemmi/model.hpp>

BOOST_PYTHON_MODULE(_SireGemmi)
{
    bp::def("sire_to_gemmi",
            &sire_to_gemmi,
            (bp::arg("mols"), bp::arg("map")),
            "Convert sire system to a gemmi structure");

    bp::def("gemmi_to_sire",
            &gemmi_to_sire,
            (bp::arg("mols"), bp::arg("map")),
            "Convert a gemmi Structure to a sire system");

    bp::def("_register_pdbx_loader",
            &register_pdbx_loader,
            "Internal function called once used to register PDBx support");

    pyboost11::converter<gemmi::Structure>();
}
