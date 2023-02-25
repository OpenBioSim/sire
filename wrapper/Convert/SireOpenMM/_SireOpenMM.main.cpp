
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "sire_openmm.h"

#include "Helpers/convertlist.hpp"

namespace bp = boost::python;

using namespace SireOpenMM;

BOOST_PYTHON_MODULE(_SireOpenMM)
{
    typedef SireMol::SelectorMol (*openmm_to_sire_function_type)(const OpenMM::System &, SireBase::PropertyMap const &);
    typedef OpenMM::System (*sire_to_openmm_function_type)(const SireMol::SelectorMol &, SireBase::PropertyMap const &);

    sire_to_openmm_function_type sire_to_openmm_function_value(&sire_to_openmm);
    openmm_to_sire_function_type openmm_to_sire_function_value(&openmm_to_sire);

    bp::def("sire_to_openmm",
            sire_to_openmm_function_value,
            (bp::arg("mols"), bp::arg("map")),
            "Convert sire molecules to an OpenMM::System");

    bp::def("openmm_to_sire",
            openmm_to_sire_function_value,
            (bp::arg("mols"), bp::arg("map")),
            "Convert an OpenMM::System to sire molecules");
}
