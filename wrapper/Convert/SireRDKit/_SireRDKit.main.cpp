
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "sire_rdkit.h"

#include "Helpers/convertlist.hpp"

namespace bp = boost::python;

BOOST_PYTHON_MODULE(_SireRDKit)
{
        typedef SireMol::SelectorMol (*rdkit_to_sire_function_type)(QList<RDKit::ROMOL_SPTR> const &, SireBase::PropertyMap const &);
        typedef QList<RDKit::ROMOL_SPTR> (*sire_to_rdkit_function_type)(SireMol::SelectorMol const &, SireBase::PropertyMap const &);

        rdkit_to_sire_function_type rdkit_to_sire_function_value(&rdkit_to_sire);
        sire_to_rdkit_function_type sire_to_rdkit_function_value(&sire_to_rdkit);

        bp::def("sire_to_rdkit",
                sire_to_rdkit_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Convert sire molecules to a list of RDKit molecules");

        bp::def("rdkit_to_sire",
                rdkit_to_sire_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Convert a list of RDKit molecules to sire molecules");

        register_list<QList<RDKit::ROMOL_SPTR>>();
}
