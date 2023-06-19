
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "sire_rdkit.h"

#include "Helpers/convertlist.hpp"

namespace bp = boost::python;

using namespace SireRDKit;

BOOST_PYTHON_MODULE(_SireRDKit)
{
        typedef SireMol::SelectorMol (*rdkit_to_sire_function_type)(QList<RDKit::ROMOL_SPTR> const &, SireBase::PropertyMap const &);
        typedef QList<RDKit::ROMOL_SPTR> (*sire_to_rdkit_function_type)(SireMol::SelectorMol const &, SireBase::PropertyMap const &);

        typedef QStringList (*rdkit_to_smiles_function_type)(QList<RDKit::ROMOL_SPTR> const &, SireBase::PropertyMap const &);
        typedef QString (*rdkit_to_smiles_function_type2)(RDKit::ROMOL_SPTR const &, SireBase::PropertyMap const &);

        typedef QStringList (*rdkit_to_smarts_function_type)(QList<RDKit::ROMOL_SPTR> const &, SireBase::PropertyMap const &);
        typedef QString (*rdkit_to_smarts_function_type2)(RDKit::ROMOL_SPTR const &, SireBase::PropertyMap const &);

        typedef QList<RDKit::ROMOL_SPTR> (*rdkit_remove_hydrogens_function_type)(QList<RDKit::ROMOL_SPTR> const &, SireBase::PropertyMap const &);
        typedef RDKit::ROMOL_SPTR (*rdkit_remove_hydrogens_function_type2)(RDKit::ROMOL_SPTR const &, SireBase::PropertyMap const &);

        typedef QList<RDKit::ROMOL_SPTR> (*smiles_to_rdkit_function_type)(QStringList const &, QStringList const &, SireBase::PropertyMap const &);
        typedef RDKit::ROMOL_SPTR (*smiles_to_rdkit_function_type2)(QString const &, QString const &, SireBase::PropertyMap const &);

        typedef QList<RDKit::ROMOL_SPTR> (*smarts_to_rdkit_function_type)(QStringList const &, QStringList const &, SireBase::PropertyMap const &);
        typedef RDKit::ROMOL_SPTR (*smarts_to_rdkit_function_type2)(QString const &, QString const &, SireBase::PropertyMap const &);

        rdkit_to_sire_function_type rdkit_to_sire_function_value(&rdkit_to_sire);
        sire_to_rdkit_function_type sire_to_rdkit_function_value(&sire_to_rdkit);

        rdkit_to_smiles_function_type rdkit_to_smiles_function_value(&rdkit_to_smiles);
        rdkit_to_smiles_function_type2 rdkit_to_smiles_function_value2(&rdkit_to_smiles);

        rdkit_to_smarts_function_type rdkit_to_smarts_function_value(&rdkit_to_smarts);
        rdkit_to_smarts_function_type2 rdkit_to_smarts_function_value2(&rdkit_to_smarts);

        rdkit_remove_hydrogens_function_type rdkit_remove_hydrogens_function_value(&rdkit_remove_hydrogens);
        rdkit_remove_hydrogens_function_type2 rdkit_remove_hydrogens_function_value2(&rdkit_remove_hydrogens);

        smiles_to_rdkit_function_type smiles_to_rdkit_function_value(&smiles_to_rdkit);
        smiles_to_rdkit_function_type2 smiles_to_rdkit_function_value2(&smiles_to_rdkit);

        smarts_to_rdkit_function_type smarts_to_rdkit_function_value(&smarts_to_rdkit);
        smarts_to_rdkit_function_type2 smarts_to_rdkit_function_value2(&smarts_to_rdkit);

        bp::def("sire_to_rdkit",
                sire_to_rdkit_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Convert sire molecules to a list of RDKit molecules");

        bp::def("rdkit_to_sire",
                rdkit_to_sire_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Convert a list of RDKit molecules to sire molecules");

        bp::def("rdkit_to_smiles",
                rdkit_to_smiles_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Convert a list of RDKit molecules to smiles strings");

        bp::def("rdkit_to_smarts",
                rdkit_to_smarts_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Convert a list of RDKit molecules to smarts strings");

        bp::def("rdkit_remove_hydrogens",
                rdkit_remove_hydrogens_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Remove unnecessary hydrogens from the passed list of RDKit molecules");

        bp::def("rdkit_to_smiles",
                rdkit_to_smiles_function_value2,
                (bp::arg("mol"), bp::arg("map")),
                "Convert a RDKit molecule to a smiles string");

        bp::def("rdkit_to_smarts",
                rdkit_to_smarts_function_value2,
                (bp::arg("mol"), bp::arg("map")),
                "Convert a RDKit molecule to a smarts string");

        bp::def("rdkit_remove_hydrogens",
                rdkit_remove_hydrogens_function_value,
                (bp::arg("mols"), bp::arg("map")),
                "Remove unnecessary hydrogens from the passed list of RDKit molecules");

        bp::def("rdkit_remove_hydrogens",
                rdkit_remove_hydrogens_function_value2,
                (bp::arg("mol"), bp::arg("map")),
                "Remove unnecessary hydrogens from the passed RDKit molecule");

        bp::def("smiles_to_rdkit",
                smiles_to_rdkit_function_value,
                (bp::arg("smiles"), bp::arg("labels"), bp::arg("map")),
                "Convert a list of smiles strings to RDKit molecules");

        bp::def("smiles_to_rdkit",
                smiles_to_rdkit_function_value2,
                (bp::arg("smiles"), bp::arg("label"), bp::arg("map")),
                "Convert a smiles string to an RDKit molecule");

        bp::def("smarts_to_rdkit",
                smarts_to_rdkit_function_value,
                (bp::arg("smarts"), bp::arg("labels"), bp::arg("map")),
                "Convert a list of smarts strings to RDKit molecules");

        bp::def("smarts_to_rdkit",
                smarts_to_rdkit_function_value2,
                (bp::arg("smarts"), bp::arg("label"), bp::arg("map")),
                "Convert a smarts string to an RDKit molecule");

        bp::def("_register_smarts_search",
                &register_smarts_search,
                "Internal function called once used to register smarts searching");

        register_list<QList<RDKit::ROMOL_SPTR>>();
}
