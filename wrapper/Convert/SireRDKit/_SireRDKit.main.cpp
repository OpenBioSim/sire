
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "SireMol/selectormol.h"

#include <QDebug>

namespace bp = boost::python;

QString sire_to_rdkit(const SireMol::SelectorMol &mols)
{
    return mols.toString();
}

SireMol::SelectorMol rdkit_to_sire(const QString &mols)
{
    qDebug() << mols;
    return SireMol::SelectorMol();
}

BOOST_PYTHON_MODULE(_SireRDKit)
{
    bp::def("sire_to_rdkit", &sire_to_rdkit);
    bp::def("rdkit_to_sire", &rdkit_to_sire);
}
