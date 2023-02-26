
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "sire_openmm.h"

#include "Helpers/convertlist.hpp"

#include <QDebug>

namespace bp = boost::python;

using namespace SireOpenMM;

/** Thanks to this page for instructions on how to convert from SWIG to Boost
 *  https://wiki.python.org/moin/boost.python/HowTo
 */
struct PySwigObject
{
    PyObject_HEAD void *ptr;
    const char *desc;
};

void *extract_swig_wrapped_pointer(PyObject *obj)
{
    char thisStr[] = "this";

    // first we need to get the this attribute from the Python Object
    if (!PyObject_HasAttrString(obj, thisStr))
        return NULL;

    PyObject *thisAttr = PyObject_GetAttrString(obj, thisStr);
    if (thisAttr == NULL)
        return NULL;

    // This Python Object is a SWIG Wrapper and contains our pointer
    void *pointer = ((PySwigObject *)thisAttr)->ptr;
    Py_DECREF(thisAttr);
    return pointer;
}

BOOST_PYTHON_MODULE(_SireOpenMM)
{
    typedef SireMol::SelectorMol (*openmm_to_sire_function_type)(const OpenMM::System &, SireBase::PropertyMap const &);
    typedef void (*sire_to_openmm_function_type)(OpenMM::System &, const SireMol::SelectorMol &, SireBase::PropertyMap const &);

    sire_to_openmm_function_type sire_to_openmm_function_value(&sire_to_openmm);
    openmm_to_sire_function_type openmm_to_sire_function_value(&openmm_to_sire);

    bp::def("_sire_to_openmm_system",
            sire_to_openmm_function_value,
            (bp::arg("system"), bp::arg("mols"), bp::arg("map")),
            "Convert sire molecules to an OpenMM::System");

    bp::def("_openmm_system_to_sire",
            openmm_to_sire_function_value,
            (bp::arg("mols"), bp::arg("map")),
            "Convert an OpenMM::System to sire molecules");

    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::System>());
}
