
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
    typedef SireMol::SelectorMol (*openmm_system_to_sire_function_type)(const OpenMM::System &, SireBase::PropertyMap const &);
    typedef CoordsAndVelocities (*sire_to_openmm_system_function_type)(OpenMM::System &, const SireMol::SelectorMol &, SireBase::PropertyMap const &);
    typedef void (*set_openmm_coordinates_and_velocities_function_type)(OpenMM::Context &, const CoordsAndVelocities &);

    openmm_system_to_sire_function_type openmm_system_to_sire_function_value(&openmm_system_to_sire);
    sire_to_openmm_system_function_type sire_to_openmm_system_function_value(&sire_to_openmm_system);
    set_openmm_coordinates_and_velocities_function_type set_openmm_coordinates_and_velocities_function_value(&set_openmm_coordinates_and_velocities);

    bp::class_<CoordsAndVelocities> CoordsAndVelocities_exposer_t("CoordsAndVelocities",
                                                                  "Internal class used to hold OpenMM coordinates and velocities data");

    bp::def("_openmm_system_to_sire",
            openmm_system_to_sire_function_value,
            (bp::arg("system"), bp::arg("map")),
            "Convert an OpenMM::System to a set of sire molecules.");

    bp::def("_sire_to_openmm_system",
            sire_to_openmm_system_function_value,
            (bp::arg("system"), bp::arg("mols"), bp::arg("map")),
            "Convert sire molecules to an OpenMM::System");

    bp::def("_set_openmm_coordinates_and_velocities",
            set_openmm_coordinates_and_velocities_function_value,
            (bp::arg("context"), bp::arg("coords_and_velocities")),
            "Set the coordinates and velocities in a context");

    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::System>());
    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::Context>());
}
