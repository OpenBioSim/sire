
#include "register_extras.h"

#include "boost/python.hpp"

#include <OpenMM.h>

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

namespace SireOpenMM
{
    void register_extras()
    {
        bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::System>());
        bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::Context>());
        bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::State>());
        bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::Integrator>());
    }
}
