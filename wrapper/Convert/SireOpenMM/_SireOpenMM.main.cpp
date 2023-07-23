
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "sire_openmm.h"

#include "lambdalever.h"

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
    typedef OpenMMMetaData (*sire_to_openmm_system_function_type)(OpenMM::System &, const SireMol::SelectorMol &, SireBase::PropertyMap const &);
    typedef void (*set_openmm_coordinates_and_velocities_function_type)(OpenMM::Context &, const OpenMMMetaData &);
    typedef SireMol::SelectorMol (*extract_coordinates_function_type)(OpenMM::State const &, SireMol::SelectorMol const &, SireBase::PropertyMap const &);
    typedef SireMol::SelectorMol (*extract_coordinates_function_type)(OpenMM::State const &, SireMol::SelectorMol const &, SireBase::PropertyMap const &);
    typedef SireMol::SelectorMol (*extract_coordinates_and_velocities_function_type)(OpenMM::State const &, SireMol::SelectorMol const &, SireBase::PropertyMap const &);
    typedef SireVol::SpacePtr (*extract_space_function_type)(OpenMM::State const &);

    openmm_system_to_sire_function_type openmm_system_to_sire_function_value(&openmm_system_to_sire);
    sire_to_openmm_system_function_type sire_to_openmm_system_function_value(&sire_to_openmm_system);
    set_openmm_coordinates_and_velocities_function_type set_openmm_coordinates_and_velocities_function_value(&set_openmm_coordinates_and_velocities);
    extract_coordinates_function_type extract_coordinates_function_value(&extract_coordinates);
    extract_coordinates_and_velocities_function_type extract_coordinates_and_velocities_function_value(&extract_coordinates_and_velocities);
    extract_space_function_type extract_space_function_value(&extract_space);

    bp::class_<OpenMMMetaData> OpenMMMetaData_exposer_t("OpenMMMetaData",
                                                        "Internal class used to hold OpenMM coordinates and velocities data");

    OpenMMMetaData_exposer_t.def(
        "index", &OpenMMMetaData::index,
        "Return the index used to locate atoms in the OpenMM system");

    OpenMMMetaData_exposer_t.def(
        "lambdaLever", &OpenMMMetaData::lambdaLever,
        "Return the lambda lever used to update the parameters in the "
        "OpenMM system according to lambda");

    bp::class_<LambdaLever, bp::bases<SireBase::Property>> LambdaLever_exposer_t(
        "LambdaLever",
        "A lever that can be used to change the parameters in an OpenMM system "
        "based on a lambda value (or collection of lambda values)");

    LambdaLever_exposer_t.def(
        "set_lambda", &LambdaLever::set_lambda,
        (bp::arg("system"), bp::arg("lambda_value")),
        "Update the parameters in the passed context using this lambda lever "
        "so that the parameters represent the system at the specified "
        "lambda value");

    LambdaLever_exposer_t.def(
        "lambda", &LambdaLever::lambda,
        "Return the symbol used for the lambda value in each stage");

    LambdaLever_exposer_t.def(
        "initial", &LambdaLever::initial,
        "Return the symbol used for the initial (lambda=0) parameter in each stage");

    LambdaLever_exposer_t.def(
        "final", &LambdaLever::final,
        "Return the symbol used for the final (lambda=1) parameter in each stage");

    LambdaLever_exposer_t.def(
        "add_lever", &LambdaLever::add_lever,
        (bp::arg("lever")),
        "Add a new lever named 'lever'");

    LambdaLever_exposer_t.def(
        "add_levers", &LambdaLever::add_levers,
        (bp::arg("levers")),
        "Add several new levers, whose names are in 'levers'");

    LambdaLever_exposer_t.def(
        "num_levers", &LambdaLever::num_levers,
        "Return the number of lambda levers");

    LambdaLever_exposer_t.def(
        "get_levers", &LambdaLever::get_levers,
        "Return the names of all of the lambda levers");

    LambdaLever_exposer_t.def(
        "num_stages", &LambdaLever::num_stages,
        "Return the number of lambda lever stages");

    LambdaLever_exposer_t.def(
        "get_stages", &LambdaLever::get_stages,
        "Return the names of all of the stages, in the "
        "order in which they will be applied");

    LambdaLever_exposer_t.def(
        "get_stage", &LambdaLever::get_stage,
        (bp::arg("lambda_value")),
        "Return name of the stage that will be applied for the specified "
        "global value of lambda");

    LambdaLever_exposer_t.def(
        "get_lambda_in_stage", &LambdaLever::get_lambda_in_stage,
        (bp::arg("lambda_value")),
        "Return the stage-local value of lambda that will be applied within "
        "the associated stage for the specified global value of lambda");

    LambdaLever_exposer_t.def(
        "add_stage", &LambdaLever::add_stage,
        (bp::arg("stage"), bp::arg("equation")),
        "Add a new stage with the specified name, "
        "and with the passed equation used as default for all "
        "parameters for that stage");

    LambdaLever_exposer_t.def(
        "clear", &LambdaLever::clear,
        "Remove all stages from this lambda lever");

    LambdaLever_exposer_t.def(
        "set_equation", &LambdaLever::set_equation,
        (bp::arg("stage"), bp::arg("lever"), bp::arg("equation")),
        "Set the equation used for the specified lever for the specified stage.");

    LambdaLever_exposer_t.def(
        "set_default_equation", &LambdaLever::set_default_equation,
        (bp::arg("stage"), bp::arg("equation")),
        "Set the default equation used for the levers in a stage for which "
        "a specified equation has not been supplied");

    LambdaLever_exposer_t.def(
        "get_equation", &LambdaLever::get_equation,
        (bp::arg("stage"), bp::arg("lever")),
        "Return the equation used to set the parameters for the "
        "specified lever in the specified stage");

    LambdaLever_exposer_t.def(
        "get_lever_values", &LambdaLever::get_lever_values,
        (bp::arg("lambda_values"), bp::arg("initial_value"), bp::arg("final_value")),
        "Return the values of all of the levers for the specified global "
        "lambda values, assuming a parameter that has the specified initial "
        "and final values");

    LambdaLever_exposer_t.def(
        "get_lever_stages", &LambdaLever::get_lever_stages,
        (bp::arg("lambda_values")),
        "Return the name of the stage associated with each of the passed "
        "global lambda values");

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

    bp::def("_openmm_extract_coordinates",
            extract_coordinates_function_value,
            (bp::arg("state"), bp::arg("mols"), bp::arg("map")),
            "Extract the coordinates from 'state' and copy then into the passed 'mols'");

    bp::def("_openmm_extract_coordinates_and_velocities",
            extract_coordinates_and_velocities_function_value,
            (bp::arg("state"), bp::arg("mols"), bp::arg("map")),
            "Extract the coordinates and velocities from 'state' and copy then into the passed 'mols'");

    bp::def("_openmm_extract_space",
            extract_space_function_value,
            (bp::arg("state")),
            "Extract and return the space from 'state'");

    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::System>());
    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::Context>());
    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::State>());
    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::Integrator>());
}
