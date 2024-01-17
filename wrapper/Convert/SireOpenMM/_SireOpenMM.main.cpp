
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include <boost/tuple/tuple.hpp>

#include "sire_openmm.h"

#include "emle.h"

#include "lambdalever.h"

#include "openmmminimise.h"

#include "Helpers/convertdict.hpp"
#include "Helpers/convertlist.hpp"
#include "Helpers/tuples.hpp"

#include <QMap>
#include <QVector>

using namespace SireOpenMM;

using boost::python::register_tuple;

namespace bp = boost::python;

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

// Copy constructor for Python bindings.
SireOpenMM::EMLEEngine __copy__(const SireOpenMM::EMLEEngine &other){ return SireOpenMM::EMLEEngine(other); }

BOOST_PYTHON_MODULE(_SireOpenMM)
{
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
        "set_lambda", &LambdaLever::setLambda,
        (bp::arg("system"), bp::arg("lambda_value")),
        "Update the parameters in the passed context using this lambda lever "
        "so that the parameters represent the system at the specified "
        "lambda value");

    LambdaLever_exposer_t.def(
        "schedule", &LambdaLever::getSchedule,
        "Return the LambdaSchedule used to control the parameters by lambda");

    LambdaLever_exposer_t.def(
        "set_schedule", &LambdaLever::setSchedule,
        "Set the LambdaSchedule used to control the parameters by lambda");

    LambdaLever_exposer_t.def(
        "get_perturbable_molecule_maps", &LambdaLever::getPerturbableMoleculeMaps,
        "Return the perturbable molecule maps for all of the perturbable molecules");

    bp::def("_openmm_system_to_sire",
            &openmm_system_to_sire,
            (bp::arg("system"), bp::arg("map")),
            "Convert an OpenMM::System to a set of sire molecules.");

    bp::def("_sire_to_openmm_system",
            &sire_to_openmm_system,
            (bp::arg("system"), bp::arg("mols"), bp::arg("map")),
            "Convert sire molecules to an OpenMM::System");

    bp::def("_set_openmm_coordinates_and_velocities",
            &set_openmm_coordinates_and_velocities,
            (bp::arg("context"), bp::arg("coords_and_velocities")),
            "Set the coordinates and velocities in a context");

    typedef SireMol::SelectorMol (*extract_coordinates_function_type1)(
        const OpenMM::State &,
        const SireMol::SelectorMol &,
        const QHash<SireMol::MolNum, SireBase::PropertyMap> &,
        const SireBase::PropertyMap &);

    typedef SireMol::SelectorM<SireMol::Atom> (*extract_coordinates_function_type2)(
        const OpenMM::State &,
        const SireMol::SelectorM<SireMol::Atom> &,
        const QHash<SireMol::MolNum, SireBase::PropertyMap> &,
        const SireBase::PropertyMap &);

    extract_coordinates_function_type1 extract_coordinates_value1(&extract_coordinates);
    extract_coordinates_function_type2 extract_coordinates_value2(&extract_coordinates);

    bp::def("_openmm_extract_coordinates",
            extract_coordinates_value1,
            (bp::arg("state"), bp::arg("mols"), bp::arg("perturbable_maps"), bp::arg("map")),
            "Extract the coordinates from 'state' and copy then into the passed 'mols'");

    bp::def("_openmm_extract_coordinates",
            extract_coordinates_value2,
            (bp::arg("state"), bp::arg("mols"), bp::arg("perturbable_maps"), bp::arg("map")),
            "Extract the coordinates from 'state' and copy then into the passed 'atoms'");

    typedef SireMol::SelectorMol (*extract_coordinates_and_velocities_function_type1)(
        const OpenMM::State &,
        const SireMol::SelectorMol &,
        const QHash<SireMol::MolNum, SireBase::PropertyMap> &,
        const SireBase::PropertyMap &);

    typedef SireMol::SelectorM<SireMol::Atom> (*extract_coordinates_and_velocities_function_type2)(
        const OpenMM::State &,
        const SireMol::SelectorM<SireMol::Atom> &,
        const QHash<SireMol::MolNum, SireBase::PropertyMap> &,
        const SireBase::PropertyMap &);

    extract_coordinates_and_velocities_function_type1 extract_coordinates_and_velocities_value1(&extract_coordinates_and_velocities);
    extract_coordinates_and_velocities_function_type2 extract_coordinates_and_velocities_value2(&extract_coordinates_and_velocities);

    bp::def("_openmm_extract_coordinates_and_velocities",
            extract_coordinates_and_velocities_value1,
            (bp::arg("state"), bp::arg("mols"), bp::arg("perturbable_maps"), bp::arg("map")),
            "Extract the coordinates and velocities from 'state' and copy then into the passed 'mols'");

    bp::def("_openmm_extract_coordinates_and_velocities",
            extract_coordinates_and_velocities_value2,
            (bp::arg("state"), bp::arg("mols"), bp::arg("perturbable_maps"), bp::arg("map")),
            "Extract the coordinates and velocities from 'state' and copy then into the passed 'atoms'");

    bp::def("_openmm_extract_space",
            &extract_space,
            (bp::arg("state")),
            "Extract and return the space from 'state'");

    bp::def("_openmm_set_context_platform_property",
            &set_context_platform_property,
            (bp::arg("context"), bp::arg("key"), bp::arg("value")),
            "Set the Platform property for the passed context.");

    bp::def("_minimise_openmm_context",
            &minimise_openmm_context,
            (bp::arg("context"), bp::arg("tolerance") = 10, bp::arg("max_iterations") = -1),
            "Minimise the passed context");

    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::System>());
    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::Context>());
    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::State>());
    bp::converter::registry::insert(&extract_swig_wrapped_pointer, bp::type_id<OpenMM::Integrator>());

    // A tuple return type container for EMLECallback. (Energy, QM forces, MM forces)
    register_tuple<boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>>();

    // Dictionary for mapping link atoms to QM and MM2 atoms.
    register_dict<QMap<int, int>>();
    register_dict<QMap<int, double>>();
    register_dict<QMap<int, QVector<int>>>();

    // A tuple for passing link atom information to EMLEEngine.
    register_tuple<boost::tuple<QMap<int, int>, QMap<int, QVector<int>>>>();

    bp::class_<QMMMForce, bp::bases<SireBase::Property>, boost::noncopyable>("QMMMForce", bp::no_init);

    bp::class_<EMLEEngine, bp::bases<SireOpenMM::QMMMForce, SireBase::Property>>("EMLEEngine",
            bp::init<bp::object, SireUnits::Dimension::Length, int, double>(
                (
                    bp::arg("py_object"),
                    bp::arg("cutoff")=SireUnits::Dimension::Length(8.0),
                    bp::arg("neighbour_list_frequency")=20,
                    bp::arg("lambda")=1.0
                ),
                    "Constructor: An engine that can be used to enable electrostatic embedding"
                    "of machine learning potentials via emle-engine."
                )
            )
            .def("getLambda", &EMLEEngine::getLambda, "Get the lambda value")
            .def("setLambda", &EMLEEngine::setLambda, "Set the lambda value")
            .def("getCutoff", &EMLEEngine::getCutoff, "Get the cutoff value")
            .def("setCutoff", &EMLEEngine::setCutoff, "Set the cutoff value")
            .def("getNeighbourListFrequency", &EMLEEngine::getNeighbourListFrequency, "Get the neighbour list frequency")
            .def("setNeighbourListFrequency", &EMLEEngine::setNeighbourListFrequency, "Set the neighbour list frequency")
            .def("getAtoms", &EMLEEngine::getAtoms, "Get QM atom indices")
            .def("setAtoms", &EMLEEngine::setAtoms, "Set the QM atom indices")
            .def("getLinkAtoms", &EMLEEngine::getLinkAtoms, "Get the link atoms")
            .def("setLinkAtoms", &EMLEEngine::setLinkAtoms, "Set the link atoms")
            .def("getNumbers", &EMLEEngine::getNumbers, "Get QM atomic numbers")
            .def("setNumbers", &EMLEEngine::setNumbers, "Set the QM atomic numbers")
            .def("getCharges", &EMLEEngine::getCharges, "Get the atomic charges")
            .def("setCharges", &EMLEEngine::setCharges, "Set the atomic charges")
            .def("what", &EMLEEngine::what, "Call the callback")
            .def("typeName", &EMLEEngine::typeName, "Call the callback")
            .def("call", &EMLEEngine::call, "Call the callback")
            .def( "__copy__", &__copy__)
            .def( "__deepcopy__", &__copy__)
            .def( "clone", &__copy__);

    bp::class_<EMLECallback>("EMLECallback",
            bp::init<bp::object, QString>(
                "Constructor: A callback wrapper class to enable electrostatic embedding"
                "of machine learning potentials via emle-engine."
                )
            )
        .def("call", &EMLECallback::call, "Call the callback");
}
