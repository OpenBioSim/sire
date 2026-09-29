
// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"

#include "sire_gemmi.h"

#include <gemmi/version.hpp>

namespace bp = boost::python;

using namespace SireGemmi;

bp::object sire_to_gemmi_state_py(const SireSystem::System &system,
                                  const SireBase::PropertyMap &map)
{
    const auto state = sire_to_gemmi_state(system, map);

    return bp::object(bp::handle<>(
        PyBytes_FromStringAndSize(state.data(), state.size())));
}

SireSystem::System gemmi_state_to_sire_py(const bp::object &state,
                                          const SireBase::PropertyMap &map)
{
    char *data = nullptr;
    Py_ssize_t size = 0;

    if (PyBytes_AsStringAndSize(state.ptr(), &data, &size) != 0)
        bp::throw_error_already_set();

    return gemmi_state_to_sire(std::string(data, size), map);
}

BOOST_PYTHON_MODULE(_SireGemmi)
{
    bp::def("_sire_to_gemmi_state",
            &sire_to_gemmi_state_py,
            (bp::arg("mols"), bp::arg("map")),
            "Convert sire system to the pickle state of a gemmi Structure");

    bp::def("_gemmi_state_to_sire",
            &gemmi_state_to_sire_py,
            (bp::arg("state"), bp::arg("map")),
            "Convert the pickle state of a gemmi Structure to a sire system");

    bp::scope().attr("_gemmi_version") = GEMMI_VERSION;

    bp::def("_register_pdbx_loader",
            &register_pdbx_loader,
            "Internal function called once used to register PDBx support");
}
