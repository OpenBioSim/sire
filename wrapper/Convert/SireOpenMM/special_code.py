###############################################
#
# This file contains special code to help
# with the wrapping of SireOpenMM classes
#
#


def fixMB(mb):
    mb.add_declaration_code('#include "./register_extras.h"')
    mb.add_registration_code("SireOpenMM::register_extras();")


def fix_OpenMM_include(c):
    c.add_declaration_code('#include "OpenMM.h"')


special_code = {
    "std::vector<OpenMM::vec3>": fix_OpenMM_include,
}
