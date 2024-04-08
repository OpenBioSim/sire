###############################################
#
# This file contains special code to help
# with the wrapping of SireUnits classes
#
#

import re


def findGlobals():
    # read in the information about this module
    lines = open("module_info", "r").readlines()
    root = lines[2].split()[1]
    sourcedir = lines[1].split()[1]

    lines = open("%s/%s/units.h" % (root, sourcedir), "r").readlines()

    FILE = open("_Units_global_variables.pyman.hpp", "w")

    print("#ifndef _Units_global_variables_hpp", file=FILE)
    print("#define _Units_global_variables_hpp", file=FILE)
    print("\nvoid register_man_global_variables();\n", file=FILE)
    print("#endif", file=FILE)

    FILE.close()

    FILE = open("_Units_global_variables.pyman.cpp", "w")

    print('\n#include "_Units_global_variables.pyman.hpp"', file=FILE)
    print("#include <boost/python.hpp>", file=FILE)
    print('#include "SireUnits/units.h"', file=FILE)
    print('#include "SireUnits/temperature.h"', file=FILE)
    print("\nusing namespace boost::python;", file=FILE)
    print("using namespace SireUnits;", file=FILE)
    print("using namespace SireUnits::Dimension;\n", file=FILE)

    print("void register_man_global_variables()", file=FILE)
    print("{", file=FILE)

    for line in lines:
        match = re.search(r"const Dimension::([\w\d\-<,>]+)\s+(\w+)", line)

        if match:
            name = match.group(2)
            print('    scope().attr("%s") = %s;\n' % (name, name), file=FILE)
        else:
            match = re.search(r"const GeneralUnit\s+(\w+)", line)

            if match:
                name = match.group(1)
                print('    scope().attr("%s") = %s;\n' % (name, name), file=FILE)

    # add Celsius and Fahrenheit manually
    print('    scope().attr("celsius") = celsius;\n', file=FILE)
    print('    scope().attr("fahrenheit") = fahrenheit;\n', file=FILE)

    print("}\n", file=FILE)


def fix_GeneralUnit(c):
    c.add_registration_code("def( bp::other<double>() + bp::self )")
    c.add_registration_code("def( bp::other<double>() - bp::self )")
    c.add_registration_code("def( bp::other<double>() * bp::self )")
    c.add_registration_code("def( bp::other<double>() / bp::self )")

    c.add_declaration_code('#include "generalunit.h"')


def fix_Temperature(c):
    c.add_registration_code("def( bp::other<double>() + bp::self )")
    c.add_registration_code("def( bp::other<double>() - bp::self )")
    c.add_registration_code("def( bp::other<double>() * bp::self )")
    c.add_registration_code("def( bp::other<double>() / bp::self )")

    c.add_declaration_code('#include "generalunit.h"')

    c.add_registration_code(
        "def( bp::self + bp::other<SireUnits::Dimension::GeneralUnit>() )"
    )
    c.add_registration_code(
        "def( bp::self - bp::other<SireUnits::Dimension::GeneralUnit>() )"
    )
    c.add_registration_code(
        "def( bp::self * bp::other<SireUnits::Dimension::GeneralUnit>() )"
    )
    c.add_registration_code(
        "def( bp::self / bp::other<SireUnits::Dimension::GeneralUnit>() )"
    )


def fix_TempBase(c):
    c.add_declaration_code('#include "generalunit.h"')


def fixMB(mb):
    mb.add_declaration_code('#include "SireUnits/temperature.h"')
    mb.add_declaration_code('#include "sireunits_dimensions.h"')
    mb.add_declaration_code('#include "generalunit.h"')
    mb.add_declaration_code('#include "_Units_global_variables.pyman.hpp"')

    mb.add_registration_code("register_SireUnits_dimensions();")
    mb.add_registration_code("register_man_global_variables();")

    # add all of the global physical constants to the module
    findGlobals()


special_code = {
    "SireUnits::Dimension::GeneralUnit": fix_GeneralUnit,
    "SireUnits::Celsius": fix_Temperature,
    "SireUnits::Fahrenheit": fix_Temperature,
    "SireUnits::Dimension::TempBase": fix_TempBase,
}

implicitly_convertible = [
    ("SireUnits::Dimension::TempBase", "SireUnits::Dimension::Temperature"),
    ("double", "SireUnits::Dimension::GeneralUnit"),
    ("QString", "SireUnits::Dimension::GeneralUnit"),
    ("QString", "SireUnits::Dimension::Mass"),
    ("QString", "SireUnits::Dimension::MolarMass"),
    ("QString", "SireUnits::Dimension::Length"),
    ("QString", "SireUnits::Dimension::Time"),
    ("QString", "SireUnits::Dimension::Charge"),
    ("QString", "SireUnits::Dimension::MolarCharge"),
    ("QString", "SireUnits::Dimension::Temperature"),
    ("QString", "SireUnits::Dimension::Angle"),
    ("QString", "SireUnits::Dimension::Area"),
    ("QString", "SireUnits::Dimension::Volume"),
    ("QString", "SireUnits::Dimension::Velocity"),
    ("QString", "SireUnits::Dimension::AngularVelocity"),
    ("QString", "SireUnits::Dimension::Acceleration"),
    ("QString", "SireUnits::Dimension::AngularAcceleration"),
    ("QString", "SireUnits::Dimension::Energy"),
    ("QString", "SireUnits::Dimension::MolarEnergy"),
    ("QString", "SireUnits::Dimension::Power"),
    ("QString", "SireUnits::Dimension::MolarPower"),
    ("QString", "SireUnits::Dimension::Density"),
    ("QString", "SireUnits::Dimension::MolarDensity"),
    ("QString", "SireUnits::Dimension::Force"),
    ("QString", "SireUnits::Dimension::Pressure"),
    ("QString", "SireUnits::Dimension::Capacitance"),
    ("QString", "SireUnits::Dimension::Current"),
    ("QString", "SireUnits::Dimension::Potential"),
    ("QString", "SireUnits::Dimension::HarmonicBondConstant"),
    ("QString", "SireUnits::Dimension::HarmonicAngleConstant"),
]
