// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License


#include "boost/python.hpp"

#include "Celsius.pypp.hpp"

#include "Fahrenheit.pypp.hpp"

#include "GeneralUnit.pypp.hpp"

#include "TempBase.pypp.hpp"

#include "Unit.pypp.hpp"

#include "_Units_free_functions.pypp.hpp"

namespace bp = boost::python;

#include "SireUnits_containers.h"

#include "SireUnits_registrars.h"

#include "SireUnits/temperature.h"

#include "sireunits_dimensions.h"

#include "generalunit.h"

#include "_Units_global_variables.pyman.hpp"

BOOST_PYTHON_MODULE(_Units){
    register_SireUnits_objects();

    register_SireUnits_containers();

    register_TempBase_class();

    register_Celsius_class();

    register_Unit_class();

    register_GeneralUnit_class();

    register_Fahrenheit_class();

    bp::implicitly_convertible< SireUnits::Dimension::TempBase, SireUnits::Dimension::Temperature >();

    bp::implicitly_convertible< double, SireUnits::Dimension::GeneralUnit >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::GeneralUnit >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Mass >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::MolarMass >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Length >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Time >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Charge >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::MolarCharge >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Temperature >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Angle >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Area >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Volume >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Velocity >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::AngularVelocity >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Acceleration >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::AngularAcceleration >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Energy >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::MolarEnergy >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Power >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::MolarPower >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Density >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::MolarDensity >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Force >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Pressure >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Capacitance >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Current >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::Potential >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::HarmonicBondConstant >();

    bp::implicitly_convertible< QString, SireUnits::Dimension::HarmonicAngleConstant >();

    register_SireUnits_dimensions();

    register_man_global_variables();

    register_free_functions();
}

