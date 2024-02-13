#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireOpenMM_properties.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "qmmm.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "qmmm.h"
void register_SireOpenMM_properties()
{
    register_property_container< SireOpenMM::QMEnginePtr, SireOpenMM::QMEngine >();
    register_property_container< SireOpenMM::QMEnginePtr, SireOpenMM::QMEngine >();
}
