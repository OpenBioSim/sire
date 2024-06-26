// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ChargeParameterName3D.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/countflops.h"

#include "SireBase/errors.h"

#include "SireBase/propertylist.h"

#include "SireBase/sparsematrix.hpp"

#include "SireError/errors.h"

#include "SireFF/errors.h"

#include "SireMaths/maths.h"

#include "SireMol/atomcoords.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "SireVol/cartesian.h"

#include "coulombpotential.h"

#include "switchingfunction.h"

#include <QDebug>

#include "coulombpotential.h"

SireMM::ChargeParameterName3D __copy__(const SireMM::ChargeParameterName3D &other){ return SireMM::ChargeParameterName3D(other); }

#include "Helpers/copy.hpp"

const char* pvt_get_name(const SireMM::ChargeParameterName3D&){ return "SireMM::ChargeParameterName3D";}

#include "Helpers/release_gil_policy.hpp"

void register_ChargeParameterName3D_class(){

    { //::SireMM::ChargeParameterName3D
        typedef bp::class_< SireMM::ChargeParameterName3D, bp::bases< SireMM::ChargeParameterName > > ChargeParameterName3D_exposer_t;
        ChargeParameterName3D_exposer_t ChargeParameterName3D_exposer = ChargeParameterName3D_exposer_t( "ChargeParameterName3D", "This class provides the default name of the properties\nthat contain the charge, LJ and 3D coordinates properties", bp::init< >("") );
        bp::scope ChargeParameterName3D_scope( ChargeParameterName3D_exposer );
        ChargeParameterName3D_exposer.def( "__copy__", &__copy__<SireMM::ChargeParameterName3D>);
        ChargeParameterName3D_exposer.def( "__deepcopy__", &__copy__<SireMM::ChargeParameterName3D>);
        ChargeParameterName3D_exposer.def( "clone", &__copy__<SireMM::ChargeParameterName3D>);
        ChargeParameterName3D_exposer.def( "__str__", &pvt_get_name);
        ChargeParameterName3D_exposer.def( "__repr__", &pvt_get_name);
    }

}
