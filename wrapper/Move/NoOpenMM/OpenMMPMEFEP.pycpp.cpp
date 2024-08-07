// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "OpenMMPMEFEP.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/forcetable.h"

#include "SireIO/amber.h"

#include "SireMM/atomljs.h"

#include "SireMM/internalff.h"

#include "SireMM/internalperturbation.h"

#include "SireMaths/constants.h"

#include "SireMaths/rangenerator.h"

#include "SireMaths/vector.h"

#include "SireMol/amberparameters.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atommasses.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/core.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/partialmolecule.h"

#include "SireMol/perturbation.h"

#include "SireMove/flexibility.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/convert.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "ensemble.h"

#include "openmmfrenergyst.h"

#include <QDebug>

#include <QTime>

#include <iostream>

#include "openmmfrenergyst.h"

SireMove::OpenMMPMEFEP __copy__(const SireMove::OpenMMPMEFEP &other){ return SireMove::OpenMMPMEFEP(other); }

const char* pvt_get_name(const SireMove::OpenMMPMEFEP&){ return "SireMove::OpenMMPMEFEP";}

void register_OpenMMPMEFEP_class(){

    { //::SireMove::OpenMMPMEFEP
        typedef bp::class_< SireMove::OpenMMPMEFEP > OpenMMPMEFEP_exposer_t;
        OpenMMPMEFEP_exposer_t OpenMMPMEFEP_exposer = OpenMMPMEFEP_exposer_t( "OpenMMPMEFEP", bp::init< >() );
        bp::scope OpenMMPMEFEP_scope( OpenMMPMEFEP_exposer );
        { //::SireMove::OpenMMPMEFEP::typeName

            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::OpenMMPMEFEP::typeName );

            OpenMMPMEFEP_exposer.def(
                "typeName"
                , typeName_function_value );

        }
        OpenMMPMEFEP_exposer.staticmethod( "typeName" );
        OpenMMPMEFEP_exposer.def( "__copy__", &__copy__);
        OpenMMPMEFEP_exposer.def( "__deepcopy__", &__copy__);
        OpenMMPMEFEP_exposer.def( "clone", &__copy__);
        OpenMMPMEFEP_exposer.def( "__str__", &pvt_get_name);
        OpenMMPMEFEP_exposer.def( "__repr__", &pvt_get_name);
    }

}
