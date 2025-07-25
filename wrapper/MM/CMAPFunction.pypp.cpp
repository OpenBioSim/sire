// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CMAPFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/console.h"

#include "SireError/errors.h"

#include "SireMol/atommapping.h"

#include "SireMol/atommatcher.h"

#include "SireMol/atomselection.h"

#include "SireMol/errors.h"

#include "SireMol/moleculedata.h"

#include "SireMol/moleculeinfodata.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "cmapfunctions.h"

#include "sireglobal.h"

#include "cmapfunctions.h"

SireMM::CMAPFunction __copy__(const SireMM::CMAPFunction &other){ return SireMM::CMAPFunction(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CMAPFunction_class(){

    { //::SireMM::CMAPFunction
        typedef bp::class_< SireMM::CMAPFunction > CMAPFunction_exposer_t;
        CMAPFunction_exposer_t CMAPFunction_exposer = CMAPFunction_exposer_t( "CMAPFunction", "This class holds a single CMAP function for a single set\nof five atoms\n", bp::init< >("") );
        bp::scope CMAPFunction_scope( CMAPFunction_exposer );
        CMAPFunction_exposer.def( bp::init< SireMol::CGAtomIdx const &, SireMol::CGAtomIdx const &, SireMol::CGAtomIdx const &, SireMol::CGAtomIdx const &, SireMol::CGAtomIdx const &, SireMM::CMAPParameter const & >(( bp::arg("atm0"), bp::arg("atm1"), bp::arg("atm2"), bp::arg("atm3"), bp::arg("atm4"), bp::arg("param") ), "") );
        CMAPFunction_exposer.def( bp::init< SireMM::CMAPFunction const & >(( bp::arg("other") ), "") );
        { //::SireMM::CMAPFunction::atom0
        
            typedef ::SireMol::CGAtomIdx const & ( ::SireMM::CMAPFunction::*atom0_function_type)(  ) const;
            atom0_function_type atom0_function_value( &::SireMM::CMAPFunction::atom0 );
            
            CMAPFunction_exposer.def( 
                "atom0"
                , atom0_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CMAPFunction::atom1
        
            typedef ::SireMol::CGAtomIdx const & ( ::SireMM::CMAPFunction::*atom1_function_type)(  ) const;
            atom1_function_type atom1_function_value( &::SireMM::CMAPFunction::atom1 );
            
            CMAPFunction_exposer.def( 
                "atom1"
                , atom1_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CMAPFunction::atom2
        
            typedef ::SireMol::CGAtomIdx const & ( ::SireMM::CMAPFunction::*atom2_function_type)(  ) const;
            atom2_function_type atom2_function_value( &::SireMM::CMAPFunction::atom2 );
            
            CMAPFunction_exposer.def( 
                "atom2"
                , atom2_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CMAPFunction::atom3
        
            typedef ::SireMol::CGAtomIdx const & ( ::SireMM::CMAPFunction::*atom3_function_type)(  ) const;
            atom3_function_type atom3_function_value( &::SireMM::CMAPFunction::atom3 );
            
            CMAPFunction_exposer.def( 
                "atom3"
                , atom3_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CMAPFunction::atom4
        
            typedef ::SireMol::CGAtomIdx const & ( ::SireMM::CMAPFunction::*atom4_function_type)(  ) const;
            atom4_function_type atom4_function_value( &::SireMM::CMAPFunction::atom4 );
            
            CMAPFunction_exposer.def( 
                "atom4"
                , atom4_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        CMAPFunction_exposer.def( bp::self != bp::self );
        { //::SireMM::CMAPFunction::operator=
        
            typedef ::SireMM::CMAPFunction & ( ::SireMM::CMAPFunction::*assign_function_type)( ::SireMM::CMAPFunction const & ) ;
            assign_function_type assign_function_value( &::SireMM::CMAPFunction::operator= );
            
            CMAPFunction_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CMAPFunction_exposer.def( bp::self == bp::self );
        { //::SireMM::CMAPFunction::parameter
        
            typedef ::SireMM::CMAPParameter const & ( ::SireMM::CMAPFunction::*parameter_function_type)(  ) const;
            parameter_function_type parameter_function_value( &::SireMM::CMAPFunction::parameter );
            
            CMAPFunction_exposer.def( 
                "parameter"
                , parameter_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::CMAPFunction::toString
        
            typedef ::QString ( ::SireMM::CMAPFunction::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::CMAPFunction::toString );
            
            CMAPFunction_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CMAPFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CMAPFunction::typeName );
            
            CMAPFunction_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CMAPFunction::what
        
            typedef char const * ( ::SireMM::CMAPFunction::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::CMAPFunction::what );
            
            CMAPFunction_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CMAPFunction_exposer.staticmethod( "typeName" );
        CMAPFunction_exposer.def( "__copy__", &__copy__<SireMM::CMAPFunction>);
        CMAPFunction_exposer.def( "__deepcopy__", &__copy__<SireMM::CMAPFunction>);
        CMAPFunction_exposer.def( "clone", &__copy__<SireMM::CMAPFunction>);
        CMAPFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CMAPFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CMAPFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CMAPFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CMAPFunction_exposer.def_pickle(sire_pickle_suite< ::SireMM::CMAPFunction >());
        CMAPFunction_exposer.def( "__str__", &__str__< ::SireMM::CMAPFunction > );
        CMAPFunction_exposer.def( "__repr__", &__str__< ::SireMM::CMAPFunction > );
    }

}
