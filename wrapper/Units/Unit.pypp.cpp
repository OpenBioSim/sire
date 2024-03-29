// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Unit.pypp.hpp"

namespace bp = boost::python;

#include "dimensions.h"

#include "generalunit.h"

#include "temperature.h"

#include "dimensions.h"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireUnits::Dimension::Unit&){ return "SireUnits::Dimension::Unit";}

#include "Helpers/release_gil_policy.hpp"

void register_Unit_class(){

    { //::SireUnits::Dimension::Unit
        typedef bp::class_< SireUnits::Dimension::Unit, boost::noncopyable > Unit_exposer_t;
        Unit_exposer_t Unit_exposer = Unit_exposer_t( "Unit", "This is the base class of all units - at its heart, this is\njust a scale factor - how many times the base unit is the\ncurrent value.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Unit_scope( Unit_exposer );
        { //::SireUnits::Dimension::Unit::convertFromInternal
        
            typedef double ( ::SireUnits::Dimension::Unit::*convertFromInternal_function_type)( double ) const;
            convertFromInternal_function_type convertFromInternal_function_value( &::SireUnits::Dimension::Unit::convertFromInternal );
            
            Unit_exposer.def( 
                "convertFromInternal"
                , convertFromInternal_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireUnits::Dimension::Unit::convertToInternal
        
            typedef double ( ::SireUnits::Dimension::Unit::*convertToInternal_function_type)( double ) const;
            convertToInternal_function_type convertToInternal_function_value( &::SireUnits::Dimension::Unit::convertToInternal );
            
            Unit_exposer.def( 
                "convertToInternal"
                , convertToInternal_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireUnits::Dimension::Unit::scaleFactor
        
            typedef double ( ::SireUnits::Dimension::Unit::*scaleFactor_function_type)(  ) const;
            scaleFactor_function_type scaleFactor_function_value( &::SireUnits::Dimension::Unit::scaleFactor );
            
            Unit_exposer.def( 
                "scaleFactor"
                , scaleFactor_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireUnits::Dimension::Unit::value
        
            typedef double ( ::SireUnits::Dimension::Unit::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireUnits::Dimension::Unit::value );
            
            Unit_exposer.def( 
                "value"
                , value_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Unit_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireUnits::Dimension::Unit >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Unit_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireUnits::Dimension::Unit >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Unit_exposer.def_pickle(sire_pickle_suite< ::SireUnits::Dimension::Unit >());
        Unit_exposer.def( "__str__", &pvt_get_name);
        Unit_exposer.def( "__repr__", &pvt_get_name);
    }

}
