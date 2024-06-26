// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "TimeProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/timeproperty.h"

#include "SireBase/variantproperty.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "timeproperty.h"

#include "timeproperty.h"

SireBase::TimeProperty __copy__(const SireBase::TimeProperty &other){ return SireBase::TimeProperty(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_TimeProperty_class(){

    { //::SireBase::TimeProperty
        typedef bp::class_< SireBase::TimeProperty, bp::bases< SireBase::Property > > TimeProperty_exposer_t;
        TimeProperty_exposer_t TimeProperty_exposer = TimeProperty_exposer_t( "TimeProperty", "This class provides a thin Property wrapper around times\n\nThis class is deprecated and only kept for compatibility with\nold S3 files.\n\nNow you should use GeneralUnitProperty for all units\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor - this constructs the integer 0") );
        bp::scope TimeProperty_scope( TimeProperty_exposer );
        TimeProperty_exposer.def( bp::init< SireUnits::Dimension::Time >(( bp::arg("value") ), "Construct from the passed length") );
        TimeProperty_exposer.def( bp::init< SireBase::TimeProperty const & >(( bp::arg("other") ), "Copy constructor") );
        TimeProperty_exposer.def( bp::init< SireBase::Property const & >(( bp::arg("other") ), "Construct from a Property") );
        TimeProperty_exposer.def( bp::self != bp::self );
        { //::SireBase::TimeProperty::operator=
        
            typedef ::SireBase::TimeProperty & ( ::SireBase::TimeProperty::*assign_function_type)( ::SireBase::TimeProperty const & ) ;
            assign_function_type assign_function_value( &::SireBase::TimeProperty::operator= );
            
            TimeProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        TimeProperty_exposer.def( bp::self == bp::self );
        { //::SireBase::TimeProperty::toString
        
            typedef ::QString ( ::SireBase::TimeProperty::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::TimeProperty::toString );
            
            TimeProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::TimeProperty::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::TimeProperty::typeName );
            
            TimeProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::TimeProperty::value
        
            typedef ::SireUnits::Dimension::Time ( ::SireBase::TimeProperty::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireBase::TimeProperty::value );
            
            TimeProperty_exposer.def( 
                "value"
                , value_function_value
                , bp::release_gil_policy()
                , "Return this number cast as a Time" );
        
        }
        TimeProperty_exposer.staticmethod( "typeName" );
        TimeProperty_exposer.def( "__copy__", &__copy__<SireBase::TimeProperty>);
        TimeProperty_exposer.def( "__deepcopy__", &__copy__<SireBase::TimeProperty>);
        TimeProperty_exposer.def( "clone", &__copy__<SireBase::TimeProperty>);
        TimeProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::TimeProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TimeProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::TimeProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TimeProperty_exposer.def_pickle(sire_pickle_suite< ::SireBase::TimeProperty >());
        TimeProperty_exposer.def( "__str__", &__str__< ::SireBase::TimeProperty > );
        TimeProperty_exposer.def( "__repr__", &__str__< ::SireBase::TimeProperty > );
    }

}
