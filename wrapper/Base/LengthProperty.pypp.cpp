// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "LengthProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/lengthproperty.h"

#include "SireBase/variantproperty.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "lengthproperty.h"

#include "lengthproperty.h"

SireBase::LengthProperty __copy__(const SireBase::LengthProperty &other){ return SireBase::LengthProperty(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_LengthProperty_class(){

    { //::SireBase::LengthProperty
        typedef bp::class_< SireBase::LengthProperty, bp::bases< SireBase::Property > > LengthProperty_exposer_t;
        LengthProperty_exposer_t LengthProperty_exposer = LengthProperty_exposer_t( "LengthProperty", "This class provides a thin Property wrapper around lengths\n\nThis class is deprecated, and only kept for compatibility with old S3\nfiles.\n\nYou should now use GeneralUnitProperty to hold units\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor - this constructs the integer 0") );
        bp::scope LengthProperty_scope( LengthProperty_exposer );
        LengthProperty_exposer.def( bp::init< SireUnits::Dimension::Length >(( bp::arg("value") ), "Construct from the passed length") );
        LengthProperty_exposer.def( bp::init< SireBase::LengthProperty const & >(( bp::arg("other") ), "Copy constructor") );
        LengthProperty_exposer.def( bp::init< SireBase::Property const & >(( bp::arg("other") ), "Construct from a VariantProperty") );
        LengthProperty_exposer.def( bp::self != bp::self );
        { //::SireBase::LengthProperty::operator=
        
            typedef ::SireBase::LengthProperty & ( ::SireBase::LengthProperty::*assign_function_type)( ::SireBase::LengthProperty const & ) ;
            assign_function_type assign_function_value( &::SireBase::LengthProperty::operator= );
            
            LengthProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        LengthProperty_exposer.def( bp::self == bp::self );
        { //::SireBase::LengthProperty::toString
        
            typedef ::QString ( ::SireBase::LengthProperty::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::LengthProperty::toString );
            
            LengthProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::LengthProperty::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::LengthProperty::typeName );
            
            LengthProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::LengthProperty::value
        
            typedef ::SireUnits::Dimension::Length ( ::SireBase::LengthProperty::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireBase::LengthProperty::value );
            
            LengthProperty_exposer.def( 
                "value"
                , value_function_value
                , bp::release_gil_policy()
                , "Return this number cast as a double" );
        
        }
        LengthProperty_exposer.staticmethod( "typeName" );
        LengthProperty_exposer.def( "__copy__", &__copy__<SireBase::LengthProperty>);
        LengthProperty_exposer.def( "__deepcopy__", &__copy__<SireBase::LengthProperty>);
        LengthProperty_exposer.def( "clone", &__copy__<SireBase::LengthProperty>);
        LengthProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::LengthProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LengthProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::LengthProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LengthProperty_exposer.def_pickle(sire_pickle_suite< ::SireBase::LengthProperty >());
        LengthProperty_exposer.def( "__str__", &__str__< ::SireBase::LengthProperty > );
        LengthProperty_exposer.def( "__repr__", &__str__< ::SireBase::LengthProperty > );
    }

}
