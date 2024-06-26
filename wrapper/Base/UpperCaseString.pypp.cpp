// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "UpperCaseString.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "stringmangler.h"

#include <QMutex>

#include "stringmangler.h"

SireBase::UpperCaseString __copy__(const SireBase::UpperCaseString &other){ return SireBase::UpperCaseString(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_UpperCaseString_class(){

    { //::SireBase::UpperCaseString
        typedef bp::class_< SireBase::UpperCaseString, bp::bases< SireBase::StringMangler, SireBase::Property > > UpperCaseString_exposer_t;
        UpperCaseString_exposer_t UpperCaseString_exposer = UpperCaseString_exposer_t( "UpperCaseString", "This mangler just trims spaces away from the beginning\nand end of the string and upper-cases the whole string\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope UpperCaseString_scope( UpperCaseString_exposer );
        UpperCaseString_exposer.def( bp::init< SireBase::UpperCaseString const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireBase::UpperCaseString::mangle
        
            typedef ::QString ( ::SireBase::UpperCaseString::*mangle_function_type)( ::QString const & ) const;
            mangle_function_type mangle_function_value( &::SireBase::UpperCaseString::mangle );
            
            UpperCaseString_exposer.def( 
                "mangle"
                , mangle_function_value
                , ( bp::arg("input") )
                , bp::release_gil_policy()
                , "Mangle the string - remove all initial and trailing spaces" );
        
        }
        UpperCaseString_exposer.def( bp::self != bp::self );
        { //::SireBase::UpperCaseString::operator=
        
            typedef ::SireBase::UpperCaseString & ( ::SireBase::UpperCaseString::*assign_function_type)( ::SireBase::UpperCaseString const & ) ;
            assign_function_type assign_function_value( &::SireBase::UpperCaseString::operator= );
            
            UpperCaseString_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        UpperCaseString_exposer.def( bp::self == bp::self );
        { //::SireBase::UpperCaseString::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::UpperCaseString::typeName );
            
            UpperCaseString_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        UpperCaseString_exposer.staticmethod( "typeName" );
        UpperCaseString_exposer.def( "__copy__", &__copy__<SireBase::UpperCaseString>);
        UpperCaseString_exposer.def( "__deepcopy__", &__copy__<SireBase::UpperCaseString>);
        UpperCaseString_exposer.def( "clone", &__copy__<SireBase::UpperCaseString>);
        UpperCaseString_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::UpperCaseString >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        UpperCaseString_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::UpperCaseString >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        UpperCaseString_exposer.def_pickle(sire_pickle_suite< ::SireBase::UpperCaseString >());
        UpperCaseString_exposer.def( "__str__", &__str__< ::SireBase::UpperCaseString > );
        UpperCaseString_exposer.def( "__repr__", &__str__< ::SireBase::UpperCaseString > );
    }

}
