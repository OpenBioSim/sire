// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "NoMangling.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "stringmangler.h"

#include <QMutex>

#include "stringmangler.h"

SireBase::NoMangling __copy__(const SireBase::NoMangling &other){ return SireBase::NoMangling(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_NoMangling_class(){

    { //::SireBase::NoMangling
        typedef bp::class_< SireBase::NoMangling, bp::bases< SireBase::StringMangler, SireBase::Property > > NoMangling_exposer_t;
        NoMangling_exposer_t NoMangling_exposer = NoMangling_exposer_t( "NoMangling", "This mangler does absolutely nothing to the string\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope NoMangling_scope( NoMangling_exposer );
        NoMangling_exposer.def( bp::init< SireBase::NoMangling const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireBase::NoMangling::mangle
        
            typedef ::QString ( ::SireBase::NoMangling::*mangle_function_type)( ::QString const & ) const;
            mangle_function_type mangle_function_value( &::SireBase::NoMangling::mangle );
            
            NoMangling_exposer.def( 
                "mangle"
                , mangle_function_value
                , ( bp::arg("input") )
                , bp::release_gil_policy()
                , "Mangle the string - remove all initial and trailing spaces" );
        
        }
        NoMangling_exposer.def( bp::self != bp::self );
        { //::SireBase::NoMangling::operator=
        
            typedef ::SireBase::NoMangling & ( ::SireBase::NoMangling::*assign_function_type)( ::SireBase::NoMangling const & ) ;
            assign_function_type assign_function_value( &::SireBase::NoMangling::operator= );
            
            NoMangling_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        NoMangling_exposer.def( bp::self == bp::self );
        { //::SireBase::NoMangling::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::NoMangling::typeName );
            
            NoMangling_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        NoMangling_exposer.staticmethod( "typeName" );
        NoMangling_exposer.def( "__copy__", &__copy__<SireBase::NoMangling>);
        NoMangling_exposer.def( "__deepcopy__", &__copy__<SireBase::NoMangling>);
        NoMangling_exposer.def( "clone", &__copy__<SireBase::NoMangling>);
        NoMangling_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::NoMangling >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NoMangling_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::NoMangling >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        NoMangling_exposer.def_pickle(sire_pickle_suite< ::SireBase::NoMangling >());
        NoMangling_exposer.def( "__str__", &__str__< ::SireBase::NoMangling > );
        NoMangling_exposer.def( "__repr__", &__str__< ::SireBase::NoMangling > );
    }

}
