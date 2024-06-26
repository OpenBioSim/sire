// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "CPUID.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"

#include "cpuid.h"

#include "tostring.h"

#include "cpuid.h"

SireBase::CPUID __copy__(const SireBase::CPUID &other){ return SireBase::CPUID(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CPUID_class(){

    { //::SireBase::CPUID
        typedef bp::class_< SireBase::CPUID, bp::bases< SireBase::Property > > CPUID_exposer_t;
        CPUID_exposer_t CPUID_exposer = CPUID_exposer_t( "CPUID", "This class obtains and displays the capabilities and ID of\nthe CPU at runtime\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope CPUID_scope( CPUID_exposer );
        CPUID_exposer.def( bp::init< SireBase::CPUID const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireBase::CPUID::brand
        
            typedef ::QString ( ::SireBase::CPUID::*brand_function_type)(  ) const;
            brand_function_type brand_function_value( &::SireBase::CPUID::brand );
            
            CPUID_exposer.def( 
                "brand"
                , brand_function_value
                , bp::release_gil_policy()
                , "Return the Brand string for this CPU" );
        
        }
        { //::SireBase::CPUID::clockSpeed
        
            typedef int ( ::SireBase::CPUID::*clockSpeed_function_type)(  ) const;
            clockSpeed_function_type clockSpeed_function_value( &::SireBase::CPUID::clockSpeed );
            
            CPUID_exposer.def( 
                "clockSpeed"
                , clockSpeed_function_value
                , bp::release_gil_policy()
                , "Return the clockspeed of this processor. A value of -1 is returned\nif this is not known" );
        
        }
        { //::SireBase::CPUID::numCores
        
            typedef int ( ::SireBase::CPUID::*numCores_function_type)(  ) const;
            numCores_function_type numCores_function_value( &::SireBase::CPUID::numCores );
            
            CPUID_exposer.def( 
                "numCores"
                , numCores_function_value
                , bp::release_gil_policy()
                , "Return the number of cores of this processor. A value of 1 is returned\nif this is not known (as we must have at least 1 core)" );
        
        }
        CPUID_exposer.def( bp::self != bp::self );
        { //::SireBase::CPUID::operator=
        
            typedef ::SireBase::CPUID & ( ::SireBase::CPUID::*assign_function_type)( ::SireBase::CPUID const & ) ;
            assign_function_type assign_function_value( &::SireBase::CPUID::operator= );
            
            CPUID_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CPUID_exposer.def( bp::self == bp::self );
        { //::SireBase::CPUID::supportableFeatures
        
            typedef ::QStringList ( ::SireBase::CPUID::*supportableFeatures_function_type)(  ) const;
            supportableFeatures_function_type supportableFeatures_function_value( &::SireBase::CPUID::supportableFeatures );
            
            CPUID_exposer.def( 
                "supportableFeatures"
                , supportableFeatures_function_value
                , bp::release_gil_policy()
                , "Return the list of all searchable supportable features" );
        
        }
        { //::SireBase::CPUID::supportedFeatures
        
            typedef ::QStringList ( ::SireBase::CPUID::*supportedFeatures_function_type)(  ) const;
            supportedFeatures_function_type supportedFeatures_function_value( &::SireBase::CPUID::supportedFeatures );
            
            CPUID_exposer.def( 
                "supportedFeatures"
                , supportedFeatures_function_value
                , bp::release_gil_policy()
                , "Return the list of all features supported on this CPU" );
        
        }
        { //::SireBase::CPUID::supports
        
            typedef bool ( ::SireBase::CPUID::*supports_function_type)( ::QString const & ) const;
            supports_function_type supports_function_value( &::SireBase::CPUID::supports );
            
            CPUID_exposer.def( 
                "supports"
                , supports_function_value
                , ( bp::arg("feature") )
                , bp::release_gil_policy()
                , "Returns whether or not the CPU supports the passed feature.\nNote that the passed feature must be one of the strings\nas returned by supportableFeatures" );
        
        }
        { //::SireBase::CPUID::supportsAVX
        
            typedef bool ( ::SireBase::CPUID::*supportsAVX_function_type)(  ) const;
            supportsAVX_function_type supportsAVX_function_value( &::SireBase::CPUID::supportsAVX );
            
            CPUID_exposer.def( 
                "supportsAVX"
                , supportsAVX_function_value
                , bp::release_gil_policy()
                , "Return whether or not this processor supports AVX vector instructions" );
        
        }
        { //::SireBase::CPUID::supportsSSE2
        
            typedef bool ( ::SireBase::CPUID::*supportsSSE2_function_type)(  ) const;
            supportsSSE2_function_type supportsSSE2_function_value( &::SireBase::CPUID::supportsSSE2 );
            
            CPUID_exposer.def( 
                "supportsSSE2"
                , supportsSSE2_function_value
                , bp::release_gil_policy()
                , "Return whether or not this processor supports SSE2 vector instructions" );
        
        }
        { //::SireBase::CPUID::toString
        
            typedef ::QString ( ::SireBase::CPUID::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::CPUID::toString );
            
            CPUID_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::CPUID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::CPUID::typeName );
            
            CPUID_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::CPUID::vendor
        
            typedef ::QString ( ::SireBase::CPUID::*vendor_function_type)(  ) const;
            vendor_function_type vendor_function_value( &::SireBase::CPUID::vendor );
            
            CPUID_exposer.def( 
                "vendor"
                , vendor_function_value
                , bp::release_gil_policy()
                , "Return the Vendor string for this CPU" );
        
        }
        { //::SireBase::CPUID::what
        
            typedef char const * ( ::SireBase::CPUID::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireBase::CPUID::what );
            
            CPUID_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CPUID_exposer.staticmethod( "typeName" );
        CPUID_exposer.def( "__copy__", &__copy__<SireBase::CPUID>);
        CPUID_exposer.def( "__deepcopy__", &__copy__<SireBase::CPUID>);
        CPUID_exposer.def( "clone", &__copy__<SireBase::CPUID>);
        CPUID_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::CPUID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CPUID_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::CPUID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CPUID_exposer.def_pickle(sire_pickle_suite< ::SireBase::CPUID >());
        CPUID_exposer.def( "__str__", &__str__< ::SireBase::CPUID > );
        CPUID_exposer.def( "__repr__", &__str__< ::SireBase::CPUID > );
    }

}
