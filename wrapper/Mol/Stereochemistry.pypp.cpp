// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Stereochemistry.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "stereochemistry.h"

#include <QDebug>

#include "stereochemistry.h"

SireMol::Stereochemistry __copy__(const SireMol::Stereochemistry &other){ return SireMol::Stereochemistry(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Stereochemistry_class(){

    { //::SireMol::Stereochemistry
        typedef bp::class_< SireMol::Stereochemistry, bp::bases< SireBase::Property > > Stereochemistry_exposer_t;
        Stereochemistry_exposer_t Stereochemistry_exposer = Stereochemistry_exposer_t( "Stereochemistry", "This class represents a bonds stereochemistry\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor (default is an undefined stereoscopy)") );
        bp::scope Stereochemistry_scope( Stereochemistry_exposer );
        Stereochemistry_exposer.def( bp::init< QString const & >(( bp::arg("s") ), "Construct from the passed string") );
        Stereochemistry_exposer.def( bp::init< SireMol::Stereochemistry const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::Stereochemistry::down
        
            typedef ::SireMol::Stereochemistry ( *down_function_type )(  );
            down_function_type down_function_value( &::SireMol::Stereochemistry::down );
            
            Stereochemistry_exposer.def( 
                "down"
                , down_function_value
                , bp::release_gil_policy()
                , "Return a down Stereochemistry" );
        
        }
        { //::SireMol::Stereochemistry::fromRDKit
        
            typedef ::SireMol::Stereochemistry ( *fromRDKit_function_type )( ::QString const & );
            fromRDKit_function_type fromRDKit_function_value( &::SireMol::Stereochemistry::fromRDKit );
            
            Stereochemistry_exposer.def( 
                "fromRDKit"
                , fromRDKit_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Construct from a string representation of a RDKit stereochemistry" );
        
        }
        { //::SireMol::Stereochemistry::fromSDF
        
            typedef ::SireMol::Stereochemistry ( *fromSDF_function_type )( int );
            fromSDF_function_type fromSDF_function_value( &::SireMol::Stereochemistry::fromSDF );
            
            Stereochemistry_exposer.def( 
                "fromSDF"
                , fromSDF_function_value
                , ( bp::arg("val") )
                , bp::release_gil_policy()
                , "Construct from the the passed SDF number" );
        
        }
        { //::SireMol::Stereochemistry::isDefined
        
            typedef bool ( ::SireMol::Stereochemistry::*isDefined_function_type)(  ) const;
            isDefined_function_type isDefined_function_value( &::SireMol::Stereochemistry::isDefined );
            
            Stereochemistry_exposer.def( 
                "isDefined"
                , isDefined_function_value
                , bp::release_gil_policy()
                , "Return whether or not the stereoscopy is defined" );
        
        }
        { //::SireMol::Stereochemistry::isDown
        
            typedef bool ( ::SireMol::Stereochemistry::*isDown_function_type)(  ) const;
            isDown_function_type isDown_function_value( &::SireMol::Stereochemistry::isDown );
            
            Stereochemistry_exposer.def( 
                "isDown"
                , isDown_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a down bond" );
        
        }
        { //::SireMol::Stereochemistry::isNotStereo
        
            typedef bool ( ::SireMol::Stereochemistry::*isNotStereo_function_type)(  ) const;
            isNotStereo_function_type isNotStereo_function_value( &::SireMol::Stereochemistry::isNotStereo );
            
            Stereochemistry_exposer.def( 
                "isNotStereo"
                , isNotStereo_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a not stereo bond" );
        
        }
        { //::SireMol::Stereochemistry::isUp
        
            typedef bool ( ::SireMol::Stereochemistry::*isUp_function_type)(  ) const;
            isUp_function_type isUp_function_value( &::SireMol::Stereochemistry::isUp );
            
            Stereochemistry_exposer.def( 
                "isUp"
                , isUp_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is an up bond" );
        
        }
        { //::SireMol::Stereochemistry::notStereo
        
            typedef ::SireMol::Stereochemistry ( *notStereo_function_type )(  );
            notStereo_function_type notStereo_function_value( &::SireMol::Stereochemistry::notStereo );
            
            Stereochemistry_exposer.def( 
                "notStereo"
                , notStereo_function_value
                , bp::release_gil_policy()
                , "Return a not stereo Stereochemistry" );
        
        }
        Stereochemistry_exposer.def( bp::self != bp::self );
        { //::SireMol::Stereochemistry::operator=
        
            typedef ::SireMol::Stereochemistry & ( ::SireMol::Stereochemistry::*assign_function_type)( ::SireMol::Stereochemistry const & ) ;
            assign_function_type assign_function_value( &::SireMol::Stereochemistry::operator= );
            
            Stereochemistry_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Stereochemistry_exposer.def( bp::self == bp::self );
        { //::SireMol::Stereochemistry::toRDKit
        
            typedef ::QString ( ::SireMol::Stereochemistry::*toRDKit_function_type)(  ) const;
            toRDKit_function_type toRDKit_function_value( &::SireMol::Stereochemistry::toRDKit );
            
            Stereochemistry_exposer.def( 
                "toRDKit"
                , toRDKit_function_value
                , bp::release_gil_policy()
                , "Return a string representation of the RDKit stereo value" );
        
        }
        { //::SireMol::Stereochemistry::toSDF
        
            typedef int ( ::SireMol::Stereochemistry::*toSDF_function_type)(  ) const;
            toSDF_function_type toSDF_function_value( &::SireMol::Stereochemistry::toSDF );
            
            Stereochemistry_exposer.def( 
                "toSDF"
                , toSDF_function_value
                , bp::release_gil_policy()
                , "Return the SDF-format value for this bond. This returns\n0 if the stereoscopy is undefined\n" );
        
        }
        { //::SireMol::Stereochemistry::toString
        
            typedef ::QString ( ::SireMol::Stereochemistry::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Stereochemistry::toString );
            
            Stereochemistry_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Stereochemistry::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Stereochemistry::typeName );
            
            Stereochemistry_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Stereochemistry::undefined
        
            typedef ::SireMol::Stereochemistry ( *undefined_function_type )(  );
            undefined_function_type undefined_function_value( &::SireMol::Stereochemistry::undefined );
            
            Stereochemistry_exposer.def( 
                "undefined"
                , undefined_function_value
                , bp::release_gil_policy()
                , "Return an undefined Stereochemistry" );
        
        }
        { //::SireMol::Stereochemistry::up
        
            typedef ::SireMol::Stereochemistry ( *up_function_type )(  );
            up_function_type up_function_value( &::SireMol::Stereochemistry::up );
            
            Stereochemistry_exposer.def( 
                "up"
                , up_function_value
                , bp::release_gil_policy()
                , "Return an up Stereochemistry" );
        
        }
        { //::SireMol::Stereochemistry::value
        
            typedef int ( ::SireMol::Stereochemistry::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireMol::Stereochemistry::value );
            
            Stereochemistry_exposer.def( 
                "value"
                , value_function_value
                , bp::release_gil_policy()
                , "Return the stereo type (uses SDF values, e.g. 0 is not stereo,\n1 is up, 6 is down. We have added -1 to mean undefined)\n" );
        
        }
        Stereochemistry_exposer.staticmethod( "down" );
        Stereochemistry_exposer.staticmethod( "fromRDKit" );
        Stereochemistry_exposer.staticmethod( "fromSDF" );
        Stereochemistry_exposer.staticmethod( "notStereo" );
        Stereochemistry_exposer.staticmethod( "typeName" );
        Stereochemistry_exposer.staticmethod( "undefined" );
        Stereochemistry_exposer.staticmethod( "up" );
        Stereochemistry_exposer.def( "__copy__", &__copy__);
        Stereochemistry_exposer.def( "__deepcopy__", &__copy__);
        Stereochemistry_exposer.def( "clone", &__copy__);
        Stereochemistry_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Stereochemistry >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Stereochemistry_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Stereochemistry >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Stereochemistry_exposer.def_pickle(sire_pickle_suite< ::SireMol::Stereochemistry >());
        Stereochemistry_exposer.def( "__str__", &__str__< ::SireMol::Stereochemistry > );
        Stereochemistry_exposer.def( "__repr__", &__str__< ::SireMol::Stereochemistry > );
    }

}