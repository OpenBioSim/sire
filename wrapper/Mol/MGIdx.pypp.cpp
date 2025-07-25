// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "MGIdx.pypp.hpp"

namespace bp = boost::python;

#include "mgidx.h"

#include "mgidx.h"

#include "mgid.h"

#include "mgidx.h"

#include "mgname.h"

#include "mgnum.h"

#include "moleculegroups.h"

SireMol::MGIdx __copy__(const SireMol::MGIdx &other){ return SireMol::MGIdx(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_MGIdx_class(){

    { //::SireMol::MGIdx
        typedef bp::class_< SireMol::MGIdx, bp::bases< SireMol::MGID, SireID::ID, SireID::IndexBase > > MGIdx_exposer_t;
        MGIdx_exposer_t MGIdx_exposer = MGIdx_exposer_t( "MGIdx", "This is an ID object that is used to index molecule groups (e.g. index\nin a list or array, or in a MoleculeGroups set).\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope MGIdx_scope( MGIdx_exposer );
        MGIdx_exposer.def( bp::init< qint32 >(( bp::arg("idx") ), "") );
        MGIdx_exposer.def( bp::init< SireMol::MGIdx const & >(( bp::arg("other") ), "") );
        { //::SireMol::MGIdx::hash
        
            typedef ::uint ( ::SireMol::MGIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::MGIdx::hash );
            
            MGIdx_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MGIdx::isNull
        
            typedef bool ( ::SireMol::MGIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::MGIdx::isNull );
            
            MGIdx_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MGIdx::map
        
            typedef ::QList< SireMol::MGNum > ( ::SireMol::MGIdx::*map_function_type)( ::SireMol::MolGroupsBase const & ) const;
            map_function_type map_function_value( &::SireMol::MGIdx::map );
            
            MGIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MGIdx::null
        
            typedef ::SireMol::MGIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::MGIdx::null );
            
            MGIdx_exposer.def( 
                "null"
                , null_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MGIdx_exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireMol::MGIdx::operator=
        
            typedef ::SireMol::MGIdx & ( ::SireMol::MGIdx::*assign_function_type)( ::SireMol::MGIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::MGIdx::operator= );
            
            MGIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::MGIdx::toString
        
            typedef ::QString ( ::SireMol::MGIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::MGIdx::toString );
            
            MGIdx_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MGIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::MGIdx::typeName );
            
            MGIdx_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MGIdx::what
        
            typedef char const * ( ::SireMol::MGIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::MGIdx::what );
            
            MGIdx_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MGIdx_exposer.staticmethod( "null" );
        MGIdx_exposer.staticmethod( "typeName" );
        MGIdx_exposer.def( "__copy__", &__copy__<SireMol::MGIdx>);
        MGIdx_exposer.def( "__deepcopy__", &__copy__<SireMol::MGIdx>);
        MGIdx_exposer.def( "clone", &__copy__<SireMol::MGIdx>);
        MGIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::MGIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MGIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::MGIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MGIdx_exposer.def_pickle(sire_pickle_suite< ::SireMol::MGIdx >());
        MGIdx_exposer.def( "__str__", &__str__< ::SireMol::MGIdx > );
        MGIdx_exposer.def( "__repr__", &__str__< ::SireMol::MGIdx > );
        MGIdx_exposer.def( "__hash__", &::SireMol::MGIdx::hash );
    }

}
