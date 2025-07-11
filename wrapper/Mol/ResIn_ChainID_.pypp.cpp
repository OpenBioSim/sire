// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ResIn_ChainID_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "atom.h"

#include "chain.h"

#include "chainid.h"

#include "chainidentifier.h"

#include "chainresid.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "groupgroupids.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "chainid.h"

SireMol::ResIn<SireMol::ChainID> __copy__(const SireMol::ResIn<SireMol::ChainID> &other){ return SireMol::ResIn<SireMol::ChainID>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_ResIn_ChainID__class(){

    { //::SireMol::ResIn< SireMol::ChainID >
        typedef bp::class_< SireMol::ResIn< SireMol::ChainID >, bp::bases< SireMol::ResID, SireID::ID > > ResIn_ChainID__exposer_t;
        ResIn_ChainID__exposer_t ResIn_ChainID__exposer = ResIn_ChainID__exposer_t( "ResIn_ChainID_", "", bp::init< >("") );
        bp::scope ResIn_ChainID__scope( ResIn_ChainID__exposer );
        ResIn_ChainID__exposer.def( bp::init< SireMol::ChainID const & >(( bp::arg("id") ), "") );
        ResIn_ChainID__exposer.def( bp::init< SireMol::ChainID const &, qint32 >(( bp::arg("id"), bp::arg("i") ), "") );
        ResIn_ChainID__exposer.def( bp::init< SireMol::ChainID const &, qint32, qint32 >(( bp::arg("id"), bp::arg("i"), bp::arg("j") ), "") );
        ResIn_ChainID__exposer.def( bp::init< SireMol::ResIn< SireMol::ChainID > const & >(( bp::arg("other") ), "") );
        { //::SireMol::ResIn< SireMol::ChainID >::hash
        
            typedef SireMol::ResIn< SireMol::ChainID > exported_class_t;
            typedef ::uint ( ::SireMol::ResIn< SireMol::ChainID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::ResIn< SireMol::ChainID >::hash );
            
            ResIn_ChainID__exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResIn< SireMol::ChainID >::isNull
        
            typedef SireMol::ResIn< SireMol::ChainID > exported_class_t;
            typedef bool ( ::SireMol::ResIn< SireMol::ChainID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::ResIn< SireMol::ChainID >::isNull );
            
            ResIn_ChainID__exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResIn< SireMol::ChainID >::map
        
            typedef SireMol::ResIn< SireMol::ChainID > exported_class_t;
            typedef ::QList< SireMol::ResIdx > ( ::SireMol::ResIn< SireMol::ChainID >::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::ResIn< SireMol::ChainID >::map );
            
            ResIn_ChainID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        ResIn_ChainID__exposer.def( bp::self != bp::self );
        ResIn_ChainID__exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireMol::ResIn< SireMol::ChainID >::operator=
        
            typedef SireMol::ResIn< SireMol::ChainID > exported_class_t;
            typedef ::SireMol::ResIn< SireMol::ChainID > & ( ::SireMol::ResIn< SireMol::ChainID >::*assign_function_type)( ::SireMol::ResIn< SireMol::ChainID > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResIn< SireMol::ChainID >::operator= );
            
            ResIn_ChainID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ResIn_ChainID__exposer.def( bp::self == bp::self );
        ResIn_ChainID__exposer.def( bp::self == bp::other< SireID::ID >() );
        { //::SireMol::ResIn< SireMol::ChainID >::toString
        
            typedef SireMol::ResIn< SireMol::ChainID > exported_class_t;
            typedef ::QString ( ::SireMol::ResIn< SireMol::ChainID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResIn< SireMol::ChainID >::toString );
            
            ResIn_ChainID__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResIn< SireMol::ChainID >::typeName
        
            typedef SireMol::ResIn< SireMol::ChainID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResIn< SireMol::ChainID >::typeName );
            
            ResIn_ChainID__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResIn< SireMol::ChainID >::what
        
            typedef SireMol::ResIn< SireMol::ChainID > exported_class_t;
            typedef char const * ( ::SireMol::ResIn< SireMol::ChainID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ResIn< SireMol::ChainID >::what );
            
            ResIn_ChainID__exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ResIn_ChainID__exposer.staticmethod( "typeName" );
        ResIn_ChainID__exposer.def( "__copy__", &__copy__<SireMol::ResIn<SireMol::ChainID>>);
        ResIn_ChainID__exposer.def( "__deepcopy__", &__copy__<SireMol::ResIn<SireMol::ChainID>>);
        ResIn_ChainID__exposer.def( "clone", &__copy__<SireMol::ResIn<SireMol::ChainID>>);
        ResIn_ChainID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResIn<SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIn_ChainID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResIn<SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIn_ChainID__exposer.def_pickle(sire_pickle_suite< ::SireMol::ResIn<SireMol::ChainID> >());
        ResIn_ChainID__exposer.def( "__str__", &__str__< ::SireMol::ResIn<SireMol::ChainID> > );
        ResIn_ChainID__exposer.def( "__repr__", &__str__< ::SireMol::ResIn<SireMol::ChainID> > );
        ResIn_ChainID__exposer.def( "__hash__", &::SireMol::ResIn<SireMol::ChainID>::hash );
    }

}
