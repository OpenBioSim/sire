// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "IDOrSet_ChainID_.pypp.hpp"

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

SireID::IDOrSet<SireMol::ChainID> __copy__(const SireID::IDOrSet<SireMol::ChainID> &other){ return SireID::IDOrSet<SireMol::ChainID>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_IDOrSet_ChainID__class(){

    { //::SireID::IDOrSet< SireMol::ChainID >
        typedef bp::class_< SireID::IDOrSet< SireMol::ChainID >, bp::bases< SireMol::ChainID, SireID::ID > > IDOrSet_ChainID__exposer_t;
        IDOrSet_ChainID__exposer_t IDOrSet_ChainID__exposer = IDOrSet_ChainID__exposer_t( "IDOrSet_ChainID_", "", bp::init< >("") );
        bp::scope IDOrSet_ChainID__scope( IDOrSet_ChainID__exposer );
        IDOrSet_ChainID__exposer.def( bp::init< SireMol::ChainID const & >(( bp::arg("id") ), "") );
        IDOrSet_ChainID__exposer.def( bp::init< SireMol::ChainID const &, SireMol::ChainID const & >(( bp::arg("id0"), bp::arg("id1") ), "") );
        IDOrSet_ChainID__exposer.def( bp::init< QList< SireMol::ChainIdentifier > const & >(( bp::arg("ids") ), "") );
        IDOrSet_ChainID__exposer.def( bp::init< SireID::IDOrSet< SireMol::ChainID > const & >(( bp::arg("ids") ), "") );
        IDOrSet_ChainID__exposer.def( bp::init< SireID::IDOrSet< SireMol::ChainID > const & >(( bp::arg("other") ), "") );
        { //::SireID::IDOrSet< SireMol::ChainID >::IDs
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef ::QSet< SireMol::ChainIdentifier > const & ( ::SireID::IDOrSet< SireMol::ChainID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDOrSet< SireMol::ChainID >::IDs );
            
            IDOrSet_ChainID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::ChainID >::hash
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef ::uint ( ::SireID::IDOrSet< SireMol::ChainID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDOrSet< SireMol::ChainID >::hash );
            
            IDOrSet_ChainID__exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::ChainID >::isNull
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef bool ( ::SireID::IDOrSet< SireMol::ChainID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDOrSet< SireMol::ChainID >::isNull );
            
            IDOrSet_ChainID__exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::ChainID >::map
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef ::QList< SireMol::ChainIdx > ( ::SireID::IDOrSet< SireMol::ChainID >::*map_function_type)( ::SireID::IDOrSet< SireMol::ChainID >::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::IDOrSet< SireMol::ChainID >::map );
            
            IDOrSet_ChainID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") )
                , bp::release_gil_policy()
                , "" );
        
        }
        IDOrSet_ChainID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDOrSet_ChainID__exposer.def( bp::self != bp::self );
        IDOrSet_ChainID__exposer.def( bp::self != bp::other< SireMol::ChainID >() );
        { //::SireID::IDOrSet< SireMol::ChainID >::operator=
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef ::SireID::IDOrSet< SireMol::ChainID > & ( ::SireID::IDOrSet< SireMol::ChainID >::*assign_function_type)( ::SireID::IDOrSet< SireMol::ChainID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDOrSet< SireMol::ChainID >::operator= );
            
            IDOrSet_ChainID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::ChainID >::operator=
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef ::SireID::IDOrSet< SireMol::ChainID > & ( ::SireID::IDOrSet< SireMol::ChainID >::*assign_function_type)( ::SireMol::ChainID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDOrSet< SireMol::ChainID >::operator= );
            
            IDOrSet_ChainID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IDOrSet_ChainID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDOrSet_ChainID__exposer.def( bp::self == bp::self );
        IDOrSet_ChainID__exposer.def( bp::self == bp::other< SireMol::ChainID >() );
        { //::SireID::IDOrSet< SireMol::ChainID >::toString
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef ::QString ( ::SireID::IDOrSet< SireMol::ChainID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDOrSet< SireMol::ChainID >::toString );
            
            IDOrSet_ChainID__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::ChainID >::typeName
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDOrSet< SireMol::ChainID >::typeName );
            
            IDOrSet_ChainID__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::ChainID >::what
        
            typedef SireID::IDOrSet< SireMol::ChainID > exported_class_t;
            typedef char const * ( ::SireID::IDOrSet< SireMol::ChainID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDOrSet< SireMol::ChainID >::what );
            
            IDOrSet_ChainID__exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IDOrSet_ChainID__exposer.staticmethod( "typeName" );
        IDOrSet_ChainID__exposer.def( "__copy__", &__copy__<SireID::IDOrSet<SireMol::ChainID>>);
        IDOrSet_ChainID__exposer.def( "__deepcopy__", &__copy__<SireID::IDOrSet<SireMol::ChainID>>);
        IDOrSet_ChainID__exposer.def( "clone", &__copy__<SireID::IDOrSet<SireMol::ChainID>>);
        IDOrSet_ChainID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDOrSet<SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDOrSet_ChainID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDOrSet<SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDOrSet_ChainID__exposer.def_pickle(sire_pickle_suite< ::SireID::IDOrSet<SireMol::ChainID> >());
        IDOrSet_ChainID__exposer.def( "__str__", &__str__< ::SireID::IDOrSet<SireMol::ChainID> > );
        IDOrSet_ChainID__exposer.def( "__repr__", &__str__< ::SireID::IDOrSet<SireMol::ChainID> > );
        IDOrSet_ChainID__exposer.def( "__hash__", &::SireID::IDOrSet<SireMol::ChainID>::hash );
    }

}
