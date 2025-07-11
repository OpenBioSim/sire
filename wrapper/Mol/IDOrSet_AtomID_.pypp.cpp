// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "IDOrSet_AtomID_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomid.h"

#include "atomidentifier.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "molatomid.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "withatoms.h"

#include <QDebug>

#include "atomid.h"

SireID::IDOrSet<SireMol::AtomID> __copy__(const SireID::IDOrSet<SireMol::AtomID> &other){ return SireID::IDOrSet<SireMol::AtomID>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_IDOrSet_AtomID__class(){

    { //::SireID::IDOrSet< SireMol::AtomID >
        typedef bp::class_< SireID::IDOrSet< SireMol::AtomID >, bp::bases< SireMol::AtomID, SireID::ID > > IDOrSet_AtomID__exposer_t;
        IDOrSet_AtomID__exposer_t IDOrSet_AtomID__exposer = IDOrSet_AtomID__exposer_t( "IDOrSet_AtomID_", "", bp::init< >("") );
        bp::scope IDOrSet_AtomID__scope( IDOrSet_AtomID__exposer );
        IDOrSet_AtomID__exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("id") ), "") );
        IDOrSet_AtomID__exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("id0"), bp::arg("id1") ), "") );
        IDOrSet_AtomID__exposer.def( bp::init< QList< SireMol::AtomIdentifier > const & >(( bp::arg("ids") ), "") );
        IDOrSet_AtomID__exposer.def( bp::init< SireID::IDOrSet< SireMol::AtomID > const & >(( bp::arg("ids") ), "") );
        IDOrSet_AtomID__exposer.def( bp::init< SireID::IDOrSet< SireMol::AtomID > const & >(( bp::arg("other") ), "") );
        { //::SireID::IDOrSet< SireMol::AtomID >::IDs
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef ::QSet< SireMol::AtomIdentifier > const & ( ::SireID::IDOrSet< SireMol::AtomID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDOrSet< SireMol::AtomID >::IDs );
            
            IDOrSet_AtomID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::AtomID >::hash
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef ::uint ( ::SireID::IDOrSet< SireMol::AtomID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDOrSet< SireMol::AtomID >::hash );
            
            IDOrSet_AtomID__exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::AtomID >::isNull
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef bool ( ::SireID::IDOrSet< SireMol::AtomID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDOrSet< SireMol::AtomID >::isNull );
            
            IDOrSet_AtomID__exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::AtomID >::map
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireID::IDOrSet< SireMol::AtomID >::*map_function_type)( ::SireID::IDOrSet< SireMol::AtomID >::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::IDOrSet< SireMol::AtomID >::map );
            
            IDOrSet_AtomID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") )
                , bp::release_gil_policy()
                , "" );
        
        }
        IDOrSet_AtomID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDOrSet_AtomID__exposer.def( bp::self != bp::self );
        IDOrSet_AtomID__exposer.def( bp::self != bp::other< SireMol::AtomID >() );
        { //::SireID::IDOrSet< SireMol::AtomID >::operator=
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef ::SireID::IDOrSet< SireMol::AtomID > & ( ::SireID::IDOrSet< SireMol::AtomID >::*assign_function_type)( ::SireID::IDOrSet< SireMol::AtomID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDOrSet< SireMol::AtomID >::operator= );
            
            IDOrSet_AtomID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::AtomID >::operator=
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef ::SireID::IDOrSet< SireMol::AtomID > & ( ::SireID::IDOrSet< SireMol::AtomID >::*assign_function_type)( ::SireMol::AtomID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDOrSet< SireMol::AtomID >::operator= );
            
            IDOrSet_AtomID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IDOrSet_AtomID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDOrSet_AtomID__exposer.def( bp::self == bp::self );
        IDOrSet_AtomID__exposer.def( bp::self == bp::other< SireMol::AtomID >() );
        { //::SireID::IDOrSet< SireMol::AtomID >::toString
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef ::QString ( ::SireID::IDOrSet< SireMol::AtomID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDOrSet< SireMol::AtomID >::toString );
            
            IDOrSet_AtomID__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::AtomID >::typeName
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDOrSet< SireMol::AtomID >::typeName );
            
            IDOrSet_AtomID__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::AtomID >::what
        
            typedef SireID::IDOrSet< SireMol::AtomID > exported_class_t;
            typedef char const * ( ::SireID::IDOrSet< SireMol::AtomID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDOrSet< SireMol::AtomID >::what );
            
            IDOrSet_AtomID__exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IDOrSet_AtomID__exposer.staticmethod( "typeName" );
        IDOrSet_AtomID__exposer.def( "__copy__", &__copy__<SireID::IDOrSet<SireMol::AtomID>>);
        IDOrSet_AtomID__exposer.def( "__deepcopy__", &__copy__<SireID::IDOrSet<SireMol::AtomID>>);
        IDOrSet_AtomID__exposer.def( "clone", &__copy__<SireID::IDOrSet<SireMol::AtomID>>);
        IDOrSet_AtomID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDOrSet<SireMol::AtomID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDOrSet_AtomID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDOrSet<SireMol::AtomID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDOrSet_AtomID__exposer.def_pickle(sire_pickle_suite< ::SireID::IDOrSet<SireMol::AtomID> >());
        IDOrSet_AtomID__exposer.def( "__str__", &__str__< ::SireID::IDOrSet<SireMol::AtomID> > );
        IDOrSet_AtomID__exposer.def( "__repr__", &__str__< ::SireID::IDOrSet<SireMol::AtomID> > );
        IDOrSet_AtomID__exposer.def( "__hash__", &::SireID::IDOrSet<SireMol::AtomID>::hash );
    }

}
