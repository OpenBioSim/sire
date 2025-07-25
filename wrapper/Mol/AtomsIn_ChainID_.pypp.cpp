// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "AtomsIn_ChainID_.pypp.hpp"

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

SireMol::AtomsIn<SireMol::ChainID> __copy__(const SireMol::AtomsIn<SireMol::ChainID> &other){ return SireMol::AtomsIn<SireMol::ChainID>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AtomsIn_ChainID__class(){

    { //::SireMol::AtomsIn< SireMol::ChainID >
        typedef bp::class_< SireMol::AtomsIn< SireMol::ChainID >, bp::bases< SireMol::AtomID, SireID::ID > > AtomsIn_ChainID__exposer_t;
        AtomsIn_ChainID__exposer_t AtomsIn_ChainID__exposer = AtomsIn_ChainID__exposer_t( "AtomsIn_ChainID_", "", bp::init< >("") );
        bp::scope AtomsIn_ChainID__scope( AtomsIn_ChainID__exposer );
        AtomsIn_ChainID__exposer.def( bp::init< SireMol::ChainID const & >(( bp::arg("id") ), "") );
        AtomsIn_ChainID__exposer.def( bp::init< SireMol::ChainID const &, qint32 >(( bp::arg("id"), bp::arg("i") ), "") );
        AtomsIn_ChainID__exposer.def( bp::init< SireMol::ChainID const &, qint32, qint32 >(( bp::arg("id"), bp::arg("i"), bp::arg("j") ), "") );
        AtomsIn_ChainID__exposer.def( bp::init< SireMol::AtomsIn< SireMol::ChainID > const & >(( bp::arg("other") ), "") );
        { //::SireMol::AtomsIn< SireMol::ChainID >::hash
        
            typedef SireMol::AtomsIn< SireMol::ChainID > exported_class_t;
            typedef ::uint ( ::SireMol::AtomsIn< SireMol::ChainID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::AtomsIn< SireMol::ChainID >::hash );
            
            AtomsIn_ChainID__exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::ChainID >::isNull
        
            typedef SireMol::AtomsIn< SireMol::ChainID > exported_class_t;
            typedef bool ( ::SireMol::AtomsIn< SireMol::ChainID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::AtomsIn< SireMol::ChainID >::isNull );
            
            AtomsIn_ChainID__exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::ChainID >::map
        
            typedef SireMol::AtomsIn< SireMol::ChainID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::AtomsIn< SireMol::ChainID >::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::AtomsIn< SireMol::ChainID >::map );
            
            AtomsIn_ChainID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomsIn_ChainID__exposer.def( bp::self != bp::self );
        AtomsIn_ChainID__exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::ChainID >::operator=
        
            typedef SireMol::AtomsIn< SireMol::ChainID > exported_class_t;
            typedef ::SireMol::AtomsIn< SireMol::ChainID > & ( ::SireMol::AtomsIn< SireMol::ChainID >::*assign_function_type)( ::SireMol::AtomsIn< SireMol::ChainID > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomsIn< SireMol::ChainID >::operator= );
            
            AtomsIn_ChainID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomsIn_ChainID__exposer.def( bp::self == bp::self );
        AtomsIn_ChainID__exposer.def( bp::self == bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::ChainID >::toString
        
            typedef SireMol::AtomsIn< SireMol::ChainID > exported_class_t;
            typedef ::QString ( ::SireMol::AtomsIn< SireMol::ChainID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomsIn< SireMol::ChainID >::toString );
            
            AtomsIn_ChainID__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::ChainID >::typeName
        
            typedef SireMol::AtomsIn< SireMol::ChainID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomsIn< SireMol::ChainID >::typeName );
            
            AtomsIn_ChainID__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::ChainID >::what
        
            typedef SireMol::AtomsIn< SireMol::ChainID > exported_class_t;
            typedef char const * ( ::SireMol::AtomsIn< SireMol::ChainID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::AtomsIn< SireMol::ChainID >::what );
            
            AtomsIn_ChainID__exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomsIn_ChainID__exposer.staticmethod( "typeName" );
        AtomsIn_ChainID__exposer.def( "__copy__", &__copy__<SireMol::AtomsIn<SireMol::ChainID>>);
        AtomsIn_ChainID__exposer.def( "__deepcopy__", &__copy__<SireMol::AtomsIn<SireMol::ChainID>>);
        AtomsIn_ChainID__exposer.def( "clone", &__copy__<SireMol::AtomsIn<SireMol::ChainID>>);
        AtomsIn_ChainID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomsIn<SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_ChainID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomsIn<SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_ChainID__exposer.def_pickle(sire_pickle_suite< ::SireMol::AtomsIn<SireMol::ChainID> >());
        AtomsIn_ChainID__exposer.def( "__str__", &__str__< ::SireMol::AtomsIn<SireMol::ChainID> > );
        AtomsIn_ChainID__exposer.def( "__repr__", &__str__< ::SireMol::AtomsIn<SireMol::ChainID> > );
        AtomsIn_ChainID__exposer.def( "__hash__", &::SireMol::AtomsIn<SireMol::ChainID>::hash );
    }

}
