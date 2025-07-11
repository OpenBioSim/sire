// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "CGName.pypp.hpp"

namespace bp = boost::python;

#include "cgname.h"

#include "cgname.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "cgid.h"

#include "cgidentifier.h"

#include "chain.h"

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

#include "cgid.h"

SireMol::CGName __copy__(const SireMol::CGName &other){ return SireMol::CGName(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CGName_class(){

    { //::SireMol::CGName
        typedef bp::class_< SireMol::CGName, bp::bases< SireMol::CGID, SireID::ID, SireID::Name > > CGName_exposer_t;
        CGName_exposer_t CGName_exposer = CGName_exposer_t( "CGName", "This class holds the name of a CutGroup.\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope CGName_scope( CGName_exposer );
        CGName_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        CGName_exposer.def( bp::init< QString const &, SireID::CaseSensitivity >(( bp::arg("name"), bp::arg("case_sensitivity") ), "") );
        CGName_exposer.def( bp::init< SireMol::CGName const & >(( bp::arg("other") ), "") );
        { //::SireMol::CGName::hash
        
            typedef ::uint ( ::SireMol::CGName::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::CGName::hash );
            
            CGName_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGName::isNull
        
            typedef bool ( ::SireMol::CGName::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::CGName::isNull );
            
            CGName_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGName::map
        
            typedef ::QList< SireMol::CGIdx > ( ::SireMol::CGName::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::CGName::map );
            
            CGName_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        CGName_exposer.def( bp::self != bp::other< SireID::ID >() );
        CGName_exposer.def( bp::self != bp::self );
        { //::SireMol::CGName::operator=
        
            typedef ::SireMol::CGName & ( ::SireMol::CGName::*assign_function_type)( ::SireMol::CGName const & ) ;
            assign_function_type assign_function_value( &::SireMol::CGName::operator= );
            
            CGName_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CGName_exposer.def( bp::self == bp::other< SireID::ID >() );
        CGName_exposer.def( bp::self == bp::self );
        { //::SireMol::CGName::toString
        
            typedef ::QString ( ::SireMol::CGName::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::CGName::toString );
            
            CGName_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGName::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CGName::typeName );
            
            CGName_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGName::what
        
            typedef char const * ( ::SireMol::CGName::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::CGName::what );
            
            CGName_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CGName_exposer.staticmethod( "typeName" );
        CGName_exposer.def( "__copy__", &__copy__<SireMol::CGName>);
        CGName_exposer.def( "__deepcopy__", &__copy__<SireMol::CGName>);
        CGName_exposer.def( "clone", &__copy__<SireMol::CGName>);
        CGName_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CGName >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGName_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CGName >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGName_exposer.def_pickle(sire_pickle_suite< ::SireMol::CGName >());
        CGName_exposer.def( "__str__", &__str__< ::SireMol::CGName > );
        CGName_exposer.def( "__repr__", &__str__< ::SireMol::CGName > );
        CGName_exposer.def( "__hash__", &::SireMol::CGName::hash );
    }

}
