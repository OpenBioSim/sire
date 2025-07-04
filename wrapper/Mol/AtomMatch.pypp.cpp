// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "AtomMatch.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/lazyevaluator.h"

#include "SireMol/core.h"

#include "SireMol/errors.h"

#include "SireMol/mover_metaid.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atommatch.h"

#include "atommatch.h"

SireMol::AtomMatch __copy__(const SireMol::AtomMatch &other){ return SireMol::AtomMatch(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_AtomMatch_class(){

    { //::SireMol::AtomMatch
        typedef bp::class_< SireMol::AtomMatch, bp::bases< SireMol::PartialMolecule, SireMol::MoleculeView, SireBase::Property > > AtomMatch_exposer_t;
        AtomMatch_exposer_t AtomMatch_exposer = AtomMatch_exposer_t( "AtomMatch", "This class holds the results of performing a match on a molecule.", bp::init< >("") );
        bp::scope AtomMatch_scope( AtomMatch_exposer );
        AtomMatch_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "") );
        AtomMatch_exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const &, QList< long long > const & >(( bp::arg("molview"), bp::arg("matches") ), "") );
        AtomMatch_exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const &, QList< QList< long long > > const & >(( bp::arg("molview"), bp::arg("matches") ), "") );
        AtomMatch_exposer.def( bp::init< SireMol::AtomMatch const & >(( bp::arg("other") ), "") );
        { //::SireMol::AtomMatch::group
        
            typedef ::SireMol::Selector< SireMol::Atom > ( ::SireMol::AtomMatch::*group_function_type)( int ) const;
            group_function_type group_function_value( &::SireMol::AtomMatch::group );
            
            AtomMatch_exposer.def( 
                "group"
                , group_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomMatch::groups
        
            typedef ::QList< SireMol::Selector< SireMol::Atom > > ( ::SireMol::AtomMatch::*groups_function_type)(  ) const;
            groups_function_type groups_function_value( &::SireMol::AtomMatch::groups );
            
            AtomMatch_exposer.def( 
                "groups"
                , groups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomMatch::nGroups
        
            typedef int ( ::SireMol::AtomMatch::*nGroups_function_type)(  ) const;
            nGroups_function_type nGroups_function_value( &::SireMol::AtomMatch::nGroups );
            
            AtomMatch_exposer.def( 
                "nGroups"
                , nGroups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomMatch_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomMatch::operator=
        
            typedef ::SireMol::AtomMatch & ( ::SireMol::AtomMatch::*assign_function_type)( ::SireMol::AtomMatch const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomMatch::operator= );
            
            AtomMatch_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomMatch_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomMatch::toString
        
            typedef ::QString ( ::SireMol::AtomMatch::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomMatch::toString );
            
            AtomMatch_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomMatch::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomMatch::typeName );
            
            AtomMatch_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomMatch::what
        
            typedef char const * ( ::SireMol::AtomMatch::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::AtomMatch::what );
            
            AtomMatch_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomMatch_exposer.staticmethod( "typeName" );
        AtomMatch_exposer.def( "__copy__", &__copy__<SireMol::AtomMatch>);
        AtomMatch_exposer.def( "__deepcopy__", &__copy__<SireMol::AtomMatch>);
        AtomMatch_exposer.def( "clone", &__copy__<SireMol::AtomMatch>);
        AtomMatch_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomMatch >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomMatch_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomMatch >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomMatch_exposer.def_pickle(sire_pickle_suite< ::SireMol::AtomMatch >());
        AtomMatch_exposer.def( "__str__", &__str__< ::SireMol::AtomMatch > );
        AtomMatch_exposer.def( "__repr__", &__str__< ::SireMol::AtomMatch > );
        AtomMatch_exposer.def( "__len__", &__len_size< ::SireMol::AtomMatch > );
    }

}
