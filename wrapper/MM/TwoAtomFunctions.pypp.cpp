// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "TwoAtomFunctions.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/symbols.h"

#include "SireError/errors.h"

#include "SireMol/atommatcher.h"

#include "SireMol/atomselection.h"

#include "SireMol/errors.h"

#include "SireMol/moleculeinfodata.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "sireglobal.h"

#include "twoatomfunctions.h"

#include "twoatomfunctions.h"

#include "SireMol/moleculedata.h"

SireMM::TwoAtomFunctions __copy__(const SireMM::TwoAtomFunctions &other){ return SireMM::TwoAtomFunctions(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_TwoAtomFunctions_class(){

    { //::SireMM::TwoAtomFunctions
        typedef bp::class_< SireMM::TwoAtomFunctions, bp::bases< SireMM::AtomFunctions, SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > TwoAtomFunctions_exposer_t;
        TwoAtomFunctions_exposer_t TwoAtomFunctions_exposer = TwoAtomFunctions_exposer_t( "TwoAtomFunctions", bp::init< >() );
        bp::scope TwoAtomFunctions_scope( TwoAtomFunctions_exposer );
        TwoAtomFunctions_exposer.def( bp::init< SireMol::MoleculeData const & >(( bp::arg("moldata") )) );
        TwoAtomFunctions_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") )) );
        TwoAtomFunctions_exposer.def( bp::init< SireMM::TwoAtomFunctions const & >(( bp::arg("other") )) );
        { //::SireMM::TwoAtomFunctions::clear
        
            typedef void ( ::SireMM::TwoAtomFunctions::*clear_function_type)( ::SireMol::AtomIdx ) ;
            clear_function_type clear_function_value( &::SireMM::TwoAtomFunctions::clear );
            
            TwoAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::clear
        
            typedef void ( ::SireMM::TwoAtomFunctions::*clear_function_type)( ::SireMol::AtomID const & ) ;
            clear_function_type clear_function_value( &::SireMM::TwoAtomFunctions::clear );
            
            TwoAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::clear
        
            typedef void ( ::SireMM::TwoAtomFunctions::*clear_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx ) ;
            clear_function_type clear_function_value( &::SireMM::TwoAtomFunctions::clear );
            
            TwoAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::clear
        
            typedef void ( ::SireMM::TwoAtomFunctions::*clear_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) ;
            clear_function_type clear_function_value( &::SireMM::TwoAtomFunctions::clear );
            
            TwoAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::clear
        
            typedef void ( ::SireMM::TwoAtomFunctions::*clear_function_type)( ::SireMol::BondID const & ) ;
            clear_function_type clear_function_value( &::SireMM::TwoAtomFunctions::clear );
            
            TwoAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value
                , ( bp::arg("bondid") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::clear
        
            typedef void ( ::SireMM::TwoAtomFunctions::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireMM::TwoAtomFunctions::clear );
            
            TwoAtomFunctions_exposer.def( 
                "clear"
                , clear_function_value );
        
        }
        { //::SireMM::TwoAtomFunctions::force
        
            typedef ::SireCAS::Expression ( ::SireMM::TwoAtomFunctions::*force_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx,::SireCAS::Symbol const & ) const;
            force_function_type force_function_value( &::SireMM::TwoAtomFunctions::force );
            
            TwoAtomFunctions_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("symbol") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::force
        
            typedef ::SireCAS::Expression ( ::SireMM::TwoAtomFunctions::*force_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireCAS::Symbol const & ) const;
            force_function_type force_function_value( &::SireMM::TwoAtomFunctions::force );
            
            TwoAtomFunctions_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("symbol") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::force
        
            typedef ::SireCAS::Expression ( ::SireMM::TwoAtomFunctions::*force_function_type)( ::SireMol::BondID const &,::SireCAS::Symbol const & ) const;
            force_function_type force_function_value( &::SireMM::TwoAtomFunctions::force );
            
            TwoAtomFunctions_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("bondid"), bp::arg("symbol") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::forces
        
            typedef ::QVector< SireMM::TwoAtomFunction > ( ::SireMM::TwoAtomFunctions::*forces_function_type)( ::SireCAS::Symbol const & ) const;
            forces_function_type forces_function_value( &::SireMM::TwoAtomFunctions::forces );
            
            TwoAtomFunctions_exposer.def( 
                "forces"
                , forces_function_value
                , ( bp::arg("symbol") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::includeOnly
        
            typedef ::SireMM::TwoAtomFunctions ( ::SireMM::TwoAtomFunctions::*includeOnly_function_type)( ::SireMol::AtomSelection const &,bool ) const;
            includeOnly_function_type includeOnly_function_value( &::SireMM::TwoAtomFunctions::includeOnly );
            
            TwoAtomFunctions_exposer.def( 
                "includeOnly"
                , includeOnly_function_value
                , ( bp::arg("selection"), bp::arg("isstrict")=(bool)(true) ) );
        
        }
        { //::SireMM::TwoAtomFunctions::isEmpty
        
            typedef bool ( ::SireMM::TwoAtomFunctions::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::TwoAtomFunctions::isEmpty );
            
            TwoAtomFunctions_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMM::TwoAtomFunctions::nFunctions
        
            typedef int ( ::SireMM::TwoAtomFunctions::*nFunctions_function_type)(  ) const;
            nFunctions_function_type nFunctions_function_value( &::SireMM::TwoAtomFunctions::nFunctions );
            
            TwoAtomFunctions_exposer.def( 
                "nFunctions"
                , nFunctions_function_value );
        
        }
        TwoAtomFunctions_exposer.def( bp::self != bp::self );
        { //::SireMM::TwoAtomFunctions::operator=
        
            typedef ::SireMM::TwoAtomFunctions & ( ::SireMM::TwoAtomFunctions::*assign_function_type)( ::SireMM::TwoAtomFunctions const & ) ;
            assign_function_type assign_function_value( &::SireMM::TwoAtomFunctions::operator= );
            
            TwoAtomFunctions_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        TwoAtomFunctions_exposer.def( bp::self == bp::self );
        { //::SireMM::TwoAtomFunctions::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::TwoAtomFunctions::*potential_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx ) const;
            potential_function_type potential_function_value( &::SireMM::TwoAtomFunctions::potential );
            
            TwoAtomFunctions_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::TwoAtomFunctions::*potential_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            potential_function_type potential_function_value( &::SireMM::TwoAtomFunctions::potential );
            
            TwoAtomFunctions_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::TwoAtomFunctions::*potential_function_type)( ::SireMol::BondID const & ) const;
            potential_function_type potential_function_value( &::SireMM::TwoAtomFunctions::potential );
            
            TwoAtomFunctions_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("bondid") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::potentials
        
            typedef ::QVector< SireMM::TwoAtomFunction > ( ::SireMM::TwoAtomFunctions::*potentials_function_type)(  ) const;
            potentials_function_type potentials_function_value( &::SireMM::TwoAtomFunctions::potentials );
            
            TwoAtomFunctions_exposer.def( 
                "potentials"
                , potentials_function_value );
        
        }
        { //::SireMM::TwoAtomFunctions::set
        
            typedef void ( ::SireMM::TwoAtomFunctions::*set_function_type)( ::SireMol::AtomIdx,::SireMol::AtomIdx,::SireCAS::Expression const & ) ;
            set_function_type set_function_value( &::SireMM::TwoAtomFunctions::set );
            
            TwoAtomFunctions_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("expression") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::set
        
            typedef void ( ::SireMM::TwoAtomFunctions::*set_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,::SireCAS::Expression const & ) ;
            set_function_type set_function_value( &::SireMM::TwoAtomFunctions::set );
            
            TwoAtomFunctions_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("expression") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::set
        
            typedef void ( ::SireMM::TwoAtomFunctions::*set_function_type)( ::SireMol::BondID const &,::SireCAS::Expression const & ) ;
            set_function_type set_function_value( &::SireMM::TwoAtomFunctions::set );
            
            TwoAtomFunctions_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("bondid"), bp::arg("expression") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::substitute
        
            typedef void ( ::SireMM::TwoAtomFunctions::*substitute_function_type)( ::SireCAS::Identities const & ) ;
            substitute_function_type substitute_function_value( &::SireMM::TwoAtomFunctions::substitute );
            
            TwoAtomFunctions_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") ) );
        
        }
        { //::SireMM::TwoAtomFunctions::toString
        
            typedef ::QString ( ::SireMM::TwoAtomFunctions::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::TwoAtomFunctions::toString );
            
            TwoAtomFunctions_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMM::TwoAtomFunctions::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::TwoAtomFunctions::typeName );
            
            TwoAtomFunctions_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        TwoAtomFunctions_exposer.staticmethod( "typeName" );
        TwoAtomFunctions_exposer.def( "__copy__", &__copy__);
        TwoAtomFunctions_exposer.def( "__deepcopy__", &__copy__);
        TwoAtomFunctions_exposer.def( "clone", &__copy__);
        TwoAtomFunctions_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::TwoAtomFunctions >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TwoAtomFunctions_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::TwoAtomFunctions >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TwoAtomFunctions_exposer.def( "__str__", &__str__< ::SireMM::TwoAtomFunctions > );
        TwoAtomFunctions_exposer.def( "__repr__", &__str__< ::SireMM::TwoAtomFunctions > );
    }

}
