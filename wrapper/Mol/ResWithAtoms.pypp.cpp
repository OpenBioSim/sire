// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "ResWithAtoms.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "withatoms.h"

#include "withatoms.h"

SireMol::ResWithAtoms __copy__(const SireMol::ResWithAtoms &other){ return SireMol::ResWithAtoms(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_ResWithAtoms_class(){

    { //::SireMol::ResWithAtoms
        typedef bp::class_< SireMol::ResWithAtoms, bp::bases< SireMol::ResID, SireID::ID > > ResWithAtoms_exposer_t;
        ResWithAtoms_exposer_t ResWithAtoms_exposer = ResWithAtoms_exposer_t( "ResWithAtoms", "This ID class identifies residues that contain atoms that\nmatch the passed AtomID\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope ResWithAtoms_scope( ResWithAtoms_exposer );
        ResWithAtoms_exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("atomid") ), "Construct from the passed AtomID") );
        ResWithAtoms_exposer.def( bp::init< SireMol::ResWithAtoms const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::ResWithAtoms::atomID
        
            typedef ::SireMol::AtomID const & ( ::SireMol::ResWithAtoms::*atomID_function_type)(  ) const;
            atomID_function_type atomID_function_value( &::SireMol::ResWithAtoms::atomID );
            
            ResWithAtoms_exposer.def( 
                "atomID"
                , atomID_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the atom ID" );
        
        }
        { //::SireMol::ResWithAtoms::hash
        
            typedef ::uint ( ::SireMol::ResWithAtoms::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::ResWithAtoms::hash );
            
            ResWithAtoms_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "Return a hash of this identifier" );
        
        }
        { //::SireMol::ResWithAtoms::isNull
        
            typedef bool ( ::SireMol::ResWithAtoms::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::ResWithAtoms::isNull );
            
            ResWithAtoms_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "Is this selection null?" );
        
        }
        { //::SireMol::ResWithAtoms::map
        
            typedef ::QList< SireMol::ResIdx > ( ::SireMol::ResWithAtoms::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::ResWithAtoms::map );
            
            ResWithAtoms_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "Map this ID to the list of indicies of residues that match this ID\nThrow: SireMol::missing_atom\nThrow: SireMol::missing_residue\nThrow: SireError::invalid_index\n" );
        
        }
        ResWithAtoms_exposer.def( bp::self != bp::other< SireID::ID >() );
        ResWithAtoms_exposer.def( bp::self != bp::self );
        { //::SireMol::ResWithAtoms::operator=
        
            typedef ::SireMol::ResWithAtoms & ( ::SireMol::ResWithAtoms::*assign_function_type)( ::SireMol::ResWithAtoms const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResWithAtoms::operator= );
            
            ResWithAtoms_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ResWithAtoms_exposer.def( bp::self == bp::other< SireID::ID >() );
        ResWithAtoms_exposer.def( bp::self == bp::self );
        { //::SireMol::ResWithAtoms::toString
        
            typedef ::QString ( ::SireMol::ResWithAtoms::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResWithAtoms::toString );
            
            ResWithAtoms_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representatio of this ID" );
        
        }
        { //::SireMol::ResWithAtoms::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResWithAtoms::typeName );
            
            ResWithAtoms_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResWithAtoms::what
        
            typedef char const * ( ::SireMol::ResWithAtoms::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ResWithAtoms::what );
            
            ResWithAtoms_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ResWithAtoms_exposer.staticmethod( "typeName" );
        ResWithAtoms_exposer.def( "__copy__", &__copy__<SireMol::ResWithAtoms>);
        ResWithAtoms_exposer.def( "__deepcopy__", &__copy__<SireMol::ResWithAtoms>);
        ResWithAtoms_exposer.def( "clone", &__copy__<SireMol::ResWithAtoms>);
        ResWithAtoms_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResWithAtoms >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResWithAtoms_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResWithAtoms >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResWithAtoms_exposer.def_pickle(sire_pickle_suite< ::SireMol::ResWithAtoms >());
        ResWithAtoms_exposer.def( "__str__", &__str__< ::SireMol::ResWithAtoms > );
        ResWithAtoms_exposer.def( "__repr__", &__str__< ::SireMol::ResWithAtoms > );
        ResWithAtoms_exposer.def( "__hash__", &::SireMol::ResWithAtoms::hash );
    }

}
