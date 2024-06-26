// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "PositionalRestraints.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "positionalrestraints.h"

#include <QDebug>

#include "positionalrestraints.h"

SireMM::PositionalRestraints __copy__(const SireMM::PositionalRestraints &other){ return SireMM::PositionalRestraints(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_PositionalRestraints_class(){

    { //::SireMM::PositionalRestraints
        typedef bp::class_< SireMM::PositionalRestraints, bp::bases< SireMM::Restraints, SireBase::Property > > PositionalRestraints_exposer_t;
        PositionalRestraints_exposer_t PositionalRestraints_exposer = PositionalRestraints_exposer_t( "PositionalRestraints", "This class provides the information for a collection of positional\nrestraints that can be added to a collection of molecues. Each\nrestraint can act on a particle or the centroid of a collection\nof particles. The restaints are spherically symmetric, and\nare either flat-bottom harmonics or harmonic potentials\n", bp::init< >("Null constructor") );
        bp::scope PositionalRestraints_scope( PositionalRestraints_exposer );
        PositionalRestraints_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        PositionalRestraints_exposer.def( bp::init< SireMM::PositionalRestraint const & >(( bp::arg("restraint") ), "") );
        PositionalRestraints_exposer.def( bp::init< QList< SireMM::PositionalRestraint > const & >(( bp::arg("restraints") ), "") );
        PositionalRestraints_exposer.def( bp::init< QString const &, SireMM::PositionalRestraint const & >(( bp::arg("name"), bp::arg("restraint") ), "") );
        PositionalRestraints_exposer.def( bp::init< QString const &, QList< SireMM::PositionalRestraint > const & >(( bp::arg("name"), bp::arg("restraints") ), "") );
        PositionalRestraints_exposer.def( bp::init< SireMM::PositionalRestraints const & >(( bp::arg("other") ), "") );
        { //::SireMM::PositionalRestraints::add
        
            typedef void ( ::SireMM::PositionalRestraints::*add_function_type)( ::SireMM::PositionalRestraint const & ) ;
            add_function_type add_function_value( &::SireMM::PositionalRestraints::add );
            
            PositionalRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraint") )
                , bp::release_gil_policy()
                , "Add a restraint onto the list" );
        
        }
        { //::SireMM::PositionalRestraints::add
        
            typedef void ( ::SireMM::PositionalRestraints::*add_function_type)( ::SireMM::PositionalRestraints const & ) ;
            add_function_type add_function_value( &::SireMM::PositionalRestraints::add );
            
            PositionalRestraints_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("restraints") )
                , bp::release_gil_policy()
                , "Add a restraint onto the list" );
        
        }
        { //::SireMM::PositionalRestraints::at
        
            typedef ::SireMM::PositionalRestraint const & ( ::SireMM::PositionalRestraints::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMM::PositionalRestraints::at );
            
            PositionalRestraints_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ith restraint" );
        
        }
        { //::SireMM::PositionalRestraints::atomRestraints
        
            typedef ::QList< SireMM::PositionalRestraint > ( ::SireMM::PositionalRestraints::*atomRestraints_function_type)(  ) const;
            atomRestraints_function_type atomRestraints_function_value( &::SireMM::PositionalRestraints::atomRestraints );
            
            PositionalRestraints_exposer.def( 
                "atomRestraints"
                , atomRestraints_function_value
                , bp::release_gil_policy()
                , "Return all of the atom restraints" );
        
        }
        { //::SireMM::PositionalRestraints::centroidRestraints
        
            typedef ::QList< SireMM::PositionalRestraint > ( ::SireMM::PositionalRestraints::*centroidRestraints_function_type)(  ) const;
            centroidRestraints_function_type centroidRestraints_function_value( &::SireMM::PositionalRestraints::centroidRestraints );
            
            PositionalRestraints_exposer.def( 
                "centroidRestraints"
                , centroidRestraints_function_value
                , bp::release_gil_policy()
                , "Return all of the centroid restraints" );
        
        }
        { //::SireMM::PositionalRestraints::count
        
            typedef int ( ::SireMM::PositionalRestraints::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::PositionalRestraints::count );
            
            PositionalRestraints_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::PositionalRestraints::hasAtomRestraints
        
            typedef bool ( ::SireMM::PositionalRestraints::*hasAtomRestraints_function_type)(  ) const;
            hasAtomRestraints_function_type hasAtomRestraints_function_value( &::SireMM::PositionalRestraints::hasAtomRestraints );
            
            PositionalRestraints_exposer.def( 
                "hasAtomRestraints"
                , hasAtomRestraints_function_value
                , bp::release_gil_policy()
                , "Return whether or not there are any atom restraints" );
        
        }
        { //::SireMM::PositionalRestraints::hasCentroidRestraints
        
            typedef bool ( ::SireMM::PositionalRestraints::*hasCentroidRestraints_function_type)(  ) const;
            hasCentroidRestraints_function_type hasCentroidRestraints_function_value( &::SireMM::PositionalRestraints::hasCentroidRestraints );
            
            PositionalRestraints_exposer.def( 
                "hasCentroidRestraints"
                , hasCentroidRestraints_function_value
                , bp::release_gil_policy()
                , "Return whether or not there are any centroid restraints" );
        
        }
        { //::SireMM::PositionalRestraints::isEmpty
        
            typedef bool ( ::SireMM::PositionalRestraints::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::PositionalRestraints::isEmpty );
            
            PositionalRestraints_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::PositionalRestraints::isNull
        
            typedef bool ( ::SireMM::PositionalRestraints::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMM::PositionalRestraints::isNull );
            
            PositionalRestraints_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMM::PositionalRestraints::nAtomRestraints
        
            typedef int ( ::SireMM::PositionalRestraints::*nAtomRestraints_function_type)(  ) const;
            nAtomRestraints_function_type nAtomRestraints_function_value( &::SireMM::PositionalRestraints::nAtomRestraints );
            
            PositionalRestraints_exposer.def( 
                "nAtomRestraints"
                , nAtomRestraints_function_value
                , bp::release_gil_policy()
                , "Return the number of atom restraints" );
        
        }
        { //::SireMM::PositionalRestraints::nCentroidRestraints
        
            typedef int ( ::SireMM::PositionalRestraints::*nCentroidRestraints_function_type)(  ) const;
            nCentroidRestraints_function_type nCentroidRestraints_function_value( &::SireMM::PositionalRestraints::nCentroidRestraints );
            
            PositionalRestraints_exposer.def( 
                "nCentroidRestraints"
                , nCentroidRestraints_function_value
                , bp::release_gil_policy()
                , "Return the number of centroid restraints" );
        
        }
        { //::SireMM::PositionalRestraints::nRestraints
        
            typedef int ( ::SireMM::PositionalRestraints::*nRestraints_function_type)(  ) const;
            nRestraints_function_type nRestraints_function_value( &::SireMM::PositionalRestraints::nRestraints );
            
            PositionalRestraints_exposer.def( 
                "nRestraints"
                , nRestraints_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        PositionalRestraints_exposer.def( bp::self != bp::self );
        PositionalRestraints_exposer.def( bp::self + bp::other< SireMM::PositionalRestraint >() );
        PositionalRestraints_exposer.def( bp::self + bp::self );
        { //::SireMM::PositionalRestraints::operator=
        
            typedef ::SireMM::PositionalRestraints & ( ::SireMM::PositionalRestraints::*assign_function_type)( ::SireMM::PositionalRestraints const & ) ;
            assign_function_type assign_function_value( &::SireMM::PositionalRestraints::operator= );
            
            PositionalRestraints_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PositionalRestraints_exposer.def( bp::self == bp::self );
        { //::SireMM::PositionalRestraints::operator[]
        
            typedef ::SireMM::PositionalRestraint const & ( ::SireMM::PositionalRestraints::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::PositionalRestraints::operator[] );
            
            PositionalRestraints_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::PositionalRestraints::restraints
        
            typedef ::QList< SireMM::PositionalRestraint > ( ::SireMM::PositionalRestraints::*restraints_function_type)(  ) const;
            restraints_function_type restraints_function_value( &::SireMM::PositionalRestraints::restraints );
            
            PositionalRestraints_exposer.def( 
                "restraints"
                , restraints_function_value
                , bp::release_gil_policy()
                , "Return all of the restraints" );
        
        }
        { //::SireMM::PositionalRestraints::size
        
            typedef int ( ::SireMM::PositionalRestraints::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::PositionalRestraints::size );
            
            PositionalRestraints_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of restraints" );
        
        }
        { //::SireMM::PositionalRestraints::toString
        
            typedef ::QString ( ::SireMM::PositionalRestraints::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::PositionalRestraints::toString );
            
            PositionalRestraints_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::PositionalRestraints::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::PositionalRestraints::typeName );
            
            PositionalRestraints_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::PositionalRestraints::what
        
            typedef char const * ( ::SireMM::PositionalRestraints::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::PositionalRestraints::what );
            
            PositionalRestraints_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        PositionalRestraints_exposer.staticmethod( "typeName" );
        PositionalRestraints_exposer.def( "__copy__", &__copy__<SireMM::PositionalRestraints>);
        PositionalRestraints_exposer.def( "__deepcopy__", &__copy__<SireMM::PositionalRestraints>);
        PositionalRestraints_exposer.def( "clone", &__copy__<SireMM::PositionalRestraints>);
        PositionalRestraints_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::PositionalRestraints >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PositionalRestraints_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::PositionalRestraints >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PositionalRestraints_exposer.def_pickle(sire_pickle_suite< ::SireMM::PositionalRestraints >());
        PositionalRestraints_exposer.def( "__str__", &__str__< ::SireMM::PositionalRestraints > );
        PositionalRestraints_exposer.def( "__repr__", &__str__< ::SireMM::PositionalRestraints > );
        PositionalRestraints_exposer.def( "__len__", &__len_size< ::SireMM::PositionalRestraints > );
    }

}
