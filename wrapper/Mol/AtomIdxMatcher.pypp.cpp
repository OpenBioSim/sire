// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomIdxMatcher.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "atom.h"

#include "atomidentifier.h"

#include "atomidx.h"

#include "atommatcher.h"

#include "atommatchers.h"

#include "atomname.h"

#include "atomselection.h"

#include "evaluator.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "mover.h"

#include "selector.hpp"

#include "tostring.h"

#include "atommatchers.h"

SireMol::AtomIdxMatcher __copy__(const SireMol::AtomIdxMatcher &other){ return SireMol::AtomIdxMatcher(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomIdxMatcher_class(){

    { //::SireMol::AtomIdxMatcher
        typedef bp::class_< SireMol::AtomIdxMatcher, bp::bases< SireMol::AtomMatcher, SireBase::Property > > AtomIdxMatcher_exposer_t;
        AtomIdxMatcher_exposer_t AtomIdxMatcher_exposer = AtomIdxMatcher_exposer_t( "AtomIdxMatcher", "This is a simple atom matcher that matches the atoms based\non their index in the molecule - e.g. it matches the first\natom in molinfo0 to the first atom in molinfo1, the second\natom in molinfo0 to the second atom in molinfo1, and the\nnth atom in molinfo0 to the nth atom in molinfo1\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope AtomIdxMatcher_scope( AtomIdxMatcher_exposer );
        AtomIdxMatcher_exposer.def( bp::init< SireMol::AtomIdxMatcher const & >(( bp::arg("arg0") ), "Copy constructor") );
        AtomIdxMatcher_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomIdxMatcher::operator=
        
            typedef ::SireMol::AtomIdxMatcher & ( ::SireMol::AtomIdxMatcher::*assign_function_type)( ::SireMol::AtomIdxMatcher const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomIdxMatcher::operator= );
            
            AtomIdxMatcher_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomIdxMatcher_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomIdxMatcher::toString
        
            typedef ::QString ( ::SireMol::AtomIdxMatcher::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomIdxMatcher::toString );
            
            AtomIdxMatcher_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::AtomIdxMatcher::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomIdxMatcher::typeName );
            
            AtomIdxMatcher_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::AtomIdxMatcher::what
        
            typedef char const * ( ::SireMol::AtomIdxMatcher::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::AtomIdxMatcher::what );
            
            AtomIdxMatcher_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        AtomIdxMatcher_exposer.staticmethod( "typeName" );
        AtomIdxMatcher_exposer.def( "__copy__", &__copy__);
        AtomIdxMatcher_exposer.def( "__deepcopy__", &__copy__);
        AtomIdxMatcher_exposer.def( "clone", &__copy__);
        AtomIdxMatcher_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomIdxMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIdxMatcher_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomIdxMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIdxMatcher_exposer.def( "__str__", &__str__< ::SireMol::AtomIdxMatcher > );
        AtomIdxMatcher_exposer.def( "__repr__", &__str__< ::SireMol::AtomIdxMatcher > );
    }

}
