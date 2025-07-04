// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ResidueBeading.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/console.h"

#include "SireBase/errors.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atombeads.h"

#include "atomidx.h"

#include "atomselection.h"

#include "beadidx.h"

#include "beading.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include <boost/noncopyable.hpp>

#include "beading.h"

SireMol::ResidueBeading __copy__(const SireMol::ResidueBeading &other){ return SireMol::ResidueBeading(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_ResidueBeading_class(){

    { //::SireMol::ResidueBeading
        typedef bp::class_< SireMol::ResidueBeading, bp::bases< SireMol::Beading, SireMol::MolViewProperty, SireBase::Property > > ResidueBeading_exposer_t;
        ResidueBeading_exposer_t ResidueBeading_exposer = ResidueBeading_exposer_t( "ResidueBeading", "This is a beading function that breaks a molecule into beads\nwith one bead per residue\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope ResidueBeading_scope( ResidueBeading_exposer );
        ResidueBeading_exposer.def( bp::init< SireMol::ResidueBeading const & >(( bp::arg("other") ), "Copy constructor") );
        ResidueBeading_exposer.def( bp::self != bp::self );
        { //::SireMol::ResidueBeading::operator=
        
            typedef ::SireMol::ResidueBeading & ( ::SireMol::ResidueBeading::*assign_function_type)( ::SireMol::ResidueBeading const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResidueBeading::operator= );
            
            ResidueBeading_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ResidueBeading_exposer.def( bp::self == bp::self );
        { //::SireMol::ResidueBeading::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResidueBeading::typeName );
            
            ResidueBeading_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ResidueBeading_exposer.staticmethod( "typeName" );
        ResidueBeading_exposer.def( "__copy__", &__copy__<SireMol::ResidueBeading>);
        ResidueBeading_exposer.def( "__deepcopy__", &__copy__<SireMol::ResidueBeading>);
        ResidueBeading_exposer.def( "clone", &__copy__<SireMol::ResidueBeading>);
        ResidueBeading_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResidueBeading >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResidueBeading_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResidueBeading >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResidueBeading_exposer.def_pickle(sire_pickle_suite< ::SireMol::ResidueBeading >());
        ResidueBeading_exposer.def( "__str__", &__str__< ::SireMol::ResidueBeading > );
        ResidueBeading_exposer.def( "__repr__", &__str__< ::SireMol::ResidueBeading > );
    }

}
