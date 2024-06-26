// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "LJNBPairs.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireMol/moleculeinfo.h"

#include "SireStream/datastream.h"

#include "cljnbpairs.h"

#include "cljnbpairs.h"

SireMM::LJNBPairs __copy__(const SireMM::LJNBPairs &other){ return SireMM::LJNBPairs(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_LJNBPairs_class(){

    { //::SireMM::LJNBPairs
        typedef bp::class_< SireMM::LJNBPairs, bp::bases< SireMM::AtomPairs<SireMM::LJScaleFactor>, SireMol::MoleculeProperty, SireMol::MolViewProperty, SireBase::Property > > LJNBPairs_exposer_t;
        LJNBPairs_exposer_t LJNBPairs_exposer = LJNBPairs_exposer_t( "LJNBPairs", "This class holds all of the non-bonded scale factors that are used\nto scale the intramolecular atom-atom Lennard-Jones\ninteractions between pairs of atoms, e.g. for most MM forcefields,\nthe scale factors for 1-1, 1-2 and 1-3 pairs are zero, the\n1-4 pairs are scaled by a LJ factor (e.g. 0.5 for OPLS)\nand the 1-5 and above pairs are not scaled (i.e. the factors equal 1)\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope LJNBPairs_scope( LJNBPairs_exposer );
        LJNBPairs_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireMM::LJScaleFactor const & > >(( bp::arg("molview"), bp::arg("default_scale")=SireMM::LJScaleFactor(1) ), "Construct for the molecule viewed in molview") );
        LJNBPairs_exposer.def( bp::init< SireMol::MoleculeInfoData const &, bp::optional< SireMM::LJScaleFactor const & > >(( bp::arg("molinfo"), bp::arg("default_scale")=SireMM::LJScaleFactor(1) ), "Construct, using default_scale for all of the atom-atom\ninteractions in the molecule molinfo") );
        LJNBPairs_exposer.def( bp::init< SireMM::CLJNBPairs const & >(( bp::arg("cljpairs") ), "Construct from the LJ scaling factors in cljpairs") );
        LJNBPairs_exposer.def( bp::init< SireMM::LJNBPairs const & >(( bp::arg("other") ), "Copy constructor") );
        LJNBPairs_exposer.def( bp::self != bp::self );
        { //::SireMM::LJNBPairs::operator=
        
            typedef ::SireMM::LJNBPairs & ( ::SireMM::LJNBPairs::*assign_function_type)( ::SireMM::LJNBPairs const & ) ;
            assign_function_type assign_function_value( &::SireMM::LJNBPairs::operator= );
            
            LJNBPairs_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMM::LJNBPairs::operator=
        
            typedef ::SireMM::LJNBPairs & ( ::SireMM::LJNBPairs::*assign_function_type)( ::SireMM::CLJNBPairs const & ) ;
            assign_function_type assign_function_value( &::SireMM::LJNBPairs::operator= );
            
            LJNBPairs_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("cljpairs") )
                , bp::return_self< >()
                , "" );
        
        }
        LJNBPairs_exposer.def( bp::self == bp::self );
        { //::SireMM::LJNBPairs::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::LJNBPairs::typeName );
            
            LJNBPairs_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        LJNBPairs_exposer.staticmethod( "typeName" );
        LJNBPairs_exposer.def( "__copy__", &__copy__<SireMM::LJNBPairs>);
        LJNBPairs_exposer.def( "__deepcopy__", &__copy__<SireMM::LJNBPairs>);
        LJNBPairs_exposer.def( "clone", &__copy__<SireMM::LJNBPairs>);
        LJNBPairs_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::LJNBPairs >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LJNBPairs_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::LJNBPairs >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LJNBPairs_exposer.def_pickle(sire_pickle_suite< ::SireMM::LJNBPairs >());
        LJNBPairs_exposer.def( "__str__", &__str__< ::SireMM::LJNBPairs > );
        LJNBPairs_exposer.def( "__repr__", &__str__< ::SireMM::LJNBPairs > );
    }

}
