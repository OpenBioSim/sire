// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "CLJProbe.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "cljprobe.h"

#include "cljprobe.h"

SireMM::CLJProbe __copy__(const SireMM::CLJProbe &other){ return SireMM::CLJProbe(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CLJProbe_class(){

    { //::SireMM::CLJProbe
        typedef bp::class_< SireMM::CLJProbe, bp::bases< SireFF::Probe, SireBase::Property > > CLJProbe_exposer_t;
        CLJProbe_exposer_t CLJProbe_exposer = CLJProbe_exposer_t( "CLJProbe", "This is a probe used to probe the coulomb+LJ field\nor potential at points in a forcefield", bp::init< >("Construct a default probe the represents a unit positive charge,\nand an OPLS united atom methane") );
        bp::scope CLJProbe_scope( CLJProbe_exposer );
        CLJProbe_exposer.def( bp::init< SireUnits::Dimension::Charge >(( bp::arg("charge") ), "Construct a probe that is just the passed charge (with dummy LJ parameters)") );
        CLJProbe_exposer.def( bp::init< SireMM::LJParameter const & >(( bp::arg("ljparam") ), "Construct a probe that is just the passed LJ parameter (zero charge)") );
        CLJProbe_exposer.def( bp::init< SireUnits::Dimension::Charge, SireMM::LJParameter const & >(( bp::arg("charge"), bp::arg("ljparam") ), "Construct a probe with the passed charge and LJ parameter") );
        CLJProbe_exposer.def( bp::init< SireMM::CoulombProbe const & >(( bp::arg("probe") ), "Construct to hold just the charge from probe (dummy LJ)") );
        CLJProbe_exposer.def( bp::init< SireMM::LJProbe const & >(( bp::arg("probe") ), "Construct to hold just the LJ parameters from probe (zero charge)") );
        CLJProbe_exposer.def( bp::init< SireFF::Probe const & >(( bp::arg("probe") ), "Construct to hold just the LJ parameters from probe (zero charge)") );
        CLJProbe_exposer.def( bp::init< SireMM::CLJProbe const & >(( bp::arg("cljprobe") ), "Construct to hold just the LJ parameters from probe (zero charge)") );
        { //::SireMM::CLJProbe::charge
        
            typedef ::SireUnits::Dimension::Charge ( ::SireMM::CLJProbe::*charge_function_type)(  ) const;
            charge_function_type charge_function_value( &::SireMM::CLJProbe::charge );
            
            CLJProbe_exposer.def( 
                "charge"
                , charge_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJProbe::lj
        
            typedef ::SireMM::LJParameter const & ( ::SireMM::CLJProbe::*lj_function_type)(  ) const;
            lj_function_type lj_function_value( &::SireMM::CLJProbe::lj );
            
            CLJProbe_exposer.def( 
                "lj"
                , lj_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        CLJProbe_exposer.def( bp::self != bp::self );
        { //::SireMM::CLJProbe::operator=
        
            typedef ::SireMM::CLJProbe & ( ::SireMM::CLJProbe::*assign_function_type)( ::SireMM::CLJProbe const & ) ;
            assign_function_type assign_function_value( &::SireMM::CLJProbe::operator= );
            
            CLJProbe_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CLJProbe_exposer.def( bp::self == bp::self );
        { //::SireMM::CLJProbe::reducedCharge
        
            typedef double ( ::SireMM::CLJProbe::*reducedCharge_function_type)(  ) const;
            reducedCharge_function_type reducedCharge_function_value( &::SireMM::CLJProbe::reducedCharge );
            
            CLJProbe_exposer.def( 
                "reducedCharge"
                , reducedCharge_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJProbe::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJProbe::typeName );
            
            CLJProbe_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CLJProbe_exposer.staticmethod( "typeName" );
        CLJProbe_exposer.def( "__copy__", &__copy__<SireMM::CLJProbe>);
        CLJProbe_exposer.def( "__deepcopy__", &__copy__<SireMM::CLJProbe>);
        CLJProbe_exposer.def( "clone", &__copy__<SireMM::CLJProbe>);
        CLJProbe_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJProbe >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJProbe_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJProbe >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJProbe_exposer.def_pickle(sire_pickle_suite< ::SireMM::CLJProbe >());
        CLJProbe_exposer.def( "__str__", &__str__< ::SireMM::CLJProbe > );
        CLJProbe_exposer.def( "__repr__", &__str__< ::SireMM::CLJProbe > );
    }

}
