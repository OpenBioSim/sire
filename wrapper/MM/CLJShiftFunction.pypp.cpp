// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "CLJShiftFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/multidouble.h"

#include "SireMaths/multifloat.h"

#include "SireMaths/multiint.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireVol/gridinfo.h"

#include "cljshiftfunction.h"

#include <QDebug>

#include <QElapsedTimer>

#include "cljshiftfunction.h"

SireMM::CLJShiftFunction __copy__(const SireMM::CLJShiftFunction &other){ return SireMM::CLJShiftFunction(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CLJShiftFunction_class(){

    { //::SireMM::CLJShiftFunction
        typedef bp::class_< SireMM::CLJShiftFunction, bp::bases< SireMM::CLJCutoffFunction, SireMM::CLJFunction, SireBase::Property > > CLJShiftFunction_exposer_t;
        CLJShiftFunction_exposer_t CLJShiftFunction_exposer = CLJShiftFunction_exposer_t( "CLJShiftFunction", "This CLJFunction calculates the intermolecular coulomb and LJ energy of the passed\nCLJAtoms using a force shifted electrostatics cutoff\n\nwe use the force shifted coulomb energy described\nin Fennell and Gezelter, J. Chem. Phys., 124, 234104, 2006\nWe use alpha=0, as I have seen that a 25 A cutoff gives stable results\nwith alpha=0, and this way we avoid changing the hamiltonian significantly\nby having an erfc function\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope CLJShiftFunction_scope( CLJShiftFunction_exposer );
        CLJShiftFunction_exposer.def( bp::init< SireUnits::Dimension::Length >(( bp::arg("cutoff") ), "Copy constructor") );
        CLJShiftFunction_exposer.def( bp::init< SireUnits::Dimension::Length, SireUnits::Dimension::Length >(( bp::arg("coul_cutoff"), bp::arg("lj_cutoff") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length >(( bp::arg("space"), bp::arg("cutoff") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length, SireUnits::Dimension::Length >(( bp::arg("space"), bp::arg("coul_cutoff"), bp::arg("lj_cutoff") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("cutoff"), bp::arg("combining_rules") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireUnits::Dimension::Length, SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("coul_cutoff"), bp::arg("lj_cutoff"), bp::arg("combining_rules") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireVol::Space const &, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("space"), bp::arg("combining_rules") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("space"), bp::arg("cutoff"), bp::arg("combining_rules") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length, SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("space"), bp::arg("coul_cutoff"), bp::arg("lj_cutoff"), bp::arg("combining_rules") ), "") );
        CLJShiftFunction_exposer.def( bp::init< SireMM::CLJShiftFunction const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::CLJShiftFunction::defaultShiftFunction
        
            typedef ::SireMM::CLJFunctionPtr ( *defaultShiftFunction_function_type )(  );
            defaultShiftFunction_function_type defaultShiftFunction_function_value( &::SireMM::CLJShiftFunction::defaultShiftFunction );
            
            CLJShiftFunction_exposer.def( 
                "defaultShiftFunction"
                , defaultShiftFunction_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CLJShiftFunction_exposer.def( bp::self != bp::self );
        { //::SireMM::CLJShiftFunction::operator=
        
            typedef ::SireMM::CLJShiftFunction & ( ::SireMM::CLJShiftFunction::*assign_function_type)( ::SireMM::CLJShiftFunction const & ) ;
            assign_function_type assign_function_value( &::SireMM::CLJShiftFunction::operator= );
            
            CLJShiftFunction_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CLJShiftFunction_exposer.def( bp::self == bp::self );
        { //::SireMM::CLJShiftFunction::supportsGridCalculation
        
            typedef bool ( ::SireMM::CLJShiftFunction::*supportsGridCalculation_function_type)(  ) const;
            supportsGridCalculation_function_type supportsGridCalculation_function_value( &::SireMM::CLJShiftFunction::supportsGridCalculation );
            
            CLJShiftFunction_exposer.def( 
                "supportsGridCalculation"
                , supportsGridCalculation_function_value
                , bp::release_gil_policy()
                , "This function does support calculations using a grid" );
        
        }
        { //::SireMM::CLJShiftFunction::toString
        
            typedef ::QString ( ::SireMM::CLJShiftFunction::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::CLJShiftFunction::toString );
            
            CLJShiftFunction_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJShiftFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJShiftFunction::typeName );
            
            CLJShiftFunction_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJShiftFunction::what
        
            typedef char const * ( ::SireMM::CLJShiftFunction::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::CLJShiftFunction::what );
            
            CLJShiftFunction_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CLJShiftFunction_exposer.staticmethod( "defaultShiftFunction" );
        CLJShiftFunction_exposer.staticmethod( "typeName" );
        CLJShiftFunction_exposer.def( "__copy__", &__copy__<SireMM::CLJShiftFunction>);
        CLJShiftFunction_exposer.def( "__deepcopy__", &__copy__<SireMM::CLJShiftFunction>);
        CLJShiftFunction_exposer.def( "clone", &__copy__<SireMM::CLJShiftFunction>);
        CLJShiftFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJShiftFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJShiftFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJShiftFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJShiftFunction_exposer.def_pickle(sire_pickle_suite< ::SireMM::CLJShiftFunction >());
        CLJShiftFunction_exposer.def( "__str__", &__str__< ::SireMM::CLJShiftFunction > );
        CLJShiftFunction_exposer.def( "__repr__", &__str__< ::SireMM::CLJShiftFunction > );
    }

}
