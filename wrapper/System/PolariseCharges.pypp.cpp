// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "PolariseCharges.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/refcountdata.h"

#include "SireError/errors.h"

#include "SireFF/potentialtable.h"

#include "SireFF/probe.h"

#include "SireMM/cljprobe.h"

#include "SireMaths/nmatrix.h"

#include "SireMaths/nvector.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomenergies.h"

#include "SireMol/atompolarisabilities.h"

#include "SireMol/atomselection.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculedata.h"

#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/convert.h"

#include "SireUnits/units.h"

#include "delta.h"

#include "polarisecharges.h"

#include "polarisecharges.h"

SireSystem::PolariseCharges __copy__(const SireSystem::PolariseCharges &other){ return SireSystem::PolariseCharges(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_PolariseCharges_class(){

    { //::SireSystem::PolariseCharges
        typedef bp::class_< SireSystem::PolariseCharges, bp::bases< SireSystem::ChargeConstraint, SireSystem::MoleculeConstraint, SireSystem::Constraint, SireBase::Property > > PolariseCharges_exposer_t;
        PolariseCharges_exposer_t PolariseCharges_exposer = PolariseCharges_exposer_t( "PolariseCharges", "This charge constraint adjusts the partial charges of contained\nmolecules to give the impression that the molecule contains\npolarisable dipoles. This is based on the method developed\nby Reynolds et al. (see ...)\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope PolariseCharges_scope( PolariseCharges_exposer );
        PolariseCharges_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a constraint that uses the total energy field and a\nsingle unit charge probe to polarise the molecules in molgroup") );
        PolariseCharges_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireFF::Probe const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("probe"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a constraint that uses the total energy field and\nthe passed probe to polarise the molecules in molgroup") );
        PolariseCharges_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireCAS::Symbol const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("field_component"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a constraint that uses the field represented by field_component\nand a single unit charge to polarise the molecules in molgroup") );
        PolariseCharges_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireCAS::Symbol const &, SireFF::Probe const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("field_component"), bp::arg("probe"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a constraint that uses the field represented by field_component\nand the passed probe to polarise the molecules in molgroup") );
        PolariseCharges_exposer.def( bp::init< SireSystem::PolariseCharges const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::PolariseCharges::convergenceLimit
        
            typedef double ( ::SireSystem::PolariseCharges::*convergenceLimit_function_type)(  ) const;
            convergenceLimit_function_type convergenceLimit_function_value( &::SireSystem::PolariseCharges::convergenceLimit );
            
            PolariseCharges_exposer.def( 
                "convergenceLimit"
                , convergenceLimit_function_value
                , bp::release_gil_policy()
                , "Return the convergence limit of the calculation" );
        
        }
        { //::SireSystem::PolariseCharges::fieldComponent
        
            typedef ::SireCAS::Symbol const & ( ::SireSystem::PolariseCharges::*fieldComponent_function_type)(  ) const;
            fieldComponent_function_type fieldComponent_function_value( &::SireSystem::PolariseCharges::fieldComponent );
            
            PolariseCharges_exposer.def( 
                "fieldComponent"
                , fieldComponent_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the component of the forcefield that is used to\ncalculate the electrostatic field on the atoms to be\npolarised" );
        
        }
        PolariseCharges_exposer.def( bp::self != bp::self );
        { //::SireSystem::PolariseCharges::operator=
        
            typedef ::SireSystem::PolariseCharges & ( ::SireSystem::PolariseCharges::*assign_function_type)( ::SireSystem::PolariseCharges const & ) ;
            assign_function_type assign_function_value( &::SireSystem::PolariseCharges::operator= );
            
            PolariseCharges_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PolariseCharges_exposer.def( bp::self == bp::self );
        { //::SireSystem::PolariseCharges::probe
        
            typedef ::SireMM::CoulombProbe const & ( ::SireSystem::PolariseCharges::*probe_function_type)(  ) const;
            probe_function_type probe_function_value( &::SireSystem::PolariseCharges::probe );
            
            PolariseCharges_exposer.def( 
                "probe"
                , probe_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the probe that is used to calculate the electrostatic\nfield on the atoms to be polarised" );
        
        }
        { //::SireSystem::PolariseCharges::selfEnergyFF
        
            typedef ::SireSystem::PolariseChargesFF ( ::SireSystem::PolariseCharges::*selfEnergyFF_function_type)(  ) const;
            selfEnergyFF_function_type selfEnergyFF_function_value( &::SireSystem::PolariseCharges::selfEnergyFF );
            
            PolariseCharges_exposer.def( 
                "selfEnergyFF"
                , selfEnergyFF_function_value
                , bp::release_gil_policy()
                , "Return the forcefield that is used to calculate the self-energy of\npolarising the charges. This must be added to any system to which\nthis constraint is applied, as maintaining the constraint\n(by polarising the charges) costs energy, which must be part\nof the system Hamiltonian" );
        
        }
        { //::SireSystem::PolariseCharges::setConvergenceLimit
        
            typedef void ( ::SireSystem::PolariseCharges::*setConvergenceLimit_function_type)( double ) ;
            setConvergenceLimit_function_type setConvergenceLimit_function_value( &::SireSystem::PolariseCharges::setConvergenceLimit );
            
            PolariseCharges_exposer.def( 
                "setConvergenceLimit"
                , setConvergenceLimit_function_value
                , ( bp::arg("limit") )
                , bp::release_gil_policy()
                , "Set the convergence limit of the calculation" );
        
        }
        { //::SireSystem::PolariseCharges::toString
        
            typedef ::QString ( ::SireSystem::PolariseCharges::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireSystem::PolariseCharges::toString );
            
            PolariseCharges_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this constraint" );
        
        }
        { //::SireSystem::PolariseCharges::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::PolariseCharges::typeName );
            
            PolariseCharges_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        PolariseCharges_exposer.staticmethod( "typeName" );
        PolariseCharges_exposer.def( "__copy__", &__copy__<SireSystem::PolariseCharges>);
        PolariseCharges_exposer.def( "__deepcopy__", &__copy__<SireSystem::PolariseCharges>);
        PolariseCharges_exposer.def( "clone", &__copy__<SireSystem::PolariseCharges>);
        PolariseCharges_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::PolariseCharges >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PolariseCharges_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::PolariseCharges >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PolariseCharges_exposer.def_pickle(sire_pickle_suite< ::SireSystem::PolariseCharges >());
        PolariseCharges_exposer.def( "__str__", &__str__< ::SireSystem::PolariseCharges > );
        PolariseCharges_exposer.def( "__repr__", &__str__< ::SireSystem::PolariseCharges > );
    }

}
