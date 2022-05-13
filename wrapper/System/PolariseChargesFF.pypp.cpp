// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "PolariseChargesFF.pypp.hpp"

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

SireSystem::PolariseChargesFF __copy__(const SireSystem::PolariseChargesFF &other){ return SireSystem::PolariseChargesFF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_PolariseChargesFF_class(){

    { //::SireSystem::PolariseChargesFF
        typedef bp::class_< SireSystem::PolariseChargesFF, bp::bases< SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > PolariseChargesFF_exposer_t;
        PolariseChargesFF_exposer_t PolariseChargesFF_exposer = PolariseChargesFF_exposer_t( "PolariseChargesFF", "This class implements the forcefield that is used to calculate\nthe self-energy of polarising the charges. This is a companion\nforcefield to the PolariseCharges constraint and is not\ndesigned to be used on its own\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope PolariseChargesFF_scope( PolariseChargesFF_exposer );
        PolariseChargesFF_exposer.def( bp::init< SireSystem::PolariseCharges const & >(( bp::arg("constraint") ), "Construct to calculate the self energy of the molecules affected\nby the passed constraint - note that this forcefield wont notice\nif molecules are added or removed from the constraint, so you must\nmake sure that you add or remove molecules from this forcefield whenever\nyou add or remove molecules from this constraint") );
        PolariseChargesFF_exposer.def( bp::init< QString const &, SireSystem::PolariseCharges const & >(( bp::arg("name"), bp::arg("constraint") ), "Construct to calculate the self energy of the molecules affected\nby the passed constraint - note that this forcefield wont notice\nif molecules are added or removed from the constraint, so you must\nmake sure that you add or remove molecules from this forcefield whenever\nyou add or remove molecules from this constraint") );
        PolariseChargesFF_exposer.def( bp::init< SireSystem::PolariseChargesFF const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::PolariseChargesFF::components
        
            typedef ::SireFF::SingleComponent const & ( ::SireSystem::PolariseChargesFF::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireSystem::PolariseChargesFF::components );
            
            PolariseChargesFF_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the components of this forcefield" );
        
        }
        { //::SireSystem::PolariseChargesFF::containsProperty
        
            typedef bool ( ::SireSystem::PolariseChargesFF::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireSystem::PolariseChargesFF::containsProperty );
            
            PolariseChargesFF_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "This forcefield doesnt contain any properties" );
        
        }
        { //::SireSystem::PolariseChargesFF::mustNowRecalculateFromScratch
        
            typedef void ( ::SireSystem::PolariseChargesFF::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireSystem::PolariseChargesFF::mustNowRecalculateFromScratch );
            
            PolariseChargesFF_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , bp::release_gil_policy()
                , "Tell the forcefield that the energy must now be recalculated\nfrom scratch" );
        
        }
        PolariseChargesFF_exposer.def( bp::self != bp::self );
        { //::SireSystem::PolariseChargesFF::operator=
        
            typedef ::SireSystem::PolariseChargesFF & ( ::SireSystem::PolariseChargesFF::*assign_function_type)( ::SireSystem::PolariseChargesFF const & ) ;
            assign_function_type assign_function_value( &::SireSystem::PolariseChargesFF::operator= );
            
            PolariseChargesFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PolariseChargesFF_exposer.def( bp::self == bp::self );
        { //::SireSystem::PolariseChargesFF::properties
        
            typedef ::SireBase::Properties const & ( ::SireSystem::PolariseChargesFF::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireSystem::PolariseChargesFF::properties );
            
            PolariseChargesFF_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "This forcefield doesnt contain any properties" );
        
        }
        { //::SireSystem::PolariseChargesFF::property
        
            typedef ::SireBase::Property const & ( ::SireSystem::PolariseChargesFF::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireSystem::PolariseChargesFF::property );
            
            PolariseChargesFF_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "This forcefield doesnt contain any properties\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireSystem::PolariseChargesFF::setProperty
        
            typedef bool ( ::SireSystem::PolariseChargesFF::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireSystem::PolariseChargesFF::setProperty );
            
            PolariseChargesFF_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") )
                , bp::release_gil_policy()
                , "You cannot set any properties of this forcefield\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireSystem::PolariseChargesFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::PolariseChargesFF::typeName );
            
            PolariseChargesFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        PolariseChargesFF_exposer.staticmethod( "typeName" );
        PolariseChargesFF_exposer.def( "__copy__", &__copy__);
        PolariseChargesFF_exposer.def( "__deepcopy__", &__copy__);
        PolariseChargesFF_exposer.def( "clone", &__copy__);
        PolariseChargesFF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::PolariseChargesFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PolariseChargesFF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::PolariseChargesFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PolariseChargesFF_exposer.def_pickle(sire_pickle_suite< ::SireSystem::PolariseChargesFF >());
        PolariseChargesFF_exposer.def( "__str__", &__str__< ::SireSystem::PolariseChargesFF > );
        PolariseChargesFF_exposer.def( "__repr__", &__str__< ::SireSystem::PolariseChargesFF > );
        PolariseChargesFF_exposer.def( "__len__", &__len_count< ::SireSystem::PolariseChargesFF > );
    }

}
