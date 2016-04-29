// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "QMFF.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/energytable.h"

#include "SireFF/fieldtable.h"

#include "SireFF/forcetable.h"

#include "SireFF/potentialtable.h"

#include "SireMM/cljprobe.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "qmff.h"

#include "qmprogram.h"

#include <QDebug>

#include "qmff.h"

Squire::QMFF __copy__(const Squire::QMFF &other){ return Squire::QMFF(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_QMFF_class(){

    { //::Squire::QMFF
        typedef bp::class_< Squire::QMFF, bp::bases< SireFF::FF3D, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > QMFF_exposer_t;
        QMFF_exposer_t QMFF_exposer = QMFF_exposer_t( "QMFF", bp::init< >() );
        bp::scope QMFF_scope( QMFF_exposer );
        QMFF_exposer.def( bp::init< QString const & >(( bp::arg("name") )) );
        QMFF_exposer.def( bp::init< Squire::QMFF const & >(( bp::arg("other") )) );
        { //::Squire::QMFF::components
        
            typedef ::Squire::QMComponent const & ( ::Squire::QMFF::*components_function_type)(  ) const;
            components_function_type components_function_value( &::Squire::QMFF::components );
            
            QMFF_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::Squire::QMFF::containsProperty
        
            typedef bool ( ::Squire::QMFF::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::Squire::QMFF::containsProperty );
            
            QMFF_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") ) );
        
        }
        { //::Squire::QMFF::energy
        
            typedef void ( ::Squire::QMFF::*energy_function_type)( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::Squire::QMFF::energy );
            
            QMFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 ) );
        
        }
        { //::Squire::QMFF::energy
        
            typedef void ( ::Squire::QMFF::*energy_function_type)( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::Squire::QMFF::energy );
            
            QMFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 ) );
        
        }
        { //::Squire::QMFF::energyCommandFile
        
            typedef ::QString ( ::Squire::QMFF::*energyCommandFile_function_type)(  ) const;
            energyCommandFile_function_type energyCommandFile_function_value( &::Squire::QMFF::energyCommandFile );
            
            QMFF_exposer.def( 
                "energyCommandFile"
                , energyCommandFile_function_value );
        
        }
        { //::Squire::QMFF::field
        
            typedef void ( ::Squire::QMFF::*field_function_type)( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::Squire::QMFF::field );
            
            QMFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 ) );
        
        }
        { //::Squire::QMFF::field
        
            typedef void ( ::Squire::QMFF::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::Squire::QMFF::field );
            
            QMFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 ) );
        
        }
        { //::Squire::QMFF::field
        
            typedef void ( ::Squire::QMFF::*field_function_type)( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::Squire::QMFF::field );
            
            QMFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::Squire::QMFF::field
        
            typedef void ( ::Squire::QMFF::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::Squire::QMFF::field );
            
            QMFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 ) );
        
        }
        { //::Squire::QMFF::fieldCommandFile
        
            typedef ::QString ( ::Squire::QMFF::*fieldCommandFile_function_type)( ::SireFF::FieldTable const & ) const;
            fieldCommandFile_function_type fieldCommandFile_function_value( &::Squire::QMFF::fieldCommandFile );
            
            QMFF_exposer.def( 
                "fieldCommandFile"
                , fieldCommandFile_function_value
                , ( bp::arg("fieldtable") ) );
        
        }
        { //::Squire::QMFF::fieldCommandFile
        
            typedef ::QString ( ::Squire::QMFF::*fieldCommandFile_function_type)( ::SireFF::FieldTable const &,::SireFF::Probe const & ) const;
            fieldCommandFile_function_type fieldCommandFile_function_value( &::Squire::QMFF::fieldCommandFile );
            
            QMFF_exposer.def( 
                "fieldCommandFile"
                , fieldCommandFile_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe") ) );
        
        }
        { //::Squire::QMFF::force
        
            typedef void ( ::Squire::QMFF::*force_function_type)( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::Squire::QMFF::force );
            
            QMFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 ) );
        
        }
        { //::Squire::QMFF::force
        
            typedef void ( ::Squire::QMFF::*force_function_type)( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::Squire::QMFF::force );
            
            QMFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 ) );
        
        }
        { //::Squire::QMFF::forceCommandFile
        
            typedef ::QString ( ::Squire::QMFF::*forceCommandFile_function_type)( ::SireFF::ForceTable const & ) const;
            forceCommandFile_function_type forceCommandFile_function_value( &::Squire::QMFF::forceCommandFile );
            
            QMFF_exposer.def( 
                "forceCommandFile"
                , forceCommandFile_function_value
                , ( bp::arg("forcetable") ) );
        
        }
        { //::Squire::QMFF::mustNowRecalculateFromScratch
        
            typedef void ( ::Squire::QMFF::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::Squire::QMFF::mustNowRecalculateFromScratch );
            
            QMFF_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value );
        
        }
        QMFF_exposer.def( bp::self != bp::self );
        { //::Squire::QMFF::operator=
        
            typedef ::Squire::QMFF & ( ::Squire::QMFF::*assign_function_type)( ::Squire::QMFF const & ) ;
            assign_function_type assign_function_value( &::Squire::QMFF::operator= );
            
            QMFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        QMFF_exposer.def( bp::self == bp::self );
        { //::Squire::QMFF::parameters
        
            typedef ::SireFF::detail::AtomicParameters3D< SireMol::Element > ( ::Squire::QMFF::*parameters_function_type)(  ) const;
            parameters_function_type parameters_function_value( &::Squire::QMFF::parameters );
            
            QMFF_exposer.def( 
                "parameters"
                , parameters_function_value );
        
        }
        { //::Squire::QMFF::potential
        
            typedef void ( ::Squire::QMFF::*potential_function_type)( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::Squire::QMFF::potential );
            
            QMFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::Squire::QMFF::potential
        
            typedef void ( ::Squire::QMFF::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::Squire::QMFF::potential );
            
            QMFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::Squire::QMFF::potential
        
            typedef void ( ::Squire::QMFF::*potential_function_type)( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::Squire::QMFF::potential );
            
            QMFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::Squire::QMFF::potential
        
            typedef void ( ::Squire::QMFF::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::Squire::QMFF::potential );
            
            QMFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 ) );
        
        }
        { //::Squire::QMFF::potentialCommandFile
        
            typedef ::QString ( ::Squire::QMFF::*potentialCommandFile_function_type)( ::SireFF::PotentialTable const & ) const;
            potentialCommandFile_function_type potentialCommandFile_function_value( &::Squire::QMFF::potentialCommandFile );
            
            QMFF_exposer.def( 
                "potentialCommandFile"
                , potentialCommandFile_function_value
                , ( bp::arg("pottable") ) );
        
        }
        { //::Squire::QMFF::potentialCommandFile
        
            typedef ::QString ( ::Squire::QMFF::*potentialCommandFile_function_type)( ::SireFF::PotentialTable const &,::SireFF::Probe const & ) const;
            potentialCommandFile_function_type potentialCommandFile_function_value( &::Squire::QMFF::potentialCommandFile );
            
            QMFF_exposer.def( 
                "potentialCommandFile"
                , potentialCommandFile_function_value
                , ( bp::arg("pottable"), bp::arg("probe") ) );
        
        }
        { //::Squire::QMFF::properties
        
            typedef ::SireBase::Properties const & ( ::Squire::QMFF::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::Squire::QMFF::properties );
            
            QMFF_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::Squire::QMFF::property
        
            typedef ::SireBase::Property const & ( ::Squire::QMFF::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::Squire::QMFF::property );
            
            QMFF_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::Squire::QMFF::quantumProgram
        
            typedef ::Squire::QMProgram const & ( ::Squire::QMFF::*quantumProgram_function_type)(  ) const;
            quantumProgram_function_type quantumProgram_function_value( &::Squire::QMFF::quantumProgram );
            
            QMFF_exposer.def( 
                "quantumProgram"
                , quantumProgram_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::Squire::QMFF::setProperty
        
            typedef bool ( ::Squire::QMFF::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::Squire::QMFF::setProperty );
            
            QMFF_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") ) );
        
        }
        { //::Squire::QMFF::setQuantumProgram
        
            typedef bool ( ::Squire::QMFF::*setQuantumProgram_function_type)( ::Squire::QMProgram const & ) ;
            setQuantumProgram_function_type setQuantumProgram_function_value( &::Squire::QMFF::setQuantumProgram );
            
            QMFF_exposer.def( 
                "setQuantumProgram"
                , setQuantumProgram_function_value
                , ( bp::arg("qmprog") ) );
        
        }
        { //::Squire::QMFF::setSpace
        
            typedef bool ( ::Squire::QMFF::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::Squire::QMFF::setSpace );
            
            QMFF_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") ) );
        
        }
        { //::Squire::QMFF::setZeroEnergy
        
            typedef bool ( ::Squire::QMFF::*setZeroEnergy_function_type)( ::SireUnits::Dimension::MolarEnergy ) ;
            setZeroEnergy_function_type setZeroEnergy_function_value( &::Squire::QMFF::setZeroEnergy );
            
            QMFF_exposer.def( 
                "setZeroEnergy"
                , setZeroEnergy_function_value
                , ( bp::arg("zero_energy") ) );
        
        }
        { //::Squire::QMFF::space
        
            typedef ::SireVol::Space const & ( ::Squire::QMFF::*space_function_type)(  ) const;
            space_function_type space_function_value( &::Squire::QMFF::space );
            
            QMFF_exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::Squire::QMFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::QMFF::typeName );
            
            QMFF_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::Squire::QMFF::what
        
            typedef char const * ( ::Squire::QMFF::*what_function_type)(  ) const;
            what_function_type what_function_value( &::Squire::QMFF::what );
            
            QMFF_exposer.def( 
                "what"
                , what_function_value );
        
        }
        { //::Squire::QMFF::zeroEnergy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::Squire::QMFF::*zeroEnergy_function_type)(  ) const;
            zeroEnergy_function_type zeroEnergy_function_value( &::Squire::QMFF::zeroEnergy );
            
            QMFF_exposer.def( 
                "zeroEnergy"
                , zeroEnergy_function_value );
        
        }
        QMFF_exposer.staticmethod( "typeName" );
        QMFF_exposer.def( "__copy__", &__copy__);
        QMFF_exposer.def( "__deepcopy__", &__copy__);
        QMFF_exposer.def( "clone", &__copy__);
        QMFF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::QMFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMFF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::QMFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMFF_exposer.def( "__str__", &__str__< ::Squire::QMFF > );
        QMFF_exposer.def( "__repr__", &__str__< ::Squire::QMFF > );
        QMFF_exposer.def( "__len__", &__len_count< ::Squire::QMFF > );
    }

}
