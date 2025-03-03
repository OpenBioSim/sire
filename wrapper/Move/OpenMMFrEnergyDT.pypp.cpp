// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "OpenMMFrEnergyDT.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/variantproperty.h"

#include "SireFF/forcetable.h"

#include "SireIO/amber.h"

#include "SireMM/atomljs.h"

#include "SireMM/internalff.h"

#include "SireMaths/constants.h"

#include "SireMaths/rangenerator.h"

#include "SireMaths/vector.h"

#include "SireMol/amberparameters.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atommasses.h"

#include "SireMol/bondid.h"

#include "SireMol/connectivity.h"

#include "SireMol/core.h"

#include "SireMol/mgname.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/partialmolecule.h"

#include "SireMove/flexibility.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/convert.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "ensemble.h"

#include "openmmfrenergydt.h"

#include <QDebug>

#include <QTime>

#include <iostream>

#include "openmmfrenergydt.h"

SireMove::OpenMMFrEnergyDT __copy__(const SireMove::OpenMMFrEnergyDT &other){ return SireMove::OpenMMFrEnergyDT(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_OpenMMFrEnergyDT_class(){

    { //::SireMove::OpenMMFrEnergyDT
        typedef bp::class_< SireMove::OpenMMFrEnergyDT, bp::bases< SireMove::Integrator, SireBase::Property > > OpenMMFrEnergyDT_exposer_t;
        OpenMMFrEnergyDT_exposer_t OpenMMFrEnergyDT_exposer = OpenMMFrEnergyDT_exposer_t( "OpenMMFrEnergyDT", "This class implements a free energy methods Using OpenMM.\n\nAuthor: Julien Michel and Gaetano Calabro\n", bp::init< bp::optional< bool > >(( bp::arg("frequent_save_velocities")=(bool)(false) ), "Constructor") );
        bp::scope OpenMMFrEnergyDT_scope( OpenMMFrEnergyDT_exposer );
        OpenMMFrEnergyDT_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireMol::MoleculeGroup const &, bp::optional< bool > >(( bp::arg("molecule_group"), bp::arg("solute_group"), bp::arg("frequent_save_velocities")=(bool)(false) ), "Constructor using the passed molecule groups") );
        OpenMMFrEnergyDT_exposer.def( bp::init< SireMove::OpenMMFrEnergyDT const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::OpenMMFrEnergyDT::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::OpenMMFrEnergyDT::*createWorkspace_function_type)( ::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::OpenMMFrEnergyDT::createWorkspace );
            
            OpenMMFrEnergyDT_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("map")=SireBase::PropertyMap() )
                , "Create an empty workspace" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::OpenMMFrEnergyDT::*createWorkspace_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::OpenMMFrEnergyDT::createWorkspace );
            
            OpenMMFrEnergyDT_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , "Create a workspace for this integrator for the molecule group molgroup" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::ensemble
        
            typedef ::SireMove::Ensemble ( ::SireMove::OpenMMFrEnergyDT::*ensemble_function_type)(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::OpenMMFrEnergyDT::ensemble );
            
            OpenMMFrEnergyDT_exposer.def( 
                "ensemble"
                , ensemble_function_value
                , bp::release_gil_policy()
                , "Return the ensemble of this integrator" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getAlchemical_value
        
            typedef double ( ::SireMove::OpenMMFrEnergyDT::*getAlchemical_value_function_type)(  ) ;
            getAlchemical_value_function_type getAlchemical_value_function_value( &::SireMove::OpenMMFrEnergyDT::getAlchemical_value );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getAlchemical_value"
                , getAlchemical_value_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getAndersen
        
            typedef bool ( ::SireMove::OpenMMFrEnergyDT::*getAndersen_function_type)(  ) ;
            getAndersen_function_type getAndersen_function_value( &::SireMove::OpenMMFrEnergyDT::getAndersen );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getAndersen"
                , getAndersen_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getAndersen_frequency
        
            typedef double ( ::SireMove::OpenMMFrEnergyDT::*getAndersen_frequency_function_type)(  ) ;
            getAndersen_frequency_function_type getAndersen_frequency_function_value( &::SireMove::OpenMMFrEnergyDT::getAndersen_frequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getAndersen_frequency"
                , getAndersen_frequency_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getBufferCoords
        
            typedef bool ( ::SireMove::OpenMMFrEnergyDT::*getBufferCoords_function_type)(  ) ;
            getBufferCoords_function_type getBufferCoords_function_value( &::SireMove::OpenMMFrEnergyDT::getBufferCoords );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getBufferCoords"
                , getBufferCoords_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getCMMremoval_frequency
        
            typedef int ( ::SireMove::OpenMMFrEnergyDT::*getCMMremoval_frequency_function_type)(  ) ;
            getCMMremoval_frequency_function_type getCMMremoval_frequency_function_value( &::SireMove::OpenMMFrEnergyDT::getCMMremoval_frequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getCMMremoval_frequency"
                , getCMMremoval_frequency_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getConstraintType
        
            typedef ::QString ( ::SireMove::OpenMMFrEnergyDT::*getConstraintType_function_type)(  ) ;
            getConstraintType_function_type getConstraintType_function_value( &::SireMove::OpenMMFrEnergyDT::getConstraintType );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getConstraintType"
                , getConstraintType_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getCoulomb_power
        
            typedef int ( ::SireMove::OpenMMFrEnergyDT::*getCoulomb_power_function_type)(  ) ;
            getCoulomb_power_function_type getCoulomb_power_function_value( &::SireMove::OpenMMFrEnergyDT::getCoulomb_power );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getCoulomb_power"
                , getCoulomb_power_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getCutoffType
        
            typedef ::QString ( ::SireMove::OpenMMFrEnergyDT::*getCutoffType_function_type)(  ) ;
            getCutoffType_function_type getCutoffType_function_value( &::SireMove::OpenMMFrEnergyDT::getCutoffType );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getCutoffType"
                , getCutoffType_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getCutoff_distance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMove::OpenMMFrEnergyDT::*getCutoff_distance_function_type)(  ) ;
            getCutoff_distance_function_type getCutoff_distance_function_value( &::SireMove::OpenMMFrEnergyDT::getCutoff_distance );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getCutoff_distance"
                , getCutoff_distance_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getDeltaAlchemical
        
            typedef double ( ::SireMove::OpenMMFrEnergyDT::*getDeltaAlchemical_function_type)(  ) ;
            getDeltaAlchemical_function_type getDeltaAlchemical_function_value( &::SireMove::OpenMMFrEnergyDT::getDeltaAlchemical );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getDeltaAlchemical"
                , getDeltaAlchemical_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getDeviceIndex
        
            typedef ::QString ( ::SireMove::OpenMMFrEnergyDT::*getDeviceIndex_function_type)(  ) ;
            getDeviceIndex_function_type getDeviceIndex_function_value( &::SireMove::OpenMMFrEnergyDT::getDeviceIndex );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getDeviceIndex"
                , getDeviceIndex_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getEnergyFrequency
        
            typedef int ( ::SireMove::OpenMMFrEnergyDT::*getEnergyFrequency_function_type)(  ) ;
            getEnergyFrequency_function_type getEnergyFrequency_function_value( &::SireMove::OpenMMFrEnergyDT::getEnergyFrequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getEnergyFrequency"
                , getEnergyFrequency_function_value
                , bp::release_gil_policy()
                , "Get the frequency of buffering coordinates" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getField_dielectric
        
            typedef double ( ::SireMove::OpenMMFrEnergyDT::*getField_dielectric_function_type)(  ) ;
            getField_dielectric_function_type getField_dielectric_function_value( &::SireMove::OpenMMFrEnergyDT::getField_dielectric );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getField_dielectric"
                , getField_dielectric_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getGradients
        
            typedef ::QVector< double > ( ::SireMove::OpenMMFrEnergyDT::*getGradients_function_type)(  ) ;
            getGradients_function_type getGradients_function_value( &::SireMove::OpenMMFrEnergyDT::getGradients );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getGradients"
                , getGradients_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getMCBarostat
        
            typedef bool ( ::SireMove::OpenMMFrEnergyDT::*getMCBarostat_function_type)(  ) ;
            getMCBarostat_function_type getMCBarostat_function_value( &::SireMove::OpenMMFrEnergyDT::getMCBarostat );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getMCBarostat"
                , getMCBarostat_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getMCBarostat_membrane

            typedef bool ( ::SireMove::OpenMMFrEnergyDT::*getMCBarostat_membrane_function_type)(  ) ;
            getMCBarostat_membrane_function_type getMCBarostat_membrane_function_value( &::SireMove::OpenMMFrEnergyDT::getMCBarostat_membrane );

            OpenMMFrEnergyDT_exposer.def(
                "getMCBarostat_membrane"
                , getMCBarostat_membrane_function_value
                , "" );
        }
        { //::SireMove::OpenMMFrEnergyDT::getMCBarostat_frequency
        
            typedef int ( ::SireMove::OpenMMFrEnergyDT::*getMCBarostat_frequency_function_type)(  ) ;
            getMCBarostat_frequency_function_type getMCBarostat_frequency_function_value( &::SireMove::OpenMMFrEnergyDT::getMCBarostat_frequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getMCBarostat_frequency"
                , getMCBarostat_frequency_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getPlatform
        
            typedef ::QString ( ::SireMove::OpenMMFrEnergyDT::*getPlatform_function_type)(  ) ;
            getPlatform_function_type getPlatform_function_value( &::SireMove::OpenMMFrEnergyDT::getPlatform );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getPlatform"
                , getPlatform_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getPressure
        
            typedef ::SireUnits::Dimension::Pressure ( ::SireMove::OpenMMFrEnergyDT::*getPressure_function_type)(  ) ;
            getPressure_function_type getPressure_function_value( &::SireMove::OpenMMFrEnergyDT::getPressure );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getPressure"
                , getPressure_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getRestraint
        
            typedef bool ( ::SireMove::OpenMMFrEnergyDT::*getRestraint_function_type)(  ) ;
            getRestraint_function_type getRestraint_function_value( &::SireMove::OpenMMFrEnergyDT::getRestraint );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getRestraint"
                , getRestraint_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getShift_delta
        
            typedef double ( ::SireMove::OpenMMFrEnergyDT::*getShift_delta_function_type)(  ) ;
            getShift_delta_function_type getShift_delta_function_value( &::SireMove::OpenMMFrEnergyDT::getShift_delta );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getShift_delta"
                , getShift_delta_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::getTemperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireMove::OpenMMFrEnergyDT::*getTemperature_function_type)(  ) ;
            getTemperature_function_type getTemperature_function_value( &::SireMove::OpenMMFrEnergyDT::getTemperature );
            
            OpenMMFrEnergyDT_exposer.def( 
                "getTemperature"
                , getTemperature_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::initialise
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*initialise_function_type)(  ) ;
            initialise_function_type initialise_function_value( &::SireMove::OpenMMFrEnergyDT::initialise );
            
            OpenMMFrEnergyDT_exposer.def( 
                "initialise"
                , initialise_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::integrate
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*integrate_function_type)( ::SireMove::IntegratorWorkspace &,::SireCAS::Symbol const &,::SireUnits::Dimension::Time,int,bool ) ;
            integrate_function_type integrate_function_value( &::SireMove::OpenMMFrEnergyDT::integrate );
            
            OpenMMFrEnergyDT_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("workspace"), bp::arg("nrg_component"), bp::arg("timestep"), bp::arg("nmoves"), bp::arg("record_stats") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::isTimeReversible
        
            typedef bool ( ::SireMove::OpenMMFrEnergyDT::*isTimeReversible_function_type)(  ) const;
            isTimeReversible_function_type isTimeReversible_function_value( &::SireMove::OpenMMFrEnergyDT::isTimeReversible );
            
            OpenMMFrEnergyDT_exposer.def( 
                "isTimeReversible"
                , isTimeReversible_function_value
                , bp::release_gil_policy()
                , "Return whether or not this integrator is time-reversible" );
        
        }
        OpenMMFrEnergyDT_exposer.def( bp::self != bp::self );
        { //::SireMove::OpenMMFrEnergyDT::operator=
        
            typedef ::SireMove::OpenMMFrEnergyDT & ( ::SireMove::OpenMMFrEnergyDT::*assign_function_type)( ::SireMove::OpenMMFrEnergyDT const & ) ;
            assign_function_type assign_function_value( &::SireMove::OpenMMFrEnergyDT::operator= );
            
            OpenMMFrEnergyDT_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        OpenMMFrEnergyDT_exposer.def( bp::self == bp::self );
        { //::SireMove::OpenMMFrEnergyDT::setAlchemical_value
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setAlchemical_value_function_type)( double ) ;
            setAlchemical_value_function_type setAlchemical_value_function_value( &::SireMove::OpenMMFrEnergyDT::setAlchemical_value );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setAlchemical_value"
                , setAlchemical_value_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the alchemical value used to calculate the free energy change via TI method" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setAndersen
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setAndersen_function_type)( bool ) ;
            setAndersen_function_type setAndersen_function_value( &::SireMove::OpenMMFrEnergyDT::setAndersen );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setAndersen"
                , setAndersen_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set Andersen thermostat" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setAndersen_frequency
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setAndersen_frequency_function_type)( double ) ;
            setAndersen_frequency_function_type setAndersen_frequency_function_value( &::SireMove::OpenMMFrEnergyDT::setAndersen_frequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setAndersen_frequency"
                , setAndersen_frequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Andersen Thermostat frequency collision" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setBufferCoords
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setBufferCoords_function_type)( bool ) ;
            setBufferCoords_function_type setBufferCoords_function_value( &::SireMove::OpenMMFrEnergyDT::setBufferCoords );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setBufferCoords"
                , setBufferCoords_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the flag to buffer coordinate during the free energy calculation" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setCMMremoval_frequency
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setCMMremoval_frequency_function_type)( int ) ;
            setCMMremoval_frequency_function_type setCMMremoval_frequency_function_value( &::SireMove::OpenMMFrEnergyDT::setCMMremoval_frequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setCMMremoval_frequency"
                , setCMMremoval_frequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Center of Mass motion removal frequency" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setConstraintType
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setConstraintType_function_type)( ::QString ) ;
            setConstraintType_function_type setConstraintType_function_value( &::SireMove::OpenMMFrEnergyDT::setConstraintType );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setConstraintType"
                , setConstraintType_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Constraint type: none, hbonds, allbonds, hangles" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setCoulomb_power
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setCoulomb_power_function_type)( int ) ;
            setCoulomb_power_function_type setCoulomb_power_function_value( &::SireMove::OpenMMFrEnergyDT::setCoulomb_power );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setCoulomb_power"
                , setCoulomb_power_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the coulomb power used in the soft core potential" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setCutoffType
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setCutoffType_function_type)( ::QString ) ;
            setCutoffType_function_type setCutoffType_function_value( &::SireMove::OpenMMFrEnergyDT::setCutoffType );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setCutoffType"
                , setCutoffType_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the cufott type: nocutoff, cutoffnonperiodic, cutoffperiodic" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setCutoff_distance
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setCutoff_distance_function_type)( ::SireUnits::Dimension::Length ) ;
            setCutoff_distance_function_type setCutoff_distance_function_value( &::SireMove::OpenMMFrEnergyDT::setCutoff_distance );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setCutoff_distance"
                , setCutoff_distance_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the cutoff distance in A" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setDeltatAlchemical
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setDeltatAlchemical_function_type)( double ) ;
            setDeltatAlchemical_function_type setDeltatAlchemical_function_value( &::SireMove::OpenMMFrEnergyDT::setDeltatAlchemical );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setDeltatAlchemical"
                , setDeltatAlchemical_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the delta alchemical used in the FEP method" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setDeviceIndex
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setDeviceIndex_function_type)( ::QString ) ;
            setDeviceIndex_function_type setDeviceIndex_function_value( &::SireMove::OpenMMFrEnergyDT::setDeviceIndex );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setDeviceIndex"
                , setDeviceIndex_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the OpenMM Platform: CUDA, OpenCL, CPU" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setEnergyFrequency
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setEnergyFrequency_function_type)( int ) ;
            setEnergyFrequency_function_type setEnergyFrequency_function_value( &::SireMove::OpenMMFrEnergyDT::setEnergyFrequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setEnergyFrequency"
                , setEnergyFrequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Center of Mass motion removal frequency" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setField_dielectric
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setField_dielectric_function_type)( double ) ;
            setField_dielectric_function_type setField_dielectric_function_value( &::SireMove::OpenMMFrEnergyDT::setField_dielectric );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setField_dielectric"
                , setField_dielectric_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the dielectric constant" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setMCBarostat
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setMCBarostat_function_type)( bool ) ;
            setMCBarostat_function_type setMCBarostat_function_value( &::SireMove::OpenMMFrEnergyDT::setMCBarostat );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setMCBarostat"
                , setMCBarostat_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set Monte Carlo Barostat onoff" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setMCBarostat_membrane

            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setMCBarostat_membrane_function_type)( bool ) ;
            setMCBarostat_membrane_function_type setMCBarostat_membrane_function_value( &::SireMove::OpenMMFrEnergyDT::setMCBarostat_membrane );

            OpenMMFrEnergyDT_exposer.def(
                "setMCBarostat_membrane"
                , setMCBarostat_membrane_function_value
                , ( bp::arg("arg0") )
                , "Set Monte Carlo Membrane Barostat onoff (only has an effect if Monte Carlo Barostat is on)" );

        }
        { //::SireMove::OpenMMFrEnergyDT::setMCBarostat_frequency
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setMCBarostat_frequency_function_type)( int ) ;
            setMCBarostat_frequency_function_type setMCBarostat_frequency_function_value( &::SireMove::OpenMMFrEnergyDT::setMCBarostat_frequency );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setMCBarostat_frequency"
                , setMCBarostat_frequency_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Monte Carlo Barostat frequency in time speps" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setPlatform
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setPlatform_function_type)( ::QString ) ;
            setPlatform_function_type setPlatform_function_value( &::SireMove::OpenMMFrEnergyDT::setPlatform );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setPlatform"
                , setPlatform_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the OpenMM Platform: CUDA, OpenCL, CPU" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setPressure
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setPressure_function_type)( ::SireUnits::Dimension::Pressure ) ;
            setPressure_function_type setPressure_function_value( &::SireMove::OpenMMFrEnergyDT::setPressure );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setPressure"
                , setPressure_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Pressure" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setRestraint
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setRestraint_function_type)( bool ) ;
            setRestraint_function_type setRestraint_function_value( &::SireMove::OpenMMFrEnergyDT::setRestraint );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setRestraint"
                , setRestraint_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Retraint mode" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setShift_delta
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setShift_delta_function_type)( double ) ;
            setShift_delta_function_type setShift_delta_function_value( &::SireMove::OpenMMFrEnergyDT::setShift_delta );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setShift_delta"
                , setShift_delta_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the shift used in the soft core potential" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::setTemperature
        
            typedef void ( ::SireMove::OpenMMFrEnergyDT::*setTemperature_function_type)( ::SireUnits::Dimension::Temperature ) ;
            setTemperature_function_type setTemperature_function_value( &::SireMove::OpenMMFrEnergyDT::setTemperature );
            
            OpenMMFrEnergyDT_exposer.def( 
                "setTemperature"
                , setTemperature_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Set the Temperature" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::toString
        
            typedef ::QString ( ::SireMove::OpenMMFrEnergyDT::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::OpenMMFrEnergyDT::toString );
            
            OpenMMFrEnergyDT_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this integrator" );
        
        }
        { //::SireMove::OpenMMFrEnergyDT::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::OpenMMFrEnergyDT::typeName );
            
            OpenMMFrEnergyDT_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        OpenMMFrEnergyDT_exposer.staticmethod( "typeName" );
        OpenMMFrEnergyDT_exposer.def( "__copy__", &__copy__);
        OpenMMFrEnergyDT_exposer.def( "__deepcopy__", &__copy__);
        OpenMMFrEnergyDT_exposer.def( "clone", &__copy__);
        OpenMMFrEnergyDT_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::OpenMMFrEnergyDT >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        OpenMMFrEnergyDT_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::OpenMMFrEnergyDT >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        OpenMMFrEnergyDT_exposer.def_pickle(sire_pickle_suite< ::SireMove::OpenMMFrEnergyDT >());
        OpenMMFrEnergyDT_exposer.def( "__str__", &__str__< ::SireMove::OpenMMFrEnergyDT > );
        OpenMMFrEnergyDT_exposer.def( "__repr__", &__str__< ::SireMove::OpenMMFrEnergyDT > );
    }

}
