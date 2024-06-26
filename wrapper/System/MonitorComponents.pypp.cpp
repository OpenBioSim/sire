// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MonitorComponents.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "monitorcomponent.h"

#include "monitorcomponents.h"

#include "system.h"

#include "tostring.h"

#include "monitorcomponents.h"

SireSystem::MonitorComponents __copy__(const SireSystem::MonitorComponents &other){ return SireSystem::MonitorComponents(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_MonitorComponents_class(){

    { //::SireSystem::MonitorComponents
        typedef bp::class_< SireSystem::MonitorComponents, bp::bases< SireSystem::SystemMonitor, SireBase::Property > > MonitorComponents_exposer_t;
        MonitorComponents_exposer_t MonitorComponents_exposer = MonitorComponents_exposer_t( "MonitorComponents", "This is a monitor that can be used to monitor large numbers\nof components of the system (of even all components). This\nmonitor is similar to MonitorComponent, but is better suited\nto situations where you want to monitor everything (or\neverything except for a small number of components), where\nit would be messy to specify lots of individual MonitorComponent\nmonitors\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope MonitorComponents_scope( MonitorComponents_exposer );
        MonitorComponents_exposer.def( bp::init< SireMaths::Accumulator const & >(( bp::arg("accumulator") ), "Construct a monitor to monitor all components using the\naccumulator accumulator") );
        MonitorComponents_exposer.def( bp::init< SireCAS::Symbol const &, bp::optional< SireMaths::Accumulator const & > >(( bp::arg("component"), bp::arg("accumulator")=SireMaths::Average() ), "Construct a monitor to monitor the components components\nusing the accumulator accumulator") );
        MonitorComponents_exposer.def( bp::init< SireCAS::Symbols const &, bp::optional< SireMaths::Accumulator const & > >(( bp::arg("components"), bp::arg("accumulator")=SireMaths::Average() ), "Construct a monitor to monitor the components components\nusing the accumulator accumulator") );
        MonitorComponents_exposer.def( bp::init< SireSystem::MonitorComponent const & >(( bp::arg("component_monitor") ), "Construct a MonitorComponents that is the (near) equivalent of\nthe MonitorComponent component_monitor") );
        MonitorComponents_exposer.def( bp::init< SireSystem::MonitorComponents const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::MonitorComponents::accumulator
        
            typedef ::SireMaths::Accumulator const & ( ::SireSystem::MonitorComponents::*accumulator_function_type)( ::SireCAS::Symbol const & ) const;
            accumulator_function_type accumulator_function_value( &::SireSystem::MonitorComponents::accumulator );
            
            MonitorComponents_exposer.def( 
                "accumulator"
                , accumulator_function_value
                , ( bp::arg("component") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the accumulator for the component component\nThrow: SireCAS::missing_symbol\n" );
        
        }
        { //::SireSystem::MonitorComponents::accumulatorTemplate
        
            typedef ::SireMaths::Accumulator const & ( ::SireSystem::MonitorComponents::*accumulatorTemplate_function_type)(  ) const;
            accumulatorTemplate_function_type accumulatorTemplate_function_value( &::SireSystem::MonitorComponents::accumulatorTemplate );
            
            MonitorComponents_exposer.def( 
                "accumulatorTemplate"
                , accumulatorTemplate_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the accumulator that is the template used for new accumulators\nthat are created when a new component is monitored" );
        
        }
        { //::SireSystem::MonitorComponents::clearStatistics
        
            typedef void ( ::SireSystem::MonitorComponents::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireSystem::MonitorComponents::clearStatistics );
            
            MonitorComponents_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "Completely clear the statistics" );
        
        }
        { //::SireSystem::MonitorComponents::excludeComponent
        
            typedef void ( ::SireSystem::MonitorComponents::*excludeComponent_function_type)( ::SireCAS::Symbol const & ) ;
            excludeComponent_function_type excludeComponent_function_value( &::SireSystem::MonitorComponents::excludeComponent );
            
            MonitorComponents_exposer.def( 
                "excludeComponent"
                , excludeComponent_function_value
                , ( bp::arg("component") )
                , bp::release_gil_policy()
                , "Make sure that the components in components are not monitored" );
        
        }
        { //::SireSystem::MonitorComponents::excludeComponent
        
            typedef void ( ::SireSystem::MonitorComponents::*excludeComponent_function_type)( ::SireCAS::Symbols const & ) ;
            excludeComponent_function_type excludeComponent_function_value( &::SireSystem::MonitorComponents::excludeComponent );
            
            MonitorComponents_exposer.def( 
                "excludeComponent"
                , excludeComponent_function_value
                , ( bp::arg("components") )
                , bp::release_gil_policy()
                , "Make sure that the components in components are not monitored" );
        
        }
        { //::SireSystem::MonitorComponents::excludeComponents
        
            typedef ::SireCAS::Symbols const & ( ::SireSystem::MonitorComponents::*excludeComponents_function_type)(  ) const;
            excludeComponents_function_type excludeComponents_function_value( &::SireSystem::MonitorComponents::excludeComponents );
            
            MonitorComponents_exposer.def( 
                "excludeComponents"
                , excludeComponents_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the components that will definitely not be monitored" );
        
        }
        { //::SireSystem::MonitorComponents::includeComponents
        
            typedef ::SireCAS::Symbols const & ( ::SireSystem::MonitorComponents::*includeComponents_function_type)(  ) const;
            includeComponents_function_type includeComponents_function_value( &::SireSystem::MonitorComponents::includeComponents );
            
            MonitorComponents_exposer.def( 
                "includeComponents"
                , includeComponents_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the components that will be monitored\n(if they exist, and not if they are excluded). If this\nis empty, then all components are monitored" );
        
        }
        { //::SireSystem::MonitorComponents::monitor
        
            typedef void ( ::SireSystem::MonitorComponents::*monitor_function_type)( ::SireSystem::System & ) ;
            monitor_function_type monitor_function_value( &::SireSystem::MonitorComponents::monitor );
            
            MonitorComponents_exposer.def( 
                "monitor"
                , monitor_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Monitor the system system - this will only accumulate symbols\nthat represent existing components - this does nothing for components\nthat dont exist in the system" );
        
        }
        { //::SireSystem::MonitorComponents::monitoredComponents
        
            typedef ::SireCAS::Symbols ( ::SireSystem::MonitorComponents::*monitoredComponents_function_type)(  ) const;
            monitoredComponents_function_type monitoredComponents_function_value( &::SireSystem::MonitorComponents::monitoredComponents );
            
            MonitorComponents_exposer.def( 
                "monitoredComponents"
                , monitoredComponents_function_value
                , bp::release_gil_policy()
                , "Return the set of symbols that have been monitored so far\n(so have valid accumulators)" );
        
        }
        MonitorComponents_exposer.def( bp::self != bp::self );
        { //::SireSystem::MonitorComponents::operator=
        
            typedef ::SireSystem::MonitorComponents & ( ::SireSystem::MonitorComponents::*assign_function_type)( ::SireSystem::MonitorComponents const & ) ;
            assign_function_type assign_function_value( &::SireSystem::MonitorComponents::operator= );
            
            MonitorComponents_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        MonitorComponents_exposer.def( bp::self == bp::self );
        { //::SireSystem::MonitorComponents::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::MonitorComponents::typeName );
            
            MonitorComponents_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MonitorComponents_exposer.staticmethod( "typeName" );
        MonitorComponents_exposer.def( "__copy__", &__copy__<SireSystem::MonitorComponents>);
        MonitorComponents_exposer.def( "__deepcopy__", &__copy__<SireSystem::MonitorComponents>);
        MonitorComponents_exposer.def( "clone", &__copy__<SireSystem::MonitorComponents>);
        MonitorComponents_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::MonitorComponents >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorComponents_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::MonitorComponents >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorComponents_exposer.def_pickle(sire_pickle_suite< ::SireSystem::MonitorComponents >());
        MonitorComponents_exposer.def( "__str__", &__str__< ::SireSystem::MonitorComponents > );
        MonitorComponents_exposer.def( "__repr__", &__str__< ::SireSystem::MonitorComponents > );
    }

}
