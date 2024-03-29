// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "EnergyTrajectory.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "energytrajectory.h"

#include "energytrajectory.h"

SireMaths::EnergyTrajectory __copy__(const SireMaths::EnergyTrajectory &other){ return SireMaths::EnergyTrajectory(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_EnergyTrajectory_class(){

    { //::SireMaths::EnergyTrajectory
        typedef bp::class_< SireMaths::EnergyTrajectory, bp::bases< SireBase::Property > > EnergyTrajectory_exposer_t;
        EnergyTrajectory_exposer_t EnergyTrajectory_exposer = EnergyTrajectory_exposer_t( "EnergyTrajectory", "This class holds the trajectory of energies, organised by\ntimestep the energy was recorded and the types of energy\n(e.g. kinetic, potential, values at different lambda windows)\n", bp::init< >("") );
        bp::scope EnergyTrajectory_scope( EnergyTrajectory_exposer );
        EnergyTrajectory_exposer.def( bp::init< SireMaths::EnergyTrajectory const & >(( bp::arg("other") ), "") );
        { //::SireMaths::EnergyTrajectory::clearProperties
        
            typedef void ( ::SireMaths::EnergyTrajectory::*clearProperties_function_type)(  ) ;
            clearProperties_function_type clearProperties_function_value( &::SireMaths::EnergyTrajectory::clearProperties );
            
            EnergyTrajectory_exposer.def( 
                "clearProperties"
                , clearProperties_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::count
        
            typedef int ( ::SireMaths::EnergyTrajectory::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMaths::EnergyTrajectory::count );
            
            EnergyTrajectory_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of time values (number of rows)" );
        
        }
        { //::SireMaths::EnergyTrajectory::energies
        
            typedef ::QVector< double > ( ::SireMaths::EnergyTrajectory::*energies_function_type)( ::QString const & ) const;
            energies_function_type energies_function_value( &::SireMaths::EnergyTrajectory::energies );
            
            EnergyTrajectory_exposer.def( 
                "energies"
                , energies_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return all of the energy values for the passed key (the energy-key column).\n  This is in the same order as the times, and is in the default internal\n  unit\n" );
        
        }
        { //::SireMaths::EnergyTrajectory::energies
        
            typedef ::QVector< double > ( ::SireMaths::EnergyTrajectory::*energies_function_type)( ::QString const &,::SireUnits::Dimension::GeneralUnit const & ) const;
            energies_function_type energies_function_value( &::SireMaths::EnergyTrajectory::energies );
            
            EnergyTrajectory_exposer.def( 
                "energies"
                , energies_function_value
                , ( bp::arg("key"), bp::arg("energy_unit") )
                , bp::release_gil_policy()
                , "Return all of the energies fro the passed key converted to the\n  passed unit\n" );
        
        }
        { //::SireMaths::EnergyTrajectory::get
        
            typedef ::QHash< QString, double > ( ::SireMaths::EnergyTrajectory::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMaths::EnergyTrajectory::get );
            
            EnergyTrajectory_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the time and energy components at the ith row.\n  Values are returned in internal units\n" );
        
        }
        { //::SireMaths::EnergyTrajectory::get
        
            typedef ::QHash< QString, double > ( ::SireMaths::EnergyTrajectory::*get_function_type)( int,::SireUnits::Dimension::GeneralUnit const &,::SireUnits::Dimension::GeneralUnit const & ) const;
            get_function_type get_function_value( &::SireMaths::EnergyTrajectory::get );
            
            EnergyTrajectory_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i"), bp::arg("time_unit"), bp::arg("energy_unit") )
                , bp::release_gil_policy()
                , "Return the time and energy components at the ith row.\n  Values are returned in the specified units\n" );
        
        }
        { //::SireMaths::EnergyTrajectory::getLabels
        
            typedef ::QHash< QString, QString > ( ::SireMaths::EnergyTrajectory::*getLabels_function_type)( int ) const;
            getLabels_function_type getLabels_function_value( &::SireMaths::EnergyTrajectory::getLabels );
            
            EnergyTrajectory_exposer.def( 
                "getLabels"
                , getLabels_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::getLabelsAsNumbers
        
            typedef ::QHash< QString, double > ( ::SireMaths::EnergyTrajectory::*getLabelsAsNumbers_function_type)( int ) const;
            getLabelsAsNumbers_function_type getLabelsAsNumbers_function_value( &::SireMaths::EnergyTrajectory::getLabelsAsNumbers );
            
            EnergyTrajectory_exposer.def( 
                "getLabelsAsNumbers"
                , getLabelsAsNumbers_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::hasProperty
        
            typedef bool ( ::SireMaths::EnergyTrajectory::*hasProperty_function_type)( ::SireBase::PropertyName const & ) ;
            hasProperty_function_type hasProperty_function_value( &::SireMaths::EnergyTrajectory::hasProperty );
            
            EnergyTrajectory_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::isEmpty
        
            typedef bool ( ::SireMaths::EnergyTrajectory::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMaths::EnergyTrajectory::isEmpty );
            
            EnergyTrajectory_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMaths::EnergyTrajectory::isNull
        
            typedef bool ( ::SireMaths::EnergyTrajectory::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMaths::EnergyTrajectory::isNull );
            
            EnergyTrajectory_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::keys
        
            typedef ::QStringList ( ::SireMaths::EnergyTrajectory::*keys_function_type)(  ) const;
            keys_function_type keys_function_value( &::SireMaths::EnergyTrajectory::keys );
            
            EnergyTrajectory_exposer.def( 
                "keys"
                , keys_function_value
                , bp::release_gil_policy()
                , "Return all of the energy keys" );
        
        }
        { //::SireMaths::EnergyTrajectory::labelKeys
        
            typedef ::QStringList ( ::SireMaths::EnergyTrajectory::*labelKeys_function_type)(  ) const;
            labelKeys_function_type labelKeys_function_value( &::SireMaths::EnergyTrajectory::labelKeys );
            
            EnergyTrajectory_exposer.def( 
                "labelKeys"
                , labelKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::labels
        
            typedef ::QVector< QString > ( ::SireMaths::EnergyTrajectory::*labels_function_type)( ::QString const & ) const;
            labels_function_type labels_function_value( &::SireMaths::EnergyTrajectory::labels );
            
            EnergyTrajectory_exposer.def( 
                "labels"
                , labels_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::labelsAsNumbers
        
            typedef ::QVector< double > ( ::SireMaths::EnergyTrajectory::*labelsAsNumbers_function_type)( ::QString const & ) const;
            labelsAsNumbers_function_type labelsAsNumbers_function_value( &::SireMaths::EnergyTrajectory::labelsAsNumbers );
            
            EnergyTrajectory_exposer.def( 
                "labelsAsNumbers"
                , labelsAsNumbers_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        EnergyTrajectory_exposer.def( bp::self != bp::self );
        { //::SireMaths::EnergyTrajectory::operator=
        
            typedef ::SireMaths::EnergyTrajectory & ( ::SireMaths::EnergyTrajectory::*assign_function_type)( ::SireMaths::EnergyTrajectory const & ) ;
            assign_function_type assign_function_value( &::SireMaths::EnergyTrajectory::operator= );
            
            EnergyTrajectory_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        EnergyTrajectory_exposer.def( bp::self == bp::self );
        { //::SireMaths::EnergyTrajectory::operator[]
        
            typedef ::QHash< QString, double > ( ::SireMaths::EnergyTrajectory::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMaths::EnergyTrajectory::operator[] );
            
            EnergyTrajectory_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::properties
        
            typedef ::SireBase::Properties const & ( ::SireMaths::EnergyTrajectory::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMaths::EnergyTrajectory::properties );
            
            EnergyTrajectory_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::property
        
            typedef ::SireBase::Property const & ( ::SireMaths::EnergyTrajectory::*property_function_type)( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMaths::EnergyTrajectory::property );
            
            EnergyTrajectory_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::propertyKeys
        
            typedef ::QStringList ( ::SireMaths::EnergyTrajectory::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMaths::EnergyTrajectory::propertyKeys );
            
            EnergyTrajectory_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::removeProperty
        
            typedef void ( ::SireMaths::EnergyTrajectory::*removeProperty_function_type)( ::QString const & ) ;
            removeProperty_function_type removeProperty_function_value( &::SireMaths::EnergyTrajectory::removeProperty );
            
            EnergyTrajectory_exposer.def( 
                "removeProperty"
                , removeProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::set
        
            typedef void ( ::SireMaths::EnergyTrajectory::*set_function_type)( ::SireUnits::Dimension::GeneralUnit const &,::QHash< QString, SireUnits::Dimension::GeneralUnit > const & ) ;
            set_function_type set_function_value( &::SireMaths::EnergyTrajectory::set );
            
            EnergyTrajectory_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("time"), bp::arg("energies") )
                , bp::release_gil_policy()
                , "Set the energies at time time to the components contained\n  in energies\n" );
        
        }
        { //::SireMaths::EnergyTrajectory::set
        
            typedef void ( ::SireMaths::EnergyTrajectory::*set_function_type)( ::SireUnits::Dimension::GeneralUnit const &,::QHash< QString, SireUnits::Dimension::GeneralUnit > const &,::QHash< QString, QString > const & ) ;
            set_function_type set_function_value( &::SireMaths::EnergyTrajectory::set );
            
            EnergyTrajectory_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("time"), bp::arg("energies"), bp::arg("labels") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::setProperty
        
            typedef void ( ::SireMaths::EnergyTrajectory::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireMaths::EnergyTrajectory::setProperty );
            
            EnergyTrajectory_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("key"), bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::size
        
            typedef int ( ::SireMaths::EnergyTrajectory::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMaths::EnergyTrajectory::size );
            
            EnergyTrajectory_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of time values (number of rows)" );
        
        }
        { //::SireMaths::EnergyTrajectory::times
        
            typedef ::QVector< double > ( ::SireMaths::EnergyTrajectory::*times_function_type)(  ) const;
            times_function_type times_function_value( &::SireMaths::EnergyTrajectory::times );
            
            EnergyTrajectory_exposer.def( 
                "times"
                , times_function_value
                , bp::release_gil_policy()
                , "Return all of the time values (the time column). This is\n  sorted from earliest to latest time, and is in the default internal\n  unit\n" );
        
        }
        { //::SireMaths::EnergyTrajectory::times
        
            typedef ::QVector< double > ( ::SireMaths::EnergyTrajectory::*times_function_type)( ::SireUnits::Dimension::GeneralUnit const & ) const;
            times_function_type times_function_value( &::SireMaths::EnergyTrajectory::times );
            
            EnergyTrajectory_exposer.def( 
                "times"
                , times_function_value
                , ( bp::arg("time_unit") )
                , bp::release_gil_policy()
                , "Return all of the times converted to the passed unit" );
        
        }
        { //::SireMaths::EnergyTrajectory::toString
        
            typedef ::QString ( ::SireMaths::EnergyTrajectory::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::EnergyTrajectory::toString );
            
            EnergyTrajectory_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::EnergyTrajectory::typeName );
            
            EnergyTrajectory_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::EnergyTrajectory::what
        
            typedef char const * ( ::SireMaths::EnergyTrajectory::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::EnergyTrajectory::what );
            
            EnergyTrajectory_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        EnergyTrajectory_exposer.staticmethod( "typeName" );
        EnergyTrajectory_exposer.def( "__copy__", &__copy__);
        EnergyTrajectory_exposer.def( "__deepcopy__", &__copy__);
        EnergyTrajectory_exposer.def( "clone", &__copy__);
        EnergyTrajectory_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::EnergyTrajectory >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        EnergyTrajectory_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::EnergyTrajectory >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        EnergyTrajectory_exposer.def_pickle(sire_pickle_suite< ::SireMaths::EnergyTrajectory >());
        EnergyTrajectory_exposer.def( "__str__", &__str__< ::SireMaths::EnergyTrajectory > );
        EnergyTrajectory_exposer.def( "__repr__", &__str__< ::SireMaths::EnergyTrajectory > );
        EnergyTrajectory_exposer.def( "__len__", &__len_size< ::SireMaths::EnergyTrajectory > );
    }

}
