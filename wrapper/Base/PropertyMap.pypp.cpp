// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "PropertyMap.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "propertymap.h"

#include <QDebug>

#include "propertymap.h"

SireBase::PropertyMap __copy__(const SireBase::PropertyMap &other){ return SireBase::PropertyMap(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_PropertyMap_class(){

    { //::SireBase::PropertyMap
        typedef bp::class_< SireBase::PropertyMap > PropertyMap_exposer_t;
        PropertyMap_exposer_t PropertyMap_exposer = PropertyMap_exposer_t( "PropertyMap", "This is the class that holds the collection of user-supplied\noptional properties and their locations to functions.\n\nThis class allows the following code to be written;\n\ncljff.add( mol, Property::set(charges,chgs) +\nProperty::set(ljs,ljparams) );\n\nThe PropertyMapPropertyName classes provide a kwargs\nlike interface for the C++ classes - indeed the python\nwrappers should allow code to be written like;\n\ncljff.add( mol, {charges : chgs, ljs : ljparams} )\n\nor\n\ncljff.add( mol, charges=charges, ljs=ljparams )\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope PropertyMap_scope( PropertyMap_exposer );
        PropertyMap_exposer.def( bp::init< QString const &, SireBase::PropertyName const & >(( bp::arg("property"), bp::arg("propname") ), "Construct a map that holds just a single PropertyName") );
        PropertyMap_exposer.def( bp::init< QHash< QString, SireBase::PropertyName > const & >(( bp::arg("propnames") ), "Construct a map that holds lots of PropertyNames") );
        PropertyMap_exposer.def( bp::init< SireBase::PropertyMap const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireBase::PropertyMap::addPrefix
        
            typedef ::SireBase::PropertyMap ( ::SireBase::PropertyMap::*addPrefix_function_type)( ::QString const &,::QStringList const & ) const;
            addPrefix_function_type addPrefix_function_value( &::SireBase::PropertyMap::addPrefix );
            
            PropertyMap_exposer.def( 
                "addPrefix"
                , addPrefix_function_value
                , ( bp::arg("prefix"), bp::arg("properties") )
                , bp::release_gil_policy()
                , "Add the passed prefix onto all of the passed properties,\n  returning a new property map that would map from\n  map[key] = X to map[prefixkey] = X\n" );
        
        }
        { //::SireBase::PropertyMap::addSuffix
        
            typedef ::SireBase::PropertyMap ( ::SireBase::PropertyMap::*addSuffix_function_type)( ::QString const &,::QStringList const & ) const;
            addSuffix_function_type addSuffix_function_value( &::SireBase::PropertyMap::addSuffix );
            
            PropertyMap_exposer.def( 
                "addSuffix"
                , addSuffix_function_value
                , ( bp::arg("suffix"), bp::arg("properties") )
                , bp::release_gil_policy()
                , "Add the passed suffix onto all of the passed properties,\n  returning a new property map that would map from\n  map[key] = X to map[keysuffix] = X\n" );
        
        }
        { //::SireBase::PropertyMap::isDefault
        
            typedef bool ( ::SireBase::PropertyMap::*isDefault_function_type)(  ) const;
            isDefault_function_type isDefault_function_value( &::SireBase::PropertyMap::isDefault );
            
            PropertyMap_exposer.def( 
                "isDefault"
                , isDefault_function_value
                , bp::release_gil_policy()
                , "Return whether or not this map is default - if it is,\nthen it doesnt specify any properties" );
        
        }
        { //::SireBase::PropertyMap::merge
        
            typedef ::SireBase::PropertyMap ( ::SireBase::PropertyMap::*merge_function_type)( ::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireBase::PropertyMap::merge );
            
            PropertyMap_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Return a PropertyMap that is the combination of this and other.\n  Keys set in other take precedence over keys in this.\n" );
        
        }
        PropertyMap_exposer.def( bp::self != bp::self );
        PropertyMap_exposer.def( bp::self + bp::self );
        { //::SireBase::PropertyMap::operator=
        
            typedef ::SireBase::PropertyMap & ( ::SireBase::PropertyMap::*assign_function_type)( ::SireBase::PropertyMap const & ) ;
            assign_function_type assign_function_value( &::SireBase::PropertyMap::operator= );
            
            PropertyMap_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PropertyMap_exposer.def( bp::self == bp::self );
        { //::SireBase::PropertyMap::operator[]
        
            typedef ::SireBase::PropertyName ( ::SireBase::PropertyMap::*__getitem___function_type)( char const * ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::PropertyMap::operator[] );
            
            PropertyMap_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireBase::PropertyMap::operator[]
        
            typedef ::SireBase::PropertyName ( ::SireBase::PropertyMap::*__getitem___function_type)( ::QString const & ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::PropertyMap::operator[] );
            
            PropertyMap_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireBase::PropertyMap::operator[]
        
            typedef ::SireBase::PropertyName ( ::SireBase::PropertyMap::*__getitem___function_type)( ::SireBase::PropertyName const & ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::PropertyMap::operator[] );
            
            PropertyMap_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireBase::PropertyMap::set
        
            typedef void ( ::SireBase::PropertyMap::*set_function_type)( ::QString const &,::SireBase::PropertyName const & ) ;
            set_function_type set_function_value( &::SireBase::PropertyMap::set );
            
            PropertyMap_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("name"), bp::arg("source") )
                , bp::release_gil_policy()
                , "Set the property called name to have the source or value\nin source. This replaces any existing source or value\nfor any existing property of this name in this map" );
        
        }
        { //::SireBase::PropertyMap::specified
        
            typedef bool ( ::SireBase::PropertyMap::*specified_function_type)( char const * ) const;
            specified_function_type specified_function_value( &::SireBase::PropertyMap::specified );
            
            PropertyMap_exposer.def( 
                "specified"
                , specified_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Return whether or not this map specifies the source or value\nof the property called name" );
        
        }
        { //::SireBase::PropertyMap::specified
        
            typedef bool ( ::SireBase::PropertyMap::*specified_function_type)( ::QString const & ) const;
            specified_function_type specified_function_value( &::SireBase::PropertyMap::specified );
            
            PropertyMap_exposer.def( 
                "specified"
                , specified_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Return whether or not this map specifies the source or value\nof the property called name" );
        
        }
        { //::SireBase::PropertyMap::specified
        
            typedef bool ( ::SireBase::PropertyMap::*specified_function_type)( ::SireBase::PropertyName const & ) const;
            specified_function_type specified_function_value( &::SireBase::PropertyMap::specified );
            
            PropertyMap_exposer.def( 
                "specified"
                , specified_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Return whether or not this map specifies the source or value\nof the property called name" );
        
        }
        { //::SireBase::PropertyMap::toDict
        
            typedef ::QHash< QString, SireBase::PropertyName > const ( ::SireBase::PropertyMap::*toDict_function_type)(  ) const;
            toDict_function_type toDict_function_value( &::SireBase::PropertyMap::toDict );
            
            PropertyMap_exposer.def( 
                "toDict"
                , toDict_function_value
                , bp::release_gil_policy()
                , "Return the raw underlying dictionary of the map" );
        
        }
        { //::SireBase::PropertyMap::toString
        
            typedef ::QString ( ::SireBase::PropertyMap::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::PropertyMap::toString );
            
            PropertyMap_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this PropertyMap" );
        
        }
        { //::SireBase::PropertyMap::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::PropertyMap::typeName );
            
            PropertyMap_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::PropertyMap::unset
        
            typedef void ( ::SireBase::PropertyMap::*unset_function_type)( ::QString const & ) ;
            unset_function_type unset_function_value( &::SireBase::PropertyMap::unset );
            
            PropertyMap_exposer.def( 
                "unset"
                , unset_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Unset the property called name. This will return it to default" );
        
        }
        { //::SireBase::PropertyMap::what
        
            typedef char const * ( ::SireBase::PropertyMap::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireBase::PropertyMap::what );
            
            PropertyMap_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        PropertyMap_exposer.staticmethod( "typeName" );
        PropertyMap_exposer.def( "__copy__", &__copy__<SireBase::PropertyMap>);
        PropertyMap_exposer.def( "__deepcopy__", &__copy__<SireBase::PropertyMap>);
        PropertyMap_exposer.def( "clone", &__copy__<SireBase::PropertyMap>);
        PropertyMap_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::PropertyMap >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PropertyMap_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::PropertyMap >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PropertyMap_exposer.def_pickle(sire_pickle_suite< ::SireBase::PropertyMap >());
        PropertyMap_exposer.def( "__str__", &__str__< ::SireBase::PropertyMap > );
        PropertyMap_exposer.def( "__repr__", &__str__< ::SireBase::PropertyMap > );
    }

}
