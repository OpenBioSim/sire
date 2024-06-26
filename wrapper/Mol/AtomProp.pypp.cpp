// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "AtomProp.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "atombeads.h"

#include "atomcharges.h"

#include "atomelements.h"

#include "atomenergies.h"

#include "atomforces.h"

#include "atommasses.h"

#include "atompolarisabilities.h"

#include "atomproperty.hpp"

#include "atompropertylist.h"

#include "atomradicals.h"

#include "atomradii.h"

#include "atomvelocities.h"

#include "atomproperty.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AtomProp_class(){

    { //::SireMol::AtomProp
        typedef bp::class_< SireMol::AtomProp, bp::bases< SireMol::MolViewProperty, SireBase::Property >, boost::noncopyable > AtomProp_exposer_t;
        AtomProp_exposer_t AtomProp_exposer = AtomProp_exposer_t( "AtomProp", "Small class used to give a common base to all\nAtomProperty classes", bp::no_init );
        bp::scope AtomProp_scope( AtomProp_exposer );
        { //::SireMol::AtomProp::assertCanConvert
        
            typedef void ( ::SireMol::AtomProp::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProp::assertCanConvert );
            
            AtomProp_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::assignFrom
        
            typedef void ( ::SireMol::AtomProp::*assignFrom_function_type)( ::SireMol::AtomVariantProperty const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProp::assignFrom );
            
            AtomProp_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::canConvert
        
            typedef bool ( ::SireMol::AtomProp::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProp::canConvert );
            
            AtomProp_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::divide
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProp::*divide_function_type)( ::QVector< SireMol::AtomSelection > const & ) const;
            divide_function_type divide_function_value( &::SireMol::AtomProp::divide );
            
            AtomProp_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("beads") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::divideByResidue
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProp::*divideByResidue_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            divideByResidue_function_type divideByResidue_function_value( &::SireMol::AtomProp::divideByResidue );
            
            AtomProp_exposer.def( 
                "divideByResidue"
                , divideByResidue_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::getAsProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProp::*getAsProperty_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::AtomProp::getAsProperty );
            
            AtomProp_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::getAsVariant
        
            typedef ::QVariant ( ::SireMol::AtomProp::*getAsVariant_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::AtomProp::getAsVariant );
            
            AtomProp_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::merge
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProp::*merge_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProp::merge );
            
            AtomProp_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProp::merge
        
            typedef ::SireBase::PropertyList ( ::SireMol::AtomProp::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProp::merge );
            
            AtomProp_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::AtomProp::operator=
        
            typedef ::SireMol::AtomProp & ( ::SireMol::AtomProp::*assign_function_type)( ::SireMol::AtomVariantProperty const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProp::operator= );
            
            AtomProp_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProp::toVariant
        
            typedef ::SireMol::AtomVariantProperty ( ::SireMol::AtomProp::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProp::toVariant );
            
            AtomProp_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomProp_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomProp_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomProp_exposer.def_pickle(sire_pickle_suite< ::SireMol::AtomProp >());
        AtomProp_exposer.def( "__str__", &__str__< ::SireMol::AtomProp > );
        AtomProp_exposer.def( "__repr__", &__str__< ::SireMol::AtomProp > );
    }

}
