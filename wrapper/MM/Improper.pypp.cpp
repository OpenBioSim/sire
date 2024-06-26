// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Improper.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireCAS/expression.h"

#include "SireCAS/symbol.h"

#include "SireCAS/values.h"

#include "SireMaths/torsion.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "fouratomfunctions.h"

#include "improper.h"

#include "selectorimproper.h"

#include <QDebug>

#include "improper.h"

SireMM::Improper __copy__(const SireMM::Improper &other){ return SireMM::Improper(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Improper_class(){

    { //::SireMM::Improper
        typedef bp::class_< SireMM::Improper, bp::bases< SireMol::MoleculeView, SireBase::Property > > Improper_exposer_t;
        Improper_exposer_t Improper_exposer = Improper_exposer_t( "Improper", "This class provides a molecule view to an improper", bp::init< >("") );
        bp::scope Improper_scope( Improper_exposer );
        Improper_exposer.def( bp::init< SireMol::Atom const &, SireMol::Atom const &, SireMol::Atom const &, SireMol::Atom const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") ), "") );
        Improper_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("molview"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") ), "") );
        Improper_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("moldata"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") ), "") );
        Improper_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::ImproperID const & >(( bp::arg("moldata"), bp::arg("improper") ), "") );
        Improper_exposer.def( bp::init< SireMM::Improper const & >(( bp::arg("other") ), "") );
        { //::SireMM::Improper::ID
        
            typedef ::SireMol::ImproperID ( ::SireMM::Improper::*ID_function_type)(  ) const;
            ID_function_type ID_function_value( &::SireMM::Improper::ID );
            
            Improper_exposer.def( 
                "ID"
                , ID_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::atom0
        
            typedef ::SireMol::Atom ( ::SireMM::Improper::*atom0_function_type)(  ) const;
            atom0_function_type atom0_function_value( &::SireMM::Improper::atom0 );
            
            Improper_exposer.def( 
                "atom0"
                , atom0_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::atom1
        
            typedef ::SireMol::Atom ( ::SireMM::Improper::*atom1_function_type)(  ) const;
            atom1_function_type atom1_function_value( &::SireMM::Improper::atom1 );
            
            Improper_exposer.def( 
                "atom1"
                , atom1_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::atom2
        
            typedef ::SireMol::Atom ( ::SireMM::Improper::*atom2_function_type)(  ) const;
            atom2_function_type atom2_function_value( &::SireMM::Improper::atom2 );
            
            Improper_exposer.def( 
                "atom2"
                , atom2_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::atom3
        
            typedef ::SireMol::Atom ( ::SireMM::Improper::*atom3_function_type)(  ) const;
            atom3_function_type atom3_function_value( &::SireMM::Improper::atom3 );
            
            Improper_exposer.def( 
                "atom3"
                , atom3_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::energy
        
            typedef ::SireUnits::Dimension::GeneralUnit ( ::SireMM::Improper::*energy_function_type)(  ) const;
            energy_function_type energy_function_value( &::SireMM::Improper::energy );
            
            Improper_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::energy
        
            typedef ::SireUnits::Dimension::GeneralUnit ( ::SireMM::Improper::*energy_function_type)( ::SireBase::PropertyMap const & ) const;
            energy_function_type energy_function_value( &::SireMM::Improper::energy );
            
            Improper_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMM::Improper::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMM::Improper::evaluate );
            
            Improper_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::hasMetadata
        
            typedef bool ( ::SireMM::Improper::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::Improper::hasMetadata );
            
            Improper_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::hasMetadata
        
            typedef bool ( ::SireMM::Improper::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::Improper::hasMetadata );
            
            Improper_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::hasProperty
        
            typedef bool ( ::SireMM::Improper::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMM::Improper::hasProperty );
            
            Improper_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::invert
        
            typedef ::SireMM::SelectorImproper ( ::SireMM::Improper::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMM::Improper::invert );
            
            Improper_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::isEmpty
        
            typedef bool ( ::SireMM::Improper::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::Improper::isEmpty );
            
            Improper_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::measure
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*measure_function_type)(  ) const;
            measure_function_type measure_function_value( &::SireMM::Improper::measure );
            
            Improper_exposer.def( 
                "measure"
                , measure_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::measure
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*measure_function_type)( ::SireBase::PropertyMap const & ) const;
            measure_function_type measure_function_value( &::SireMM::Improper::measure );
            
            Improper_exposer.def( 
                "measure"
                , measure_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::metadataKeys
        
            typedef ::QStringList ( ::SireMM::Improper::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::Improper::metadataKeys );
            
            Improper_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::metadataKeys
        
            typedef ::QStringList ( ::SireMM::Improper::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::Improper::metadataKeys );
            
            Improper_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::move
        
            typedef ::SireMol::Mover< SireMM::Improper > ( ::SireMM::Improper::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMM::Improper::move );
            
            Improper_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Improper_exposer.def( bp::self != bp::self );
        { //::SireMM::Improper::operator=
        
            typedef ::SireMM::Improper & ( ::SireMM::Improper::*assign_function_type)( ::SireMM::Improper const & ) ;
            assign_function_type assign_function_value( &::SireMM::Improper::operator= );
            
            Improper_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Improper_exposer.def( bp::self == bp::self );
        { //::SireMM::Improper::phi
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*phi_function_type)(  ) const;
            phi_function_type phi_function_value( &::SireMM::Improper::phi );
            
            Improper_exposer.def( 
                "phi"
                , phi_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::phi
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*phi_function_type)( ::SireBase::PropertyMap const & ) const;
            phi_function_type phi_function_value( &::SireMM::Improper::phi );
            
            Improper_exposer.def( 
                "phi"
                , phi_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::Improper::*potential_function_type)(  ) const;
            potential_function_type potential_function_value( &::SireMM::Improper::potential );
            
            Improper_exposer.def( 
                "potential"
                , potential_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::Improper::*potential_function_type)( ::SireBase::PropertyMap const & ) const;
            potential_function_type potential_function_value( &::SireMM::Improper::potential );
            
            Improper_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::properties
        
            typedef ::SireBase::Properties ( ::SireMM::Improper::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::Improper::properties );
            
            Improper_exposer.def( 
                "properties"
                , properties_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::property
        
            typedef ::SireBase::Property const & ( ::SireMM::Improper::*property_function_type)( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMM::Improper::property );
            
            Improper_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::Improper::property
        
            typedef ::SireBase::Property const & ( ::SireMM::Improper::*property_function_type)( ::SireBase::PropertyName const &,::SireBase::Property const & ) const;
            property_function_type property_function_value( &::SireMM::Improper::property );
            
            Improper_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key"), bp::arg("default_value") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::Improper::propertyKeys
        
            typedef ::QStringList ( ::SireMM::Improper::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMM::Improper::propertyKeys );
            
            Improper_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::selectedAll
        
            typedef bool ( ::SireMM::Improper::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMM::Improper::selectedAll );
            
            Improper_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMM::Improper::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMM::Improper::selection );
            
            Improper_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::selector
        
            typedef ::SireMM::SelectorImproper ( ::SireMM::Improper::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMM::Improper::selector );
            
            Improper_exposer.def( 
                "selector"
                , selector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::Improper::size );
            
            Improper_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*size_function_type)( ::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMM::Improper::size );
            
            Improper_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::theta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*theta_function_type)(  ) const;
            theta_function_type theta_function_value( &::SireMM::Improper::theta );
            
            Improper_exposer.def( 
                "theta"
                , theta_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::theta
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Improper::*theta_function_type)( ::SireBase::PropertyMap const & ) const;
            theta_function_type theta_function_value( &::SireMM::Improper::theta );
            
            Improper_exposer.def( 
                "theta"
                , theta_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::Improper::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMM::Improper::toSelector );
            
            Improper_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::toString
        
            typedef ::QString ( ::SireMM::Improper::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::Improper::toString );
            
            Improper_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::Improper::typeName );
            
            Improper_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Improper::what
        
            typedef char const * ( ::SireMM::Improper::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::Improper::what );
            
            Improper_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Improper_exposer.staticmethod( "typeName" );
        Improper_exposer.def( "__copy__", &__copy__<SireMM::Improper>);
        Improper_exposer.def( "__deepcopy__", &__copy__<SireMM::Improper>);
        Improper_exposer.def( "clone", &__copy__<SireMM::Improper>);
        Improper_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::Improper >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Improper_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::Improper >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Improper_exposer.def_pickle(sire_pickle_suite< ::SireMM::Improper >());
        Improper_exposer.def( "__str__", &__str__< ::SireMM::Improper > );
        Improper_exposer.def( "__repr__", &__str__< ::SireMM::Improper > );
        Improper_exposer.def( "__len__", &__len_size< ::SireMM::Improper > );
    }

}
