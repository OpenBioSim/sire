// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Dihedral.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireCAS/expression.h"

#include "SireCAS/symbol.h"

#include "SireCAS/values.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "dihedral.h"

#include "fouratomfunctions.h"

#include "selectordihedral.h"

#include <QDebug>

#include "dihedral.h"

SireMM::Dihedral __copy__(const SireMM::Dihedral &other){ return SireMM::Dihedral(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Dihedral_class(){

    { //::SireMM::Dihedral
        typedef bp::class_< SireMM::Dihedral, bp::bases< SireMol::MoleculeView, SireBase::Property > > Dihedral_exposer_t;
        Dihedral_exposer_t Dihedral_exposer = Dihedral_exposer_t( "Dihedral", "This class provides a molecule view to a dihedral", bp::init< >("") );
        bp::scope Dihedral_scope( Dihedral_exposer );
        Dihedral_exposer.def( bp::init< SireMol::Atom const &, SireMol::Atom const &, SireMol::Atom const &, SireMol::Atom const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") ), "") );
        Dihedral_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("molview"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") ), "") );
        Dihedral_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("moldata"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") ), "") );
        Dihedral_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::DihedralID const & >(( bp::arg("moldata"), bp::arg("dihedral") ), "") );
        Dihedral_exposer.def( bp::init< SireMM::Dihedral const & >(( bp::arg("other") ), "") );
        { //::SireMM::Dihedral::ID
        
            typedef ::SireMol::DihedralID ( ::SireMM::Dihedral::*ID_function_type)(  ) const;
            ID_function_type ID_function_value( &::SireMM::Dihedral::ID );
            
            Dihedral_exposer.def( 
                "ID"
                , ID_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::atom0
        
            typedef ::SireMol::Atom ( ::SireMM::Dihedral::*atom0_function_type)(  ) const;
            atom0_function_type atom0_function_value( &::SireMM::Dihedral::atom0 );
            
            Dihedral_exposer.def( 
                "atom0"
                , atom0_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::atom1
        
            typedef ::SireMol::Atom ( ::SireMM::Dihedral::*atom1_function_type)(  ) const;
            atom1_function_type atom1_function_value( &::SireMM::Dihedral::atom1 );
            
            Dihedral_exposer.def( 
                "atom1"
                , atom1_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::atom2
        
            typedef ::SireMol::Atom ( ::SireMM::Dihedral::*atom2_function_type)(  ) const;
            atom2_function_type atom2_function_value( &::SireMM::Dihedral::atom2 );
            
            Dihedral_exposer.def( 
                "atom2"
                , atom2_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::atom3
        
            typedef ::SireMol::Atom ( ::SireMM::Dihedral::*atom3_function_type)(  ) const;
            atom3_function_type atom3_function_value( &::SireMM::Dihedral::atom3 );
            
            Dihedral_exposer.def( 
                "atom3"
                , atom3_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::energy
        
            typedef ::SireUnits::Dimension::GeneralUnit ( ::SireMM::Dihedral::*energy_function_type)(  ) const;
            energy_function_type energy_function_value( &::SireMM::Dihedral::energy );
            
            Dihedral_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::energy
        
            typedef ::SireUnits::Dimension::GeneralUnit ( ::SireMM::Dihedral::*energy_function_type)( ::SireBase::PropertyMap const & ) const;
            energy_function_type energy_function_value( &::SireMM::Dihedral::energy );
            
            Dihedral_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMM::Dihedral::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMM::Dihedral::evaluate );
            
            Dihedral_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::hasMetadata
        
            typedef bool ( ::SireMM::Dihedral::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::Dihedral::hasMetadata );
            
            Dihedral_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::hasMetadata
        
            typedef bool ( ::SireMM::Dihedral::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::Dihedral::hasMetadata );
            
            Dihedral_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::hasProperty
        
            typedef bool ( ::SireMM::Dihedral::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMM::Dihedral::hasProperty );
            
            Dihedral_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::invert
        
            typedef ::SireMM::SelectorDihedral ( ::SireMM::Dihedral::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMM::Dihedral::invert );
            
            Dihedral_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::isEmpty
        
            typedef bool ( ::SireMM::Dihedral::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::Dihedral::isEmpty );
            
            Dihedral_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::measure
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Dihedral::*measure_function_type)(  ) const;
            measure_function_type measure_function_value( &::SireMM::Dihedral::measure );
            
            Dihedral_exposer.def( 
                "measure"
                , measure_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::measure
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Dihedral::*measure_function_type)( ::SireBase::PropertyMap const & ) const;
            measure_function_type measure_function_value( &::SireMM::Dihedral::measure );
            
            Dihedral_exposer.def( 
                "measure"
                , measure_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::metadataKeys
        
            typedef ::QStringList ( ::SireMM::Dihedral::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::Dihedral::metadataKeys );
            
            Dihedral_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::metadataKeys
        
            typedef ::QStringList ( ::SireMM::Dihedral::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::Dihedral::metadataKeys );
            
            Dihedral_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::move
        
            typedef ::SireMol::Mover< SireMM::Dihedral > ( ::SireMM::Dihedral::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMM::Dihedral::move );
            
            Dihedral_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Dihedral_exposer.def( bp::self != bp::self );
        { //::SireMM::Dihedral::operator=
        
            typedef ::SireMM::Dihedral & ( ::SireMM::Dihedral::*assign_function_type)( ::SireMM::Dihedral const & ) ;
            assign_function_type assign_function_value( &::SireMM::Dihedral::operator= );
            
            Dihedral_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Dihedral_exposer.def( bp::self == bp::self );
        { //::SireMM::Dihedral::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::Dihedral::*potential_function_type)(  ) const;
            potential_function_type potential_function_value( &::SireMM::Dihedral::potential );
            
            Dihedral_exposer.def( 
                "potential"
                , potential_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::potential
        
            typedef ::SireCAS::Expression ( ::SireMM::Dihedral::*potential_function_type)( ::SireBase::PropertyMap const & ) const;
            potential_function_type potential_function_value( &::SireMM::Dihedral::potential );
            
            Dihedral_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::properties
        
            typedef ::SireBase::Properties ( ::SireMM::Dihedral::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::Dihedral::properties );
            
            Dihedral_exposer.def( 
                "properties"
                , properties_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::property
        
            typedef ::SireBase::Property const & ( ::SireMM::Dihedral::*property_function_type)( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMM::Dihedral::property );
            
            Dihedral_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::Dihedral::property
        
            typedef ::SireBase::Property const & ( ::SireMM::Dihedral::*property_function_type)( ::SireBase::PropertyName const &,::SireBase::Property const & ) const;
            property_function_type property_function_value( &::SireMM::Dihedral::property );
            
            Dihedral_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key"), bp::arg("default_value") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMM::Dihedral::propertyKeys
        
            typedef ::QStringList ( ::SireMM::Dihedral::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMM::Dihedral::propertyKeys );
            
            Dihedral_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::selectedAll
        
            typedef bool ( ::SireMM::Dihedral::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMM::Dihedral::selectedAll );
            
            Dihedral_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMM::Dihedral::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMM::Dihedral::selection );
            
            Dihedral_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::selector
        
            typedef ::SireMM::SelectorDihedral ( ::SireMM::Dihedral::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMM::Dihedral::selector );
            
            Dihedral_exposer.def( 
                "selector"
                , selector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Dihedral::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::Dihedral::size );
            
            Dihedral_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMM::Dihedral::*size_function_type)( ::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMM::Dihedral::size );
            
            Dihedral_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::Dihedral::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMM::Dihedral::toSelector );
            
            Dihedral_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::toString
        
            typedef ::QString ( ::SireMM::Dihedral::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::Dihedral::toString );
            
            Dihedral_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::Dihedral::typeName );
            
            Dihedral_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::Dihedral::what
        
            typedef char const * ( ::SireMM::Dihedral::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::Dihedral::what );
            
            Dihedral_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Dihedral_exposer.staticmethod( "typeName" );
        Dihedral_exposer.def( "__copy__", &__copy__<SireMM::Dihedral>);
        Dihedral_exposer.def( "__deepcopy__", &__copy__<SireMM::Dihedral>);
        Dihedral_exposer.def( "clone", &__copy__<SireMM::Dihedral>);
        Dihedral_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::Dihedral >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Dihedral_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::Dihedral >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Dihedral_exposer.def_pickle(sire_pickle_suite< ::SireMM::Dihedral >());
        Dihedral_exposer.def( "__str__", &__str__< ::SireMM::Dihedral > );
        Dihedral_exposer.def( "__repr__", &__str__< ::SireMM::Dihedral > );
        Dihedral_exposer.def( "__len__", &__len_size< ::SireMM::Dihedral > );
    }

}
