// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Beads.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "beads.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "beads.h"

SireMol::Beads __copy__(const SireMol::Beads &other){ return SireMol::Beads(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Beads_class(){

    { //::SireMol::Beads
        typedef bp::class_< SireMol::Beads, bp::bases< SireMol::MoleculeView, SireBase::Property > > Beads_exposer_t;
        Beads_exposer_t Beads_exposer = Beads_exposer_t( "Beads", "This class is a view of all of the beads (for a specific Beading)\nin a molecule\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope Beads_scope( Beads_exposer );
        Beads_exposer.def( bp::init< SireMol::MoleculeData const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("moldata"), bp::arg("map")=SireBase::PropertyMap() ), "Construct the beads view of the passed molecule, using\nthe passed (optional) property map to find the beading property") );
        Beads_exposer.def( bp::init< SireMol::Beads const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::Beads::at
        
            typedef ::SireMol::Bead ( ::SireMol::Beads::*at_function_type)( ::SireMol::BeadIdx ) const;
            at_function_type at_function_value( &::SireMol::Beads::at );
            
            Beads_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("beadidx") )
                , bp::release_gil_policy()
                , "Return the bead at index beadidx\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::Beads::atomIdxs
        
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::Beads::*atomIdxs_function_type)(  ) const;
            atomIdxs_function_type atomIdxs_function_value( &::SireMol::Beads::atomIdxs );
            
            Beads_exposer.def( 
                "atomIdxs"
                , atomIdxs_function_value
                , bp::release_gil_policy()
                , "Return the indicies of all of the atoms in all of the beads" );
        
        }
        { //::SireMol::Beads::atomProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::Beads::*atomProperty_function_type)( ::SireBase::PropertyName const & ) const;
            atomProperty_function_type atomProperty_function_value( &::SireMol::Beads::atomProperty );
            
            Beads_exposer.def( 
                "atomProperty"
                , atomProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the atom properties for all of the atoms in the beads, in\nBeadIdxIndex order" );
        
        }
        { //::SireMol::Beads::bead
        
            typedef ::SireMol::Bead ( ::SireMol::Beads::*bead_function_type)( ::SireMol::BeadIdx ) const;
            bead_function_type bead_function_value( &::SireMol::Beads::bead );
            
            Beads_exposer.def( 
                "bead"
                , bead_function_value
                , ( bp::arg("beadidx") )
                , bp::release_gil_policy()
                , "Return the bead at index beadidx\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::Beads::beading
        
            typedef ::SireMol::Beading const & ( ::SireMol::Beads::*beading_function_type)(  ) const;
            beading_function_type beading_function_value( &::SireMol::Beads::beading );
            
            Beads_exposer.def( 
                "beading"
                , beading_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the beading function used to create the beads" );
        
        }
        { //::SireMol::Beads::contains
        
            typedef bool ( ::SireMol::Beads::*contains_function_type)( ::SireMol::AtomIdx ) const;
            contains_function_type contains_function_value( &::SireMol::Beads::contains );
            
            Beads_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomidx") )
                , bp::release_gil_policy()
                , "Return whether any of the beads contains the atom with index atomidx" );
        
        }
        { //::SireMol::Beads::contains
        
            typedef bool ( ::SireMol::Beads::*contains_function_type)( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Beads::contains );
            
            Beads_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not any of the beads contains the atom with ID atomid" );
        
        }
        { //::SireMol::Beads::count
        
            typedef int ( ::SireMol::Beads::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::Beads::count );
            
            Beads_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of beads" );
        
        }
        { //::SireMol::Beads::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::Beads::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Beads::evaluate );
            
            Beads_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "Return an evaluator for all of the beads" );
        
        }
        { //::SireMol::Beads::hasMetadata
        
            typedef bool ( ::SireMol::Beads::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Beads::hasMetadata );
            
            Beads_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "At the moment, the Beads object has no properties or metadata" );
        
        }
        { //::SireMol::Beads::hasMetadata
        
            typedef bool ( ::SireMol::Beads::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Beads::hasMetadata );
            
            Beads_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "At the moment, the Beads object has no properties or metadata" );
        
        }
        { //::SireMol::Beads::hasProperty
        
            typedef bool ( ::SireMol::Beads::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Beads::hasProperty );
            
            Beads_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "At the moment, the Beads object has no properties or metadata" );
        
        }
        { //::SireMol::Beads::intersects
        
            typedef bool ( ::SireMol::Beads::*intersects_function_type)( ::SireMol::AtomID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Beads::intersects );
            
            Beads_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not any of the beads contains the atom with ID atomid" );
        
        }
        { //::SireMol::Beads::isEmpty
        
            typedef bool ( ::SireMol::Beads::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Beads::isEmpty );
            
            Beads_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMol::Beads::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Beads::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Beads::metadataKeys );
            
            Beads_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "At the moment, the Beads object has no properties or metadata" );
        
        }
        { //::SireMol::Beads::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Beads::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Beads::metadataKeys );
            
            Beads_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "At the moment, the Beads object has no properties or metadata" );
        
        }
        { //::SireMol::Beads::move
        
            typedef ::SireMol::Mover< SireMol::Beads > ( ::SireMol::Beads::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Beads::move );
            
            Beads_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "Return a mover that acts of all of the beads" );
        
        }
        { //::SireMol::Beads::nAtoms
        
            typedef int ( ::SireMol::Beads::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Beads::nAtoms );
            
            Beads_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in the beads" );
        
        }
        { //::SireMol::Beads::nBeads
        
            typedef int ( ::SireMol::Beads::*nBeads_function_type)(  ) const;
            nBeads_function_type nBeads_function_value( &::SireMol::Beads::nBeads );
            
            Beads_exposer.def( 
                "nBeads"
                , nBeads_function_value
                , bp::release_gil_policy()
                , "Return the number of beads" );
        
        }
        Beads_exposer.def( bp::self != bp::self );
        { //::SireMol::Beads::operator=
        
            typedef ::SireMol::Beads & ( ::SireMol::Beads::*assign_function_type)( ::SireMol::Beads const & ) ;
            assign_function_type assign_function_value( &::SireMol::Beads::operator= );
            
            Beads_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Beads_exposer.def( bp::self == bp::self );
        { //::SireMol::Beads::operator[]
        
            typedef ::SireMol::Bead ( ::SireMol::Beads::*__getitem___function_type)( ::SireMol::BeadIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Beads::operator[] );
            
            Beads_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("beadidx") )
                , "" );
        
        }
        { //::SireMol::Beads::propertyKeys
        
            typedef ::QStringList ( ::SireMol::Beads::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Beads::propertyKeys );
            
            Beads_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "At the moment, the Beads object has no properties or metadata" );
        
        }
        { //::SireMol::Beads::selectedAll
        
            typedef bool ( ::SireMol::Beads::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Beads::selectedAll );
            
            Beads_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "Return whether or not these beads contain all atoms" );
        
        }
        { //::SireMol::Beads::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::Beads::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Beads::selection );
            
            Beads_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "Return the selection of atoms that are part of the beads" );
        
        }
        { //::SireMol::Beads::size
        
            typedef int ( ::SireMol::Beads::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::Beads::size );
            
            Beads_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of beads" );
        
        }
        { //::SireMol::Beads::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMol::Beads::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMol::Beads::toSelector );
            
            Beads_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Beads::toString
        
            typedef ::QString ( ::SireMol::Beads::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Beads::toString );
            
            Beads_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of these beads" );
        
        }
        { //::SireMol::Beads::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Beads::typeName );
            
            Beads_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Beads::update
        
            typedef void ( ::SireMol::Beads::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::Beads::update );
            
            Beads_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Update these beads with the new molecule data" );
        
        }
        Beads_exposer.staticmethod( "typeName" );
        Beads_exposer.def( "__copy__", &__copy__<SireMol::Beads>);
        Beads_exposer.def( "__deepcopy__", &__copy__<SireMol::Beads>);
        Beads_exposer.def( "clone", &__copy__<SireMol::Beads>);
        Beads_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Beads >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Beads_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Beads >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Beads_exposer.def_pickle(sire_pickle_suite< ::SireMol::Beads >());
        Beads_exposer.def( "__str__", &__str__< ::SireMol::Beads > );
        Beads_exposer.def( "__repr__", &__str__< ::SireMol::Beads > );
        Beads_exposer.def( "__len__", &__len_size< ::SireMol::Beads > );
    }

}
