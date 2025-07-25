// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Bead.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "bead.h"

#include "beadeditor.h"

#include "beads.h"

#include "chain.h"

#include "cutgroup.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "bead.h"

#include "beadproperty.hpp"

#include "SireMol/core.h"

const QString& get_Metadata_SireMol_BeadStringProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< QString >(metakey); }

const QString& get_Metadata_SireMol_BeadStringProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QString >(key, metakey); }

const qint64& get_Metadata_SireMol_BeadIntProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< qint64 >(metakey); }

const qint64& get_Metadata_SireMol_BeadIntProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< qint64 >(key, metakey); }

const double& get_Metadata_SireMol_BeadFloatProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< double >(metakey); }

const double& get_Metadata_SireMol_BeadFloatProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< double >(key, metakey); }

const QVariant& get_Metadata_SireMol_BeadVariantProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< QVariant >(metakey); }

const QVariant& get_Metadata_SireMol_BeadVariantProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QVariant >(key, metakey); }

const SireBase::PropertyPtr& get_Metadata_SireMol_BeadPropertyProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< SireBase::PropertyPtr >(metakey); }

const SireBase::PropertyPtr& get_Metadata_SireMol_BeadPropertyProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< SireBase::PropertyPtr >(key, metakey); }

SireMol::Bead __copy__(const SireMol::Bead &other){ return SireMol::Bead(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Bead_class(){

    { //::SireMol::Bead
        typedef bp::class_< SireMol::Bead, bp::bases< SireMol::MoleculeView, SireBase::Property > > Bead_exposer_t;
        Bead_exposer_t Bead_exposer = Bead_exposer_t( "Bead", "A Bead is a group of atoms (defined using a SireMol::Beading function)\nwithin a molecule. Beads can be used for coarse-graining, or for\nimplementing group-based cutoffs\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope Bead_scope( Bead_exposer );
        Bead_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::BeadIdx const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("moldata"), bp::arg("bead"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a bead view of the passed molecule, of the bead with index bead,\nusing the passed property map to find the required beading property") );
        Bead_exposer.def( bp::init< SireMol::Bead const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::Bead::assertContainsMetadata
        
            typedef void ( ::SireMol::Bead::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Bead::assertContainsMetadata );
            
            Bead_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that this bead contains the metadata with key metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Bead::assertContainsMetadata
        
            typedef void ( ::SireMol::Bead::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Bead::assertContainsMetadata );
            
            Bead_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that this bead contains the metadata property with key key:metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Bead::assertContainsProperty
        
            typedef void ( ::SireMol::Bead::*assertContainsProperty_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsProperty_function_type assertContainsProperty_function_value( &::SireMol::Bead::assertContainsProperty );
            
            Bead_exposer.def( 
                "assertContainsProperty"
                , assertContainsProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Assert that this bead contains the property with key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Bead::atom
        
            typedef ::SireMol::Atom ( ::SireMol::Bead::*atom_function_type)( int ) const;
            atom_function_type atom_function_value( &::SireMol::Bead::atom );
            
            Bead_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the ith atom in this bead\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::Bead::atomIdxs
        
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::Bead::*atomIdxs_function_type)(  ) const;
            atomIdxs_function_type atomIdxs_function_value( &::SireMol::Bead::atomIdxs );
            
            Bead_exposer.def( 
                "atomIdxs"
                , atomIdxs_function_value
                , bp::release_gil_policy()
                , "Return the list of atom indexes of the atoms in this bead" );
        
        }
        { //::SireMol::Bead::beading
        
            typedef ::SireMol::Beading const & ( ::SireMol::Bead::*beading_function_type)(  ) const;
            beading_function_type beading_function_value( &::SireMol::Bead::beading );
            
            Bead_exposer.def( 
                "beading"
                , beading_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the beading function used to bead up the molecule" );
        
        }
        { //::SireMol::Bead::beads
        
            typedef ::SireMol::Beads ( ::SireMol::Bead::*beads_function_type)(  ) const;
            beads_function_type beads_function_value( &::SireMol::Bead::beads );
            
            Bead_exposer.def( 
                "beads"
                , beads_function_value
                , bp::release_gil_policy()
                , "Return the set of all beads" );
        
        }
        { //::SireMol::Bead::contains
        
            typedef bool ( ::SireMol::Bead::*contains_function_type)( ::SireMol::AtomIdx ) const;
            contains_function_type contains_function_value( &::SireMol::Bead::contains );
            
            Bead_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomidx") )
                , bp::release_gil_policy()
                , "Return whether or not this bead contains the atom with index atomidx" );
        
        }
        { //::SireMol::Bead::contains
        
            typedef bool ( ::SireMol::Bead::*contains_function_type)( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Bead::contains );
            
            Bead_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not this bead contains the atom with ID atomid" );
        
        }
        { //::SireMol::Bead::edit
        
            typedef ::SireMol::BeadEditor ( ::SireMol::Bead::*edit_function_type)(  ) const;
            edit_function_type edit_function_value( &::SireMol::Bead::edit );
            
            Bead_exposer.def( 
                "edit"
                , edit_function_value
                , bp::release_gil_policy()
                , "Return the editor for this bead" );
        
        }
        { //::SireMol::Bead::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::Bead::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Bead::evaluate );
            
            Bead_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "Return the evaluator for this bead" );
        
        }
        { //::SireMol::Bead::hasMetadata
        
            typedef bool ( ::SireMol::Bead::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Bead::hasMetadata );
            
            Bead_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not this bead had some metadata at metakey" );
        
        }
        { //::SireMol::Bead::hasMetadata
        
            typedef bool ( ::SireMol::Bead::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Bead::hasMetadata );
            
            Bead_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not this bead has some metadata at key:metakey" );
        
        }
        { //::SireMol::Bead::hasProperty
        
            typedef bool ( ::SireMol::Bead::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Bead::hasProperty );
            
            Bead_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return whether or not this bead has a property called key" );
        
        }
        { //::SireMol::Bead::index
        
            typedef ::SireMol::BeadIdx ( ::SireMol::Bead::*index_function_type)(  ) const;
            index_function_type index_function_value( &::SireMol::Bead::index );
            
            Bead_exposer.def( 
                "index"
                , index_function_value
                , bp::release_gil_policy()
                , "Return the index of this bead" );
        
        }
        { //::SireMol::Bead::intersects
        
            typedef bool ( ::SireMol::Bead::*intersects_function_type)( ::SireMol::AtomID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Bead::intersects );
            
            Bead_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not this bead contains the atom with ID atomid" );
        
        }
        { //::SireMol::Bead::isEmpty
        
            typedef bool ( ::SireMol::Bead::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Bead::isEmpty );
            
            Bead_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return if this is an empty bead" );
        
        }
        { //::SireMol::Bead::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Bead::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Bead::metadataKeys );
            
            Bead_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "Return a list of all of the metadata properties associated with this bead" );
        
        }
        { //::SireMol::Bead::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Bead::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Bead::metadataKeys );
            
            Bead_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return a list of all of the metadata properties of the key key that\nare associated with this bead" );
        
        }
        { //::SireMol::Bead::move
        
            typedef ::SireMol::Mover< SireMol::Bead > ( ::SireMol::Bead::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Bead::move );
            
            Bead_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "Return the mover for this bead" );
        
        }
        { //::SireMol::Bead::nAtoms
        
            typedef int ( ::SireMol::Bead::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Bead::nAtoms );
            
            Bead_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in this bead" );
        
        }
        { //::SireMol::Bead::nViews
        
            typedef int ( ::SireMol::Bead::*nViews_function_type)(  ) const;
            nViews_function_type nViews_function_value( &::SireMol::Bead::nViews );
            
            Bead_exposer.def( 
                "nViews"
                , nViews_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Bead_exposer.def( bp::self != bp::self );
        { //::SireMol::Bead::operator=
        
            typedef ::SireMol::Bead & ( ::SireMol::Bead::*assign_function_type)( ::SireMol::Bead const & ) ;
            assign_function_type assign_function_value( &::SireMol::Bead::operator= );
            
            Bead_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Bead_exposer.def( bp::self == bp::self );
        { //::SireMol::Bead::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMol::Bead::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Bead::operator[] );
            
            Bead_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::Bead::propertyKeys
        
            typedef ::QStringList ( ::SireMol::Bead::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Bead::propertyKeys );
            
            Bead_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "Return a list of all of the properties associated with this bead" );
        
        }
        { //::SireMol::Bead::selectedAll
        
            typedef bool ( ::SireMol::Bead::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Bead::selectedAll );
            
            Bead_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "Return whether or not this bead includes all of the atoms in the molecule" );
        
        }
        { //::SireMol::Bead::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::Bead::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Bead::selection );
            
            Bead_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "Return the selection of atoms that are part of this bead" );
        
        }
        { //::SireMol::Bead::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMol::Bead::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMol::Bead::toSelector );
            
            Bead_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Bead::toString
        
            typedef ::QString ( ::SireMol::Bead::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Bead::toString );
            
            Bead_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this bead" );
        
        }
        { //::SireMol::Bead::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Bead::typeName );
            
            Bead_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Bead::update
        
            typedef void ( ::SireMol::Bead::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::Bead::update );
            
            Bead_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Update this bead using the passed molecule data" );
        
        }
        Bead_exposer.staticmethod( "typeName" );
        Bead_exposer.def( "_get_property_SireMol_BeadStringProperty", &SireMol::Bead::property< QString >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadStringProperty", get_Metadata_SireMol_BeadStringProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadStringProperty", &get_Metadata_SireMol_BeadStringProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_property_SireMol_BeadIntProperty", &SireMol::Bead::property< qint64 >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadIntProperty", get_Metadata_SireMol_BeadIntProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadIntProperty", &get_Metadata_SireMol_BeadIntProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_property_SireMol_BeadFloatProperty", &SireMol::Bead::property< double >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadFloatProperty", get_Metadata_SireMol_BeadFloatProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadFloatProperty", &get_Metadata_SireMol_BeadFloatProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_property_SireMol_BeadVariantProperty", &SireMol::Bead::property< QVariant >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadVariantProperty", get_Metadata_SireMol_BeadVariantProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadVariantProperty", &get_Metadata_SireMol_BeadVariantProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_property_SireMol_BeadPropertyProperty", &SireMol::Bead::property< SireBase::PropertyPtr >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadPropertyProperty", get_Metadata_SireMol_BeadPropertyProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadPropertyProperty", &get_Metadata_SireMol_BeadPropertyProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "__copy__", &__copy__<SireMol::Bead>);
        Bead_exposer.def( "__deepcopy__", &__copy__<SireMol::Bead>);
        Bead_exposer.def( "clone", &__copy__<SireMol::Bead>);
        Bead_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Bead >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Bead_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Bead >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Bead_exposer.def_pickle(sire_pickle_suite< ::SireMol::Bead >());
        Bead_exposer.def( "__str__", &__str__< ::SireMol::Bead > );
        Bead_exposer.def( "__repr__", &__str__< ::SireMol::Bead > );
        Bead_exposer.def( "__len__", &__len_size< ::SireMol::Bead > );
    }

}
