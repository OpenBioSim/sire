// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "CutGroup.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "cgatomidx.h"

#include "cgeditor.h"

#include "cutgroup.h"

#include "evaluator.h"

#include "groupatomids.h"

#include "molecule.h"

#include "mover.hpp"

#include "mover_metaid.h"

#include "selector.hpp"

#include "selectorm.hpp"

#include <QDebug>

#include "cutgroup.h"

#include "cgproperty.hpp"

#include "SireMol/core.h"

const QString& get_Metadata_SireMol_CGStringProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< QString >(metakey); }

const QString& get_Metadata_SireMol_CGStringProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QString >(key, metakey); }

const qint64& get_Metadata_SireMol_CGIntProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< qint64 >(metakey); }

const qint64& get_Metadata_SireMol_CGIntProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< qint64 >(key, metakey); }

const double& get_Metadata_SireMol_CGFloatProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< double >(metakey); }

const double& get_Metadata_SireMol_CGFloatProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< double >(key, metakey); }

const QVariant& get_Metadata_SireMol_CGVariantProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< QVariant >(metakey); }

const QVariant& get_Metadata_SireMol_CGVariantProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QVariant >(key, metakey); }

const SireBase::PropertyPtr& get_Metadata_SireMol_CGPropertyProperty_function1(const SireMol::CutGroup &atom,
                                   const QString &metakey){ return atom.metadata< SireBase::PropertyPtr >(metakey); }

const SireBase::PropertyPtr& get_Metadata_SireMol_CGPropertyProperty_function2(const SireMol::CutGroup &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< SireBase::PropertyPtr >(key, metakey); }

SireMol::CutGroup __copy__(const SireMol::CutGroup &other){ return SireMol::CutGroup(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_CutGroup_class(){

    { //::SireMol::CutGroup
        typedef bp::class_< SireMol::CutGroup, bp::bases< SireMol::MoleculeView, SireBase::Property > > CutGroup_exposer_t;
        CutGroup_exposer_t CutGroup_exposer = CutGroup_exposer_t( "CutGroup", "A CutGroup is a logical grouping of Atoms into a single group that is\nconsidered for intermolecular non-bonded cutting, and for periodic boundaries.\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope CutGroup_scope( CutGroup_exposer );
        CutGroup_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::CGID const & >(( bp::arg("moldata"), bp::arg("cgid") ), "Construct the CutGroup at ID cgid in the molecule whose data\nis in moldata\nThrow: SireMol::missing_CutGroup\nThrow: SireMol::duplicate_CutGroup\nThrow: SireError::invalid_index\n") );
        CutGroup_exposer.def( bp::init< SireMol::CutGroup const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::CutGroup::assertContainsMetadata
        
            typedef void ( ::SireMol::CutGroup::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::CutGroup::assertContainsMetadata );
            
            CutGroup_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that this CutGroup has an CGProperty piece of metadata\nat metakey metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::CutGroup::assertContainsMetadata
        
            typedef void ( ::SireMol::CutGroup::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::CutGroup::assertContainsMetadata );
            
            CutGroup_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that the property at key key has an CGProperty\npiece of metadata at metakey metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::CutGroup::assertContainsProperty
        
            typedef void ( ::SireMol::CutGroup::*assertContainsProperty_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsProperty_function_type assertContainsProperty_function_value( &::SireMol::CutGroup::assertContainsProperty );
            
            CutGroup_exposer.def( 
                "assertContainsProperty"
                , assertContainsProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Assert that this CutGroup has an CGProperty at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::CutGroup::atomIdxs
        
            typedef ::QList< SireMol::AtomIdx > const & ( ::SireMol::CutGroup::*atomIdxs_function_type)(  ) const;
            atomIdxs_function_type atomIdxs_function_value( &::SireMol::CutGroup::atomIdxs );
            
            CutGroup_exposer.def( 
                "atomIdxs"
                , atomIdxs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the indicies of the atoms in this CutGroup, in the\norder they appear in this CutGroup" );
        
        }
        { //::SireMol::CutGroup::contains
        
            typedef bool ( ::SireMol::CutGroup::*contains_function_type)( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMol::CutGroup::contains );
            
            CutGroup_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not this CutGroup contains all of\nthe atoms that match the ID atomid" );
        
        }
        { //::SireMol::CutGroup::contains
        
            typedef bool ( ::SireMol::CutGroup::*contains_function_type)( ::SireMol::AtomIdx ) const;
            contains_function_type contains_function_value( &::SireMol::CutGroup::contains );
            
            CutGroup_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomidx") )
                , bp::release_gil_policy()
                , "Return whether or not this CutGroup contains the atom\nat index atomidx in the molecule" );
        
        }
        { //::SireMol::CutGroup::edit
        
            typedef ::SireMol::CGEditor ( ::SireMol::CutGroup::*edit_function_type)(  ) const;
            edit_function_type edit_function_value( &::SireMol::CutGroup::edit );
            
            CutGroup_exposer.def( 
                "edit"
                , edit_function_value
                , bp::release_gil_policy()
                , "Return an editor that can edit this CutGroup" );
        
        }
        { //::SireMol::CutGroup::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::CutGroup::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::CutGroup::evaluate );
            
            CutGroup_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "Return an evaluator that can evaluate properties\nof this CutGroup" );
        
        }
        { //::SireMol::CutGroup::hasMetadata
        
            typedef bool ( ::SireMol::CutGroup::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::CutGroup::hasMetadata );
            
            CutGroup_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not there is a CGProperty at metakey metakey" );
        
        }
        { //::SireMol::CutGroup::hasMetadata
        
            typedef bool ( ::SireMol::CutGroup::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::CutGroup::hasMetadata );
            
            CutGroup_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether the metadata at metakey metakey for the property\nat key key is a CGProperty\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::CutGroup::hasProperty
        
            typedef bool ( ::SireMol::CutGroup::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::CutGroup::hasProperty );
            
            CutGroup_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return whether or not there is a CGProperty at key key" );
        
        }
        { //::SireMol::CutGroup::index
        
            typedef ::SireMol::CGIdx ( ::SireMol::CutGroup::*index_function_type)(  ) const;
            index_function_type index_function_value( &::SireMol::CutGroup::index );
            
            CutGroup_exposer.def( 
                "index"
                , index_function_value
                , bp::release_gil_policy()
                , "Return the index of this CutGroup in the molecule" );
        
        }
        { //::SireMol::CutGroup::intersects
        
            typedef bool ( ::SireMol::CutGroup::*intersects_function_type)( ::SireMol::AtomID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::CutGroup::intersects );
            
            CutGroup_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not this CutGroup contains some of\nthe atoms that match the ID atomid" );
        
        }
        { //::SireMol::CutGroup::invert
        
            typedef ::SireMol::Selector< SireMol::CutGroup > ( ::SireMol::CutGroup::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMol::CutGroup::invert );
            
            CutGroup_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "Return a selector that has everything except this view" );
        
        }
        { //::SireMol::CutGroup::isEmpty
        
            typedef bool ( ::SireMol::CutGroup::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::CutGroup::isEmpty );
            
            CutGroup_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Is this CutGroup empty?" );
        
        }
        { //::SireMol::CutGroup::metadataKeys
        
            typedef ::QStringList ( ::SireMol::CutGroup::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::CutGroup::metadataKeys );
            
            CutGroup_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "Return the metakeys of all CGProperty metadata" );
        
        }
        { //::SireMol::CutGroup::metadataKeys
        
            typedef ::QStringList ( ::SireMol::CutGroup::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::CutGroup::metadataKeys );
            
            CutGroup_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the metakeys of all CGProperty metadata for\nthe property at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::CutGroup::move
        
            typedef ::SireMol::Mover< SireMol::CutGroup > ( ::SireMol::CutGroup::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::CutGroup::move );
            
            CutGroup_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "Return an object that can move a copy of this CutGroup" );
        
        }
        { //::SireMol::CutGroup::nAtoms
        
            typedef int ( ::SireMol::CutGroup::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::CutGroup::nAtoms );
            
            CutGroup_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in this CutGroup" );
        
        }
        { //::SireMol::CutGroup::name
        
            typedef ::SireMol::CGName const & ( ::SireMol::CutGroup::*name_function_type)(  ) const;
            name_function_type name_function_value( &::SireMol::CutGroup::name );
            
            CutGroup_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the name of this CutGroup" );
        
        }
        { //::SireMol::CutGroup::number
        
            typedef ::SireMol::CGIdx ( ::SireMol::CutGroup::*number_function_type)(  ) const;
            number_function_type number_function_value( &::SireMol::CutGroup::number );
            
            CutGroup_exposer.def( 
                "number"
                , number_function_value
                , bp::release_gil_policy()
                , "Return the number of this CutGroup (same as the index)" );
        
        }
        CutGroup_exposer.def( bp::self != bp::self );
        CutGroup_exposer.def( bp::self + bp::other< SireMol::SelectorM< SireMol::CutGroup > >() );
        CutGroup_exposer.def( bp::self + bp::other< SireMol::Selector< SireMol::CutGroup > >() );
        CutGroup_exposer.def( bp::self + bp::self );
        CutGroup_exposer.def( bp::self - bp::other< SireMol::SelectorM< SireMol::CutGroup > >() );
        CutGroup_exposer.def( bp::self - bp::other< SireMol::Selector< SireMol::CutGroup > >() );
        CutGroup_exposer.def( bp::self - bp::self );
        { //::SireMol::CutGroup::operator=
        
            typedef ::SireMol::CutGroup & ( ::SireMol::CutGroup::*assign_function_type)( ::SireMol::CutGroup const & ) ;
            assign_function_type assign_function_value( &::SireMol::CutGroup::operator= );
            
            CutGroup_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CutGroup_exposer.def( bp::self == bp::self );
        { //::SireMol::CutGroup::propertyAsProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::CutGroup::*propertyAsProperty_function_type)( ::SireBase::PropertyName const & ) const;
            propertyAsProperty_function_type propertyAsProperty_function_value( &::SireMol::CutGroup::propertyAsProperty );
            
            CutGroup_exposer.def( 
                "propertyAsProperty"
                , propertyAsProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the specified property as a PropertyPtr" );
        
        }
        { //::SireMol::CutGroup::propertyAsVariant
        
            typedef ::QVariant ( ::SireMol::CutGroup::*propertyAsVariant_function_type)( ::SireBase::PropertyName const & ) const;
            propertyAsVariant_function_type propertyAsVariant_function_value( &::SireMol::CutGroup::propertyAsVariant );
            
            CutGroup_exposer.def( 
                "propertyAsVariant"
                , propertyAsVariant_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the specified property as a QVariant" );
        
        }
        { //::SireMol::CutGroup::propertyKeys
        
            typedef ::QStringList ( ::SireMol::CutGroup::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::CutGroup::propertyKeys );
            
            CutGroup_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "Return the keys of all CGProperty properties" );
        
        }
        { //::SireMol::CutGroup::selectedAll
        
            typedef bool ( ::SireMol::CutGroup::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::CutGroup::selectedAll );
            
            CutGroup_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "Is this CutGroup the whole molecule?" );
        
        }
        { //::SireMol::CutGroup::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::CutGroup::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::CutGroup::selection );
            
            CutGroup_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "Return the atoms that are in this CutGroup" );
        
        }
        { //::SireMol::CutGroup::selector
        
            typedef ::SireMol::Selector< SireMol::CutGroup > ( ::SireMol::CutGroup::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMol::CutGroup::selector );
            
            CutGroup_exposer.def( 
                "selector"
                , selector_function_value
                , bp::release_gil_policy()
                , "Return a selector that can change the selection of CutGroups" );
        
        }
        { //::SireMol::CutGroup::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMol::CutGroup::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMol::CutGroup::toSelector );
            
            CutGroup_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CutGroup::toString
        
            typedef ::QString ( ::SireMol::CutGroup::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::CutGroup::toString );
            
            CutGroup_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this CutGroup" );
        
        }
        { //::SireMol::CutGroup::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CutGroup::typeName );
            
            CutGroup_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CutGroup_exposer.staticmethod( "typeName" );
        CutGroup_exposer.def( "_get_property_SireMol_CGStringProperty", &SireMol::CutGroup::property< QString >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGStringProperty", get_Metadata_SireMol_CGStringProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGStringProperty", &get_Metadata_SireMol_CGStringProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_property_SireMol_CGIntProperty", &SireMol::CutGroup::property< qint64 >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGIntProperty", get_Metadata_SireMol_CGIntProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGIntProperty", &get_Metadata_SireMol_CGIntProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_property_SireMol_CGFloatProperty", &SireMol::CutGroup::property< double >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGFloatProperty", get_Metadata_SireMol_CGFloatProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGFloatProperty", &get_Metadata_SireMol_CGFloatProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_property_SireMol_CGVariantProperty", &SireMol::CutGroup::property< QVariant >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGVariantProperty", get_Metadata_SireMol_CGVariantProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGVariantProperty", &get_Metadata_SireMol_CGVariantProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_property_SireMol_CGPropertyProperty", &SireMol::CutGroup::property< SireBase::PropertyPtr >, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGPropertyProperty", get_Metadata_SireMol_CGPropertyProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "_get_metadata_SireMol_CGPropertyProperty", &get_Metadata_SireMol_CGPropertyProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        CutGroup_exposer.def( "__copy__", &__copy__<SireMol::CutGroup>);
        CutGroup_exposer.def( "__deepcopy__", &__copy__<SireMol::CutGroup>);
        CutGroup_exposer.def( "clone", &__copy__<SireMol::CutGroup>);
        CutGroup_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CutGroup >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CutGroup_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CutGroup >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CutGroup_exposer.def_pickle(sire_pickle_suite< ::SireMol::CutGroup >());
        CutGroup_exposer.def( "__str__", &__str__< ::SireMol::CutGroup > );
        CutGroup_exposer.def( "__repr__", &__str__< ::SireMol::CutGroup > );
        CutGroup_exposer.def( "__len__", &__len_size< ::SireMol::CutGroup > );
    }

}
