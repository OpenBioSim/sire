// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "SelectorBond.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/slice.h"

#include "SireCAS/symbol.h"

#include "SireCAS/values.h"

#include "SireID/index.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/selector.hpp"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "selectorbond.h"

#include "twoatomfunctions.h"

#include <QDebug>

#include "selectorbond.h"

SireMM::SelectorBond __copy__(const SireMM::SelectorBond &other){ return SireMM::SelectorBond(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_SelectorBond_class(){

    { //::SireMM::SelectorBond
        typedef bp::class_< SireMM::SelectorBond, bp::bases< SireMol::MoleculeView, SireBase::Property > > SelectorBond_exposer_t;
        SelectorBond_exposer_t SelectorBond_exposer = SelectorBond_exposer_t( "SelectorBond", "This provides a Selector<T>-style interface for multiple bonds", bp::init< >("") );
        bp::scope SelectorBond_scope( SelectorBond_exposer );
        SelectorBond_exposer.def( bp::init< SireMM::Bond const & >(( bp::arg("bond") ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeData const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::AtomID const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::BondID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("bond"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeView const &, QList< SireMol::BondID > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("bonds"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::AtomID const &, SireMol::AtomID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("atom0"), bp::arg("atom1"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atoms"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMol::Selector< SireMol::Atom > const &, SireMol::Selector< SireMol::Atom > const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atoms0"), bp::arg("atoms1"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        SelectorBond_exposer.def( bp::init< SireMM::SelectorBond const & >(( bp::arg("other") ), "") );
        { //::SireMM::SelectorBond::IDs
        
            typedef ::QList< SireMol::BondID > ( ::SireMM::SelectorBond::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireMM::SelectorBond::IDs );
            
            SelectorBond_exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::add
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*add_function_type)( ::SireMM::Bond const & ) const;
            add_function_type add_function_value( &::SireMM::SelectorBond::add );
            
            SelectorBond_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("bond") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::add
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*add_function_type)( ::SireMM::SelectorBond const & ) const;
            add_function_type add_function_value( &::SireMM::SelectorBond::add );
            
            SelectorBond_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::count
        
            typedef int ( ::SireMM::SelectorBond::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::SelectorBond::count );
            
            SelectorBond_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::energies
        
            typedef ::QList< SireUnits::Dimension::GeneralUnit > ( ::SireMM::SelectorBond::*energies_function_type)(  ) const;
            energies_function_type energies_function_value( &::SireMM::SelectorBond::energies );
            
            SelectorBond_exposer.def( 
                "energies"
                , energies_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::energies
        
            typedef ::QList< SireUnits::Dimension::GeneralUnit > ( ::SireMM::SelectorBond::*energies_function_type)( ::SireBase::PropertyMap const & ) const;
            energies_function_type energies_function_value( &::SireMM::SelectorBond::energies );
            
            SelectorBond_exposer.def( 
                "energies"
                , energies_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::energy
        
            typedef ::SireUnits::Dimension::GeneralUnit ( ::SireMM::SelectorBond::*energy_function_type)(  ) const;
            energy_function_type energy_function_value( &::SireMM::SelectorBond::energy );
            
            SelectorBond_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::energy
        
            typedef ::SireUnits::Dimension::GeneralUnit ( ::SireMM::SelectorBond::*energy_function_type)( ::SireBase::PropertyMap const & ) const;
            energy_function_type energy_function_value( &::SireMM::SelectorBond::energy );
            
            SelectorBond_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMM::SelectorBond::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMM::SelectorBond::evaluate );
            
            SelectorBond_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::hasMetadata
        
            typedef bool ( ::SireMM::SelectorBond::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::SelectorBond::hasMetadata );
            
            SelectorBond_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::hasMetadata
        
            typedef bool ( ::SireMM::SelectorBond::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMM::SelectorBond::hasMetadata );
            
            SelectorBond_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::hasProperty
        
            typedef bool ( ::SireMM::SelectorBond::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMM::SelectorBond::hasProperty );
            
            SelectorBond_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::intersection
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*intersection_function_type)( ::SireMM::SelectorBond const & ) const;
            intersection_function_type intersection_function_value( &::SireMM::SelectorBond::intersection );
            
            SelectorBond_exposer.def( 
                "intersection"
                , intersection_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::invert
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*invert_function_type)( ::SireBase::PropertyMap const & ) const;
            invert_function_type invert_function_value( &::SireMM::SelectorBond::invert );
            
            SelectorBond_exposer.def( 
                "invert"
                , invert_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::invert
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMM::SelectorBond::invert );
            
            SelectorBond_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::isEmpty
        
            typedef bool ( ::SireMM::SelectorBond::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::SelectorBond::isEmpty );
            
            SelectorBond_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::isSelector
        
            typedef bool ( ::SireMM::SelectorBond::*isSelector_function_type)(  ) const;
            isSelector_function_type isSelector_function_value( &::SireMM::SelectorBond::isSelector );
            
            SelectorBond_exposer.def( 
                "isSelector"
                , isSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::lengths
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 1, 0, 0, 0, 0, 0 > > ( ::SireMM::SelectorBond::*lengths_function_type)(  ) const;
            lengths_function_type lengths_function_value( &::SireMM::SelectorBond::lengths );
            
            SelectorBond_exposer.def( 
                "lengths"
                , lengths_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::lengths
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 1, 0, 0, 0, 0, 0 > > ( ::SireMM::SelectorBond::*lengths_function_type)( ::SireBase::PropertyMap const & ) const;
            lengths_function_type lengths_function_value( &::SireMM::SelectorBond::lengths );
            
            SelectorBond_exposer.def( 
                "lengths"
                , lengths_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::measures
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 1, 0, 0, 0, 0, 0 > > ( ::SireMM::SelectorBond::*measures_function_type)(  ) const;
            measures_function_type measures_function_value( &::SireMM::SelectorBond::measures );
            
            SelectorBond_exposer.def( 
                "measures"
                , measures_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::measures
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 1, 0, 0, 0, 0, 0 > > ( ::SireMM::SelectorBond::*measures_function_type)( ::SireBase::PropertyMap const & ) const;
            measures_function_type measures_function_value( &::SireMM::SelectorBond::measures );
            
            SelectorBond_exposer.def( 
                "measures"
                , measures_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::metadataKeys
        
            typedef ::QStringList ( ::SireMM::SelectorBond::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::SelectorBond::metadataKeys );
            
            SelectorBond_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::metadataKeys
        
            typedef ::QStringList ( ::SireMM::SelectorBond::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMM::SelectorBond::metadataKeys );
            
            SelectorBond_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::move
        
            typedef ::SireMol::Mover< SireMM::SelectorBond > ( ::SireMM::SelectorBond::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMM::SelectorBond::move );
            
            SelectorBond_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::nViews
        
            typedef int ( ::SireMM::SelectorBond::*nViews_function_type)(  ) const;
            nViews_function_type nViews_function_value( &::SireMM::SelectorBond::nViews );
            
            SelectorBond_exposer.def( 
                "nViews"
                , nViews_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SelectorBond_exposer.def( bp::self != bp::self );
        { //::SireMM::SelectorBond::operator()
        
            typedef ::SireMM::Bond ( ::SireMM::SelectorBond::*__call___function_type)( int ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorBond::operator() );
            
            SelectorBond_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator()
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorBond::operator() );
            
            SelectorBond_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator()
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*__call___function_type)( ::SireBase::Slice const & ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorBond::operator() );
            
            SelectorBond_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator()
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*__call___function_type)( ::QList< long long > const & ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorBond::operator() );
            
            SelectorBond_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator()
        
            typedef ::SireMM::SelectorBond ( ::SireMM::SelectorBond::*__call___function_type)( ::SireMol::BondID const & ) const;
            __call___function_type __call___function_value( &::SireMM::SelectorBond::operator() );
            
            SelectorBond_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("bond") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator=
        
            typedef ::SireMM::SelectorBond & ( ::SireMM::SelectorBond::*assign_function_type)( ::SireMM::SelectorBond const & ) ;
            assign_function_type assign_function_value( &::SireMM::SelectorBond::operator= );
            
            SelectorBond_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("bond") )
                , bp::return_self< >()
                , "" );
        
        }
        SelectorBond_exposer.def( bp::self == bp::self );
        { //::SireMM::SelectorBond::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorBond::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorBond::operator[] );
            
            SelectorBond_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorBond::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorBond::operator[] );
            
            SelectorBond_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorBond::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorBond::operator[] );
            
            SelectorBond_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::operator[]
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorBond::*__getitem___function_type)( ::SireMol::BondID const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::SelectorBond::operator[] );
            
            SelectorBond_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("bond") )
                , "" );
        
        }
        { //::SireMM::SelectorBond::potentials
        
            typedef ::QList< SireCAS::Expression > ( ::SireMM::SelectorBond::*potentials_function_type)(  ) const;
            potentials_function_type potentials_function_value( &::SireMM::SelectorBond::potentials );
            
            SelectorBond_exposer.def( 
                "potentials"
                , potentials_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::potentials
        
            typedef ::QList< SireCAS::Expression > ( ::SireMM::SelectorBond::*potentials_function_type)( ::SireBase::PropertyMap const & ) const;
            potentials_function_type potentials_function_value( &::SireMM::SelectorBond::potentials );
            
            SelectorBond_exposer.def( 
                "potentials"
                , potentials_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::properties
        
            typedef ::QList< SireBase::Properties > ( ::SireMM::SelectorBond::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::SelectorBond::properties );
            
            SelectorBond_exposer.def( 
                "properties"
                , properties_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::property
        
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMM::SelectorBond::*property_function_type)( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMM::SelectorBond::property );
            
            SelectorBond_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::property
        
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMM::SelectorBond::*property_function_type)( ::SireBase::PropertyName const &,::SireBase::Property const & ) const;
            property_function_type property_function_value( &::SireMM::SelectorBond::property );
            
            SelectorBond_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key"), bp::arg("default_value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::propertyKeys
        
            typedef ::QStringList ( ::SireMM::SelectorBond::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMM::SelectorBond::propertyKeys );
            
            SelectorBond_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::selectedAll
        
            typedef bool ( ::SireMM::SelectorBond::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMM::SelectorBond::selectedAll );
            
            SelectorBond_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMM::SelectorBond::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMM::SelectorBond::selection );
            
            SelectorBond_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::size
        
            typedef int ( ::SireMM::SelectorBond::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMM::SelectorBond::size );
            
            SelectorBond_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::toList
        
            typedef ::QList< SireBase::PropPtr< SireMol::MoleculeView > > ( ::SireMM::SelectorBond::*toList_function_type)(  ) const;
            toList_function_type toList_function_value( &::SireMM::SelectorBond::toList );
            
            SelectorBond_exposer.def( 
                "toList"
                , toList_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMM::SelectorBond::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMM::SelectorBond::toSelector );
            
            SelectorBond_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::toString
        
            typedef ::QString ( ::SireMM::SelectorBond::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::SelectorBond::toString );
            
            SelectorBond_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::SelectorBond::typeName );
            
            SelectorBond_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SelectorBond::what
        
            typedef char const * ( ::SireMM::SelectorBond::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::SelectorBond::what );
            
            SelectorBond_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SelectorBond_exposer.staticmethod( "typeName" );
        SelectorBond_exposer.def( "__copy__", &__copy__<SireMM::SelectorBond>);
        SelectorBond_exposer.def( "__deepcopy__", &__copy__<SireMM::SelectorBond>);
        SelectorBond_exposer.def( "clone", &__copy__<SireMM::SelectorBond>);
        SelectorBond_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::SelectorBond >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SelectorBond_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::SelectorBond >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SelectorBond_exposer.def_pickle(sire_pickle_suite< ::SireMM::SelectorBond >());
        SelectorBond_exposer.def( "__str__", &__str__< ::SireMM::SelectorBond > );
        SelectorBond_exposer.def( "__repr__", &__str__< ::SireMM::SelectorBond > );
        SelectorBond_exposer.def( "__len__", &__len_size< ::SireMM::SelectorBond > );
    }

}
