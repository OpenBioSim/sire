// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "AtomIntegerArrayProperty.pypp.hpp"

namespace bp = boost::python;

#include "atompropertylist.h"

#include "atompropertylist.h"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

#include "SireMol/atomidxmapping.h"

SireMol::AtomProperty<SireBase::IntegerArrayProperty> __copy__(const SireMol::AtomProperty<SireBase::IntegerArrayProperty> &other){ return SireMol::AtomProperty<SireBase::IntegerArrayProperty>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_AtomIntegerArrayProperty_class(){

    { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >
        typedef bp::class_< SireMol::AtomProperty< SireBase::IntegerArrayProperty >, bp::bases< SireMol::AtomProp, SireMol::MolViewProperty, SireBase::Property > > AtomIntegerArrayProperty_exposer_t;
        AtomIntegerArrayProperty_exposer_t AtomIntegerArrayProperty_exposer = AtomIntegerArrayProperty_exposer_t( "AtomIntegerArrayProperty", "", bp::init< >("") );
        bp::scope AtomIntegerArrayProperty_scope( AtomIntegerArrayProperty_exposer );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireMol::MoleculeInfo const &, SireBase::IntegerArrayProperty const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireMol::MoleculeView const &, SireBase::IntegerArrayProperty const & >(( bp::arg("molview"), bp::arg("default_value") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const &, SireBase::IntegerArrayProperty const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireBase::IntegerArrayProperty const & >(( bp::arg("value") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireBase::PackedArray2D< SireBase::IntegerArrayProperty > const & >(( bp::arg("values") ), "") );
        AtomIntegerArrayProperty_exposer.def( bp::init< SireMol::AtomProperty< SireBase::IntegerArrayProperty > const & >(( bp::arg("other") ), "") );
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::array
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::IntegerArrayProperty > const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::array );
            
            AtomIntegerArrayProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::assertCanConvert
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::assertCanConvert );
            
            AtomIntegerArrayProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::assignFrom
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*assignFrom_function_type)( ::SireMol::AtomProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::assignFrom );
            
            AtomIntegerArrayProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::at
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::IntegerArrayProperty >::Array const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*at_function_type)( ::SireMol::CGIdx ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::at );
            
            AtomIntegerArrayProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::at
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::IntegerArrayProperty const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::at );
            
            AtomIntegerArrayProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::at
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::IntegerArrayProperty const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*at_function_type)( ::SireMol::CGAtomIdx const & ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::at );
            
            AtomIntegerArrayProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::canConvert
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::canConvert );
            
            AtomIntegerArrayProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::copyFrom
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*copyFrom_function_type)( ::QVector< SireBase::IntegerArrayProperty > const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::copyFrom );
            
            AtomIntegerArrayProperty_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::copyFrom
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*copyFrom_function_type)( ::QVector< SireBase::IntegerArrayProperty > const &,::SireMol::AtomSelection const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::copyFrom );
            
            AtomIntegerArrayProperty_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values"), bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::count
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::count );
            
            AtomIntegerArrayProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::divide
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*divide_function_type)( ::QVector< SireMol::AtomSelection > const & ) const;
            divide_function_type divide_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::divide );
            
            AtomIntegerArrayProperty_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("beads") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::divideByResidue
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*divideByResidue_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            divideByResidue_function_type divideByResidue_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::divideByResidue );
            
            AtomIntegerArrayProperty_exposer.def( 
                "divideByResidue"
                , divideByResidue_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::fromVariant
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::IntegerArrayProperty > ( *fromVariant_function_type )( ::SireMol::AtomProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::fromVariant );
            
            AtomIntegerArrayProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::get
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::IntegerArrayProperty >::Array const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*get_function_type)( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::get );
            
            AtomIntegerArrayProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::get
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::IntegerArrayProperty const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::get );
            
            AtomIntegerArrayProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::get
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::IntegerArrayProperty const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*get_function_type)( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::get );
            
            AtomIntegerArrayProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::getAsProperty
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*getAsProperty_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::getAsProperty );
            
            AtomIntegerArrayProperty_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::getAsVariant
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QVariant ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*getAsVariant_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::getAsVariant );
            
            AtomIntegerArrayProperty_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::isCompatibleWith );
            
            AtomIntegerArrayProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfo const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::isCompatibleWith );
            
            AtomIntegerArrayProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::isEmpty
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::isEmpty );
            
            AtomIntegerArrayProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::matchToSelection
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::IntegerArrayProperty > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*matchToSelection_function_type)( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::matchToSelection );
            
            AtomIntegerArrayProperty_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::merge
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*merge_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::merge );
            
            AtomIntegerArrayProperty_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::merge
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PropertyList ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::merge );
            
            AtomIntegerArrayProperty_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::nAtoms
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::nAtoms );
            
            AtomIntegerArrayProperty_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::nAtoms
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*nAtoms_function_type)( ::SireMol::CGIdx ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::nAtoms );
            
            AtomIntegerArrayProperty_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("cgidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::nCutGroups
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::nCutGroups );
            
            AtomIntegerArrayProperty_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomIntegerArrayProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator=
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::IntegerArrayProperty > & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*assign_function_type)( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator= );
            
            AtomIntegerArrayProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomIntegerArrayProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[]
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::IntegerArrayProperty >::Array const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*__getitem___function_type)( ::SireMol::CGIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[] );
            
            AtomIntegerArrayProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[]
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::IntegerArrayProperty const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[] );
            
            AtomIntegerArrayProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[]
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireBase::IntegerArrayProperty const & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*__getitem___function_type)( ::SireMol::CGAtomIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[] );
            
            AtomIntegerArrayProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[]
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QList< SireBase::IntegerArrayProperty > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[] );
            
            AtomIntegerArrayProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[]
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QList< SireBase::IntegerArrayProperty > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::operator[] );
            
            AtomIntegerArrayProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::set
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::IntegerArrayProperty > & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*set_function_type)( int,::SireBase::IntegerArrayProperty const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::set );
            
            AtomIntegerArrayProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::set
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::IntegerArrayProperty > & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireBase::IntegerArrayProperty const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::set );
            
            AtomIntegerArrayProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgatomidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::set
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::IntegerArrayProperty > & ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*set_function_type)( ::SireMol::CGIdx,::QVector< SireBase::IntegerArrayProperty > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::set );
            
            AtomIntegerArrayProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("values") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::size
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::size );
            
            AtomIntegerArrayProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toList
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QList< SireBase::IntegerArrayProperty > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*toList_function_type)(  ) const;
            toList_function_type toList_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toList );
            
            AtomIntegerArrayProperty_exposer.def( 
                "toList"
                , toList_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toList
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QList< SireBase::IntegerArrayProperty > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*toList_function_type)( ::SireMol::AtomSelection const & ) const;
            toList_function_type toList_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toList );
            
            AtomIntegerArrayProperty_exposer.def( 
                "toList"
                , toList_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toString
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QString ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toString );
            
            AtomIntegerArrayProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toVariant
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::SireMol::AtomProperty< QVariant > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toVariant );
            
            AtomIntegerArrayProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toVector
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QVector< SireBase::IntegerArrayProperty > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toVector );
            
            AtomIntegerArrayProperty_exposer.def( 
                "toVector"
                , toVector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toVector
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef ::QVector< SireBase::IntegerArrayProperty > ( ::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::*toVector_function_type)( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::toVector );
            
            AtomIntegerArrayProperty_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::typeName
        
            typedef SireMol::AtomProperty< SireBase::IntegerArrayProperty > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomProperty< SireBase::IntegerArrayProperty >::typeName );
            
            AtomIntegerArrayProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomIntegerArrayProperty_exposer.staticmethod( "fromVariant" );
        AtomIntegerArrayProperty_exposer.staticmethod( "typeName" );
        AtomIntegerArrayProperty_exposer.def( "__copy__", &__copy__<SireMol::AtomProperty<SireBase::IntegerArrayProperty>>);
        AtomIntegerArrayProperty_exposer.def( "__deepcopy__", &__copy__<SireMol::AtomProperty<SireBase::IntegerArrayProperty>>);
        AtomIntegerArrayProperty_exposer.def( "clone", &__copy__<SireMol::AtomProperty<SireBase::IntegerArrayProperty>>);
        AtomIntegerArrayProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProperty<SireBase::IntegerArrayProperty> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIntegerArrayProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProperty<SireBase::IntegerArrayProperty> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIntegerArrayProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::AtomProperty<SireBase::IntegerArrayProperty> >());
        AtomIntegerArrayProperty_exposer.def( "__str__", &__str__< ::SireMol::AtomProperty<SireBase::IntegerArrayProperty> > );
        AtomIntegerArrayProperty_exposer.def( "__repr__", &__str__< ::SireMol::AtomProperty<SireBase::IntegerArrayProperty> > );
        AtomIntegerArrayProperty_exposer.def( "__len__", &__len_size< ::SireMol::AtomProperty<SireBase::IntegerArrayProperty> > );
    }

}
