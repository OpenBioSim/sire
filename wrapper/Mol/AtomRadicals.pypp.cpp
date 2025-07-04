// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "AtomRadicals.pypp.hpp"

namespace bp = boost::python;

#include "atomradicals.h"

#include "atomradicals.h"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

#include "SireMol/atomidxmapping.h"

SireMol::AtomProperty<SireMol::Radical> __copy__(const SireMol::AtomProperty<SireMol::Radical> &other){ return SireMol::AtomProperty<SireMol::Radical>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_AtomRadicals_class(){

    { //::SireMol::AtomProperty< SireMol::Radical >
        typedef bp::class_< SireMol::AtomProperty< SireMol::Radical >, bp::bases< SireMol::AtomProp, SireMol::MolViewProperty, SireBase::Property > > AtomRadicals_exposer_t;
        AtomRadicals_exposer_t AtomRadicals_exposer = AtomRadicals_exposer_t( "AtomRadicals", "", bp::init< >("") );
        bp::scope AtomRadicals_scope( AtomRadicals_exposer );
        AtomRadicals_exposer.def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") ), "") );
        AtomRadicals_exposer.def( bp::init< SireMol::MoleculeInfo const &, SireMol::Radical const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomRadicals_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "") );
        AtomRadicals_exposer.def( bp::init< SireMol::MoleculeView const &, SireMol::Radical const & >(( bp::arg("molview"), bp::arg("default_value") ), "") );
        AtomRadicals_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        AtomRadicals_exposer.def( bp::init< SireMol::MoleculeInfoData const &, SireMol::Radical const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomRadicals_exposer.def( bp::init< SireMol::Radical const & >(( bp::arg("value") ), "") );
        AtomRadicals_exposer.def( bp::init< SireBase::PackedArray2D< SireMol::Radical > const & >(( bp::arg("values") ), "") );
        AtomRadicals_exposer.def( bp::init< SireMol::AtomProperty< SireMol::Radical > const & >(( bp::arg("other") ), "") );
        { //::SireMol::AtomProperty< SireMol::Radical >::array
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMol::Radical > const & ( ::SireMol::AtomProperty< SireMol::Radical >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::AtomProperty< SireMol::Radical >::array );
            
            AtomRadicals_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::assertCanConvert
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Radical >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProperty< SireMol::Radical >::assertCanConvert );
            
            AtomRadicals_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::assignFrom
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Radical >::*assignFrom_function_type)( ::SireMol::AtomProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProperty< SireMol::Radical >::assignFrom );
            
            AtomRadicals_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::at
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMol::Radical >::Array const & ( ::SireMol::AtomProperty< SireMol::Radical >::*at_function_type)( ::SireMol::CGIdx ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireMol::Radical >::at );
            
            AtomRadicals_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::at
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::Radical const & ( ::SireMol::AtomProperty< SireMol::Radical >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireMol::Radical >::at );
            
            AtomRadicals_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::at
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::Radical const & ( ::SireMol::AtomProperty< SireMol::Radical >::*at_function_type)( ::SireMol::CGAtomIdx const & ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireMol::Radical >::at );
            
            AtomRadicals_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::canConvert
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMol::Radical >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProperty< SireMol::Radical >::canConvert );
            
            AtomRadicals_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::copyFrom
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Radical >::*copyFrom_function_type)( ::QVector< SireMol::Radical > const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireMol::Radical >::copyFrom );
            
            AtomRadicals_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::copyFrom
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMol::Radical >::*copyFrom_function_type)( ::QVector< SireMol::Radical > const &,::SireMol::AtomSelection const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireMol::Radical >::copyFrom );
            
            AtomRadicals_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values"), bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::count
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Radical >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::AtomProperty< SireMol::Radical >::count );
            
            AtomRadicals_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::divide
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireMol::Radical >::*divide_function_type)( ::QVector< SireMol::AtomSelection > const & ) const;
            divide_function_type divide_function_value( &::SireMol::AtomProperty< SireMol::Radical >::divide );
            
            AtomRadicals_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("beads") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::divideByResidue
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireMol::Radical >::*divideByResidue_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            divideByResidue_function_type divideByResidue_function_value( &::SireMol::AtomProperty< SireMol::Radical >::divideByResidue );
            
            AtomRadicals_exposer.def( 
                "divideByResidue"
                , divideByResidue_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::fromVariant
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Radical > ( *fromVariant_function_type )( ::SireMol::AtomProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::AtomProperty< SireMol::Radical >::fromVariant );
            
            AtomRadicals_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::get
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMol::Radical >::Array const & ( ::SireMol::AtomProperty< SireMol::Radical >::*get_function_type)( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireMol::Radical >::get );
            
            AtomRadicals_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::get
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::Radical const & ( ::SireMol::AtomProperty< SireMol::Radical >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireMol::Radical >::get );
            
            AtomRadicals_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::get
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::Radical const & ( ::SireMol::AtomProperty< SireMol::Radical >::*get_function_type)( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireMol::Radical >::get );
            
            AtomRadicals_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::getAsProperty
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireMol::Radical >::*getAsProperty_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::AtomProperty< SireMol::Radical >::getAsProperty );
            
            AtomRadicals_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::getAsVariant
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QVariant ( ::SireMol::AtomProperty< SireMol::Radical >::*getAsVariant_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::AtomProperty< SireMol::Radical >::getAsVariant );
            
            AtomRadicals_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMol::Radical >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireMol::Radical >::isCompatibleWith );
            
            AtomRadicals_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMol::Radical >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfo const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireMol::Radical >::isCompatibleWith );
            
            AtomRadicals_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::isEmpty
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMol::Radical >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::AtomProperty< SireMol::Radical >::isEmpty );
            
            AtomRadicals_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::matchToSelection
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Radical > ( ::SireMol::AtomProperty< SireMol::Radical >::*matchToSelection_function_type)( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMol::AtomProperty< SireMol::Radical >::matchToSelection );
            
            AtomRadicals_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::merge
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireMol::Radical >::*merge_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireMol::Radical >::merge );
            
            AtomRadicals_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::merge
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PropertyList ( ::SireMol::AtomProperty< SireMol::Radical >::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireMol::Radical >::merge );
            
            AtomRadicals_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::nAtoms
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Radical >::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireMol::Radical >::nAtoms );
            
            AtomRadicals_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::nAtoms
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Radical >::*nAtoms_function_type)( ::SireMol::CGIdx ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireMol::Radical >::nAtoms );
            
            AtomRadicals_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("cgidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::nCutGroups
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Radical >::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::AtomProperty< SireMol::Radical >::nCutGroups );
            
            AtomRadicals_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomRadicals_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomProperty< SireMol::Radical >::operator=
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Radical > & ( ::SireMol::AtomProperty< SireMol::Radical >::*assign_function_type)( ::SireMol::AtomProperty< SireMol::Radical > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProperty< SireMol::Radical >::operator= );
            
            AtomRadicals_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomRadicals_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomProperty< SireMol::Radical >::operator[]
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMol::Radical >::Array const & ( ::SireMol::AtomProperty< SireMol::Radical >::*__getitem___function_type)( ::SireMol::CGIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMol::Radical >::operator[] );
            
            AtomRadicals_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::operator[]
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::Radical const & ( ::SireMol::AtomProperty< SireMol::Radical >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMol::Radical >::operator[] );
            
            AtomRadicals_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::operator[]
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::Radical const & ( ::SireMol::AtomProperty< SireMol::Radical >::*__getitem___function_type)( ::SireMol::CGAtomIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMol::Radical >::operator[] );
            
            AtomRadicals_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::operator[]
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QList< SireMol::Radical > ( ::SireMol::AtomProperty< SireMol::Radical >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMol::Radical >::operator[] );
            
            AtomRadicals_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::operator[]
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QList< SireMol::Radical > ( ::SireMol::AtomProperty< SireMol::Radical >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMol::Radical >::operator[] );
            
            AtomRadicals_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::set
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Radical > & ( ::SireMol::AtomProperty< SireMol::Radical >::*set_function_type)( int,::SireMol::Radical const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireMol::Radical >::set );
            
            AtomRadicals_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::set
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Radical > & ( ::SireMol::AtomProperty< SireMol::Radical >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireMol::Radical const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireMol::Radical >::set );
            
            AtomRadicals_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgatomidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::set
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMol::Radical > & ( ::SireMol::AtomProperty< SireMol::Radical >::*set_function_type)( ::SireMol::CGIdx,::QVector< SireMol::Radical > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireMol::Radical >::set );
            
            AtomRadicals_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("values") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::size
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMol::Radical >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::AtomProperty< SireMol::Radical >::size );
            
            AtomRadicals_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::toList
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QList< SireMol::Radical > ( ::SireMol::AtomProperty< SireMol::Radical >::*toList_function_type)(  ) const;
            toList_function_type toList_function_value( &::SireMol::AtomProperty< SireMol::Radical >::toList );
            
            AtomRadicals_exposer.def( 
                "toList"
                , toList_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::toList
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QList< SireMol::Radical > ( ::SireMol::AtomProperty< SireMol::Radical >::*toList_function_type)( ::SireMol::AtomSelection const & ) const;
            toList_function_type toList_function_value( &::SireMol::AtomProperty< SireMol::Radical >::toList );
            
            AtomRadicals_exposer.def( 
                "toList"
                , toList_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::toString
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QString ( ::SireMol::AtomProperty< SireMol::Radical >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomProperty< SireMol::Radical >::toString );
            
            AtomRadicals_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::toVariant
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::SireMol::AtomProperty< QVariant > ( ::SireMol::AtomProperty< SireMol::Radical >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProperty< SireMol::Radical >::toVariant );
            
            AtomRadicals_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::toVector
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QVector< SireMol::Radical > ( ::SireMol::AtomProperty< SireMol::Radical >::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireMol::Radical >::toVector );
            
            AtomRadicals_exposer.def( 
                "toVector"
                , toVector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::toVector
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef ::QVector< SireMol::Radical > ( ::SireMol::AtomProperty< SireMol::Radical >::*toVector_function_type)( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireMol::Radical >::toVector );
            
            AtomRadicals_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMol::Radical >::typeName
        
            typedef SireMol::AtomProperty< SireMol::Radical > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomProperty< SireMol::Radical >::typeName );
            
            AtomRadicals_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomRadicals_exposer.staticmethod( "fromVariant" );
        AtomRadicals_exposer.staticmethod( "typeName" );
        AtomRadicals_exposer.def( "__copy__", &__copy__<SireMol::AtomProperty<SireMol::Radical>>);
        AtomRadicals_exposer.def( "__deepcopy__", &__copy__<SireMol::AtomProperty<SireMol::Radical>>);
        AtomRadicals_exposer.def( "clone", &__copy__<SireMol::AtomProperty<SireMol::Radical>>);
        AtomRadicals_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProperty<SireMol::Radical> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomRadicals_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProperty<SireMol::Radical> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomRadicals_exposer.def_pickle(sire_pickle_suite< ::SireMol::AtomProperty<SireMol::Radical> >());
        AtomRadicals_exposer.def( "__str__", &__str__< ::SireMol::AtomProperty<SireMol::Radical> > );
        AtomRadicals_exposer.def( "__repr__", &__str__< ::SireMol::AtomProperty<SireMol::Radical> > );
        AtomRadicals_exposer.def( "__len__", &__len_size< ::SireMol::AtomProperty<SireMol::Radical> > );
    }

}
