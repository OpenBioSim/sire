// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "CGVariantProperty.pypp.hpp"

namespace bp = boost::python;

#include "cgproperty.hpp"

#include "cgproperty.hpp"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

#include "SireMol/atomidxmapping.h"

SireMol::CGProperty<QVariant> __copy__(const SireMol::CGProperty<QVariant> &other){ return SireMol::CGProperty<QVariant>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_CGVariantProperty_class(){

    { //::SireMol::CGProperty< QVariant >
        typedef bp::class_< SireMol::CGProperty< QVariant >, bp::bases< SireMol::CGProp, SireMol::MolViewProperty, SireBase::Property > > CGVariantProperty_exposer_t;
        CGVariantProperty_exposer_t CGVariantProperty_exposer = CGVariantProperty_exposer_t( "CGVariantProperty", "", bp::init< >("") );
        bp::scope CGVariantProperty_scope( CGVariantProperty_exposer );
        CGVariantProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        CGVariantProperty_exposer.def( bp::init< QVector< QVariant > const & >(( bp::arg("values") ), "") );
        CGVariantProperty_exposer.def( bp::init< SireMol::CGProperty< QVariant > const & >(( bp::arg("other") ), "") );
        { //::SireMol::CGProperty< QVariant >::array
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVector< QVariant > const & ( ::SireMol::CGProperty< QVariant >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::CGProperty< QVariant >::array );
            
            CGVariantProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::assertCanConvert
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef void ( ::SireMol::CGProperty< QVariant >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::CGProperty< QVariant >::assertCanConvert );
            
            CGVariantProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::assignFrom
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef void ( ::SireMol::CGProperty< QVariant >::*assignFrom_function_type)( ::SireMol::CGProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::CGProperty< QVariant >::assignFrom );
            
            CGVariantProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("variant") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::at
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::CGProperty< QVariant >::*at_function_type)( ::SireMol::CGIdx const & ) const;
            at_function_type at_function_value( &::SireMol::CGProperty< QVariant >::at );
            
            CGVariantProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::at
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::CGProperty< QVariant >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::CGProperty< QVariant >::at );
            
            CGVariantProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::canConvert
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef bool ( ::SireMol::CGProperty< QVariant >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::CGProperty< QVariant >::canConvert );
            
            CGVariantProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::count
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef int ( ::SireMol::CGProperty< QVariant >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::CGProperty< QVariant >::count );
            
            CGVariantProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::fromVariant
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::SireMol::CGProperty< QVariant > ( *fromVariant_function_type )( ::SireMol::CGProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::CGProperty< QVariant >::fromVariant );
            
            CGVariantProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::get
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::CGProperty< QVariant >::*get_function_type)( ::SireMol::CGIdx const & ) const;
            get_function_type get_function_value( &::SireMol::CGProperty< QVariant >::get );
            
            CGVariantProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::get
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::CGProperty< QVariant >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::CGProperty< QVariant >::get );
            
            CGVariantProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::getAsProperty
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::CGProperty< QVariant >::*getAsProperty_function_type)( ::SireMol::CGIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::CGProperty< QVariant >::getAsProperty );
            
            CGVariantProperty_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::getAsVariant
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVariant ( ::SireMol::CGProperty< QVariant >::*getAsVariant_function_type)( ::SireMol::CGIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::CGProperty< QVariant >::getAsVariant );
            
            CGVariantProperty_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::isCompatibleWith
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef bool ( ::SireMol::CGProperty< QVariant >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::CGProperty< QVariant >::isCompatibleWith );
            
            CGVariantProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::isEmpty
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef bool ( ::SireMol::CGProperty< QVariant >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::CGProperty< QVariant >::isEmpty );
            
            CGVariantProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::merge
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::SireBase::PropertyList ( ::SireMol::CGProperty< QVariant >::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::CGProperty< QVariant >::merge );
            
            CGVariantProperty_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::nCutGroups
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef int ( ::SireMol::CGProperty< QVariant >::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::CGProperty< QVariant >::nCutGroups );
            
            CGVariantProperty_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CGVariantProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::CGProperty< QVariant >::operator=
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::SireMol::CGProperty< QVariant > & ( ::SireMol::CGProperty< QVariant >::*assign_function_type)( ::SireMol::CGProperty< QVariant > const & ) ;
            assign_function_type assign_function_value( &::SireMol::CGProperty< QVariant >::operator= );
            
            CGVariantProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CGVariantProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::CGProperty< QVariant >::operator[]
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::CGProperty< QVariant >::*__getitem___function_type)( ::SireMol::CGIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::CGProperty< QVariant >::operator[] );
            
            CGVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::operator[]
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::CGProperty< QVariant >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::CGProperty< QVariant >::operator[] );
            
            CGVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::operator[]
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QList< QVariant > ( ::SireMol::CGProperty< QVariant >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::CGProperty< QVariant >::operator[] );
            
            CGVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::operator[]
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QList< QVariant > ( ::SireMol::CGProperty< QVariant >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::CGProperty< QVariant >::operator[] );
            
            CGVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::set
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::SireMol::CGProperty< QVariant > & ( ::SireMol::CGProperty< QVariant >::*set_function_type)( ::SireMol::CGIdx,::QVariant const & ) ;
            set_function_type set_function_value( &::SireMol::CGProperty< QVariant >::set );
            
            CGVariantProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::size
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef int ( ::SireMol::CGProperty< QVariant >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::CGProperty< QVariant >::size );
            
            CGVariantProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::toString
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::QString ( ::SireMol::CGProperty< QVariant >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::CGProperty< QVariant >::toString );
            
            CGVariantProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::toVariant
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef ::SireMol::CGProperty< QVariant > ( ::SireMol::CGProperty< QVariant >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::CGProperty< QVariant >::toVariant );
            
            CGVariantProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGProperty< QVariant >::typeName
        
            typedef SireMol::CGProperty< QVariant > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CGProperty< QVariant >::typeName );
            
            CGVariantProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CGVariantProperty_exposer.staticmethod( "fromVariant" );
        CGVariantProperty_exposer.staticmethod( "typeName" );
        CGVariantProperty_exposer.def( "__copy__", &__copy__<SireMol::CGProperty<QVariant>>);
        CGVariantProperty_exposer.def( "__deepcopy__", &__copy__<SireMol::CGProperty<QVariant>>);
        CGVariantProperty_exposer.def( "clone", &__copy__<SireMol::CGProperty<QVariant>>);
        CGVariantProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CGProperty<QVariant> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGVariantProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CGProperty<QVariant> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGVariantProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::CGProperty<QVariant> >());
        CGVariantProperty_exposer.def( "__str__", &__str__< ::SireMol::CGProperty<QVariant> > );
        CGVariantProperty_exposer.def( "__repr__", &__str__< ::SireMol::CGProperty<QVariant> > );
        CGVariantProperty_exposer.def( "__len__", &__len_size< ::SireMol::CGProperty<QVariant> > );
    }

}
