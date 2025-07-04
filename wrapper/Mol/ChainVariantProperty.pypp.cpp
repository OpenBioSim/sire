// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ChainVariantProperty.pypp.hpp"

namespace bp = boost::python;

#include "chainproperty.hpp"

#include "chainproperty.hpp"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

#include "SireMol/atomidxmapping.h"

SireMol::ChainProperty<QVariant> __copy__(const SireMol::ChainProperty<QVariant> &other){ return SireMol::ChainProperty<QVariant>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_ChainVariantProperty_class(){

    { //::SireMol::ChainProperty< QVariant >
        typedef bp::class_< SireMol::ChainProperty< QVariant >, bp::bases< SireMol::ChainProp, SireMol::MolViewProperty, SireBase::Property > > ChainVariantProperty_exposer_t;
        ChainVariantProperty_exposer_t ChainVariantProperty_exposer = ChainVariantProperty_exposer_t( "ChainVariantProperty", "", bp::init< >("") );
        bp::scope ChainVariantProperty_scope( ChainVariantProperty_exposer );
        ChainVariantProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        ChainVariantProperty_exposer.def( bp::init< QVector< QVariant > const & >(( bp::arg("values") ), "") );
        ChainVariantProperty_exposer.def( bp::init< SireMol::ChainProperty< QVariant > const & >(( bp::arg("other") ), "") );
        { //::SireMol::ChainProperty< QVariant >::array
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVector< QVariant > const & ( ::SireMol::ChainProperty< QVariant >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::ChainProperty< QVariant >::array );
            
            ChainVariantProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::assertCanConvert
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< QVariant >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ChainProperty< QVariant >::assertCanConvert );
            
            ChainVariantProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::assignFrom
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< QVariant >::*assignFrom_function_type)( ::SireMol::ChainProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ChainProperty< QVariant >::assignFrom );
            
            ChainVariantProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::at
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::ChainProperty< QVariant >::*at_function_type)( ::SireMol::ChainIdx const & ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< QVariant >::at );
            
            ChainVariantProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::at
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::ChainProperty< QVariant >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< QVariant >::at );
            
            ChainVariantProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::canConvert
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< QVariant >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ChainProperty< QVariant >::canConvert );
            
            ChainVariantProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::count
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< QVariant >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::ChainProperty< QVariant >::count );
            
            ChainVariantProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::fromVariant
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > ( *fromVariant_function_type )( ::SireMol::ChainProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::ChainProperty< QVariant >::fromVariant );
            
            ChainVariantProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::get
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::ChainProperty< QVariant >::*get_function_type)( ::SireMol::ChainIdx const & ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< QVariant >::get );
            
            ChainVariantProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::get
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::ChainProperty< QVariant >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< QVariant >::get );
            
            ChainVariantProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::getAsProperty
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::ChainProperty< QVariant >::*getAsProperty_function_type)( ::SireMol::ChainIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::ChainProperty< QVariant >::getAsProperty );
            
            ChainVariantProperty_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::getAsVariant
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVariant ( ::SireMol::ChainProperty< QVariant >::*getAsVariant_function_type)( ::SireMol::ChainIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::ChainProperty< QVariant >::getAsVariant );
            
            ChainVariantProperty_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::isCompatibleWith
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< QVariant >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::ChainProperty< QVariant >::isCompatibleWith );
            
            ChainVariantProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::isEmpty
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< QVariant >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::ChainProperty< QVariant >::isEmpty );
            
            ChainVariantProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::merge
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::SireBase::PropertyList ( ::SireMol::ChainProperty< QVariant >::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::ChainProperty< QVariant >::merge );
            
            ChainVariantProperty_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::nChains
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< QVariant >::*nChains_function_type)(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::ChainProperty< QVariant >::nChains );
            
            ChainVariantProperty_exposer.def( 
                "nChains"
                , nChains_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ChainVariantProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::ChainProperty< QVariant >::operator=
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > & ( ::SireMol::ChainProperty< QVariant >::*assign_function_type)( ::SireMol::ChainProperty< QVariant > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChainProperty< QVariant >::operator= );
            
            ChainVariantProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ChainVariantProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::ChainProperty< QVariant >::operator[]
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::ChainProperty< QVariant >::*__getitem___function_type)( ::SireMol::ChainIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< QVariant >::operator[] );
            
            ChainVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::operator[]
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QVariant const & ( ::SireMol::ChainProperty< QVariant >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< QVariant >::operator[] );
            
            ChainVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::operator[]
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QList< QVariant > ( ::SireMol::ChainProperty< QVariant >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< QVariant >::operator[] );
            
            ChainVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::operator[]
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QList< QVariant > ( ::SireMol::ChainProperty< QVariant >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< QVariant >::operator[] );
            
            ChainVariantProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::set
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > & ( ::SireMol::ChainProperty< QVariant >::*set_function_type)( ::SireMol::ChainIdx,::QVariant const & ) ;
            set_function_type set_function_value( &::SireMol::ChainProperty< QVariant >::set );
            
            ChainVariantProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("chainidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::size
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< QVariant >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::ChainProperty< QVariant >::size );
            
            ChainVariantProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::toString
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::QString ( ::SireMol::ChainProperty< QVariant >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChainProperty< QVariant >::toString );
            
            ChainVariantProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::toVariant
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > ( ::SireMol::ChainProperty< QVariant >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ChainProperty< QVariant >::toVariant );
            
            ChainVariantProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< QVariant >::typeName
        
            typedef SireMol::ChainProperty< QVariant > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChainProperty< QVariant >::typeName );
            
            ChainVariantProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ChainVariantProperty_exposer.staticmethod( "fromVariant" );
        ChainVariantProperty_exposer.staticmethod( "typeName" );
        ChainVariantProperty_exposer.def( "__copy__", &__copy__<SireMol::ChainProperty<QVariant>>);
        ChainVariantProperty_exposer.def( "__deepcopy__", &__copy__<SireMol::ChainProperty<QVariant>>);
        ChainVariantProperty_exposer.def( "clone", &__copy__<SireMol::ChainProperty<QVariant>>);
        ChainVariantProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChainProperty<QVariant> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainVariantProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChainProperty<QVariant> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainVariantProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::ChainProperty<QVariant> >());
        ChainVariantProperty_exposer.def( "__str__", &__str__< ::SireMol::ChainProperty<QVariant> > );
        ChainVariantProperty_exposer.def( "__repr__", &__str__< ::SireMol::ChainProperty<QVariant> > );
        ChainVariantProperty_exposer.def( "__len__", &__len_size< ::SireMol::ChainProperty<QVariant> > );
    }

}
