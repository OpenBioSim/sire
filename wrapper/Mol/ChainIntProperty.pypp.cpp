// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ChainIntProperty.pypp.hpp"

namespace bp = boost::python;

#include "chainproperty.hpp"

#include "chainproperty.hpp"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

#include "SireMol/atomidxmapping.h"

SireMol::ChainProperty<long long> __copy__(const SireMol::ChainProperty<long long> &other){ return SireMol::ChainProperty<long long>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_ChainIntProperty_class(){

    { //::SireMol::ChainProperty< long long >
        typedef bp::class_< SireMol::ChainProperty< long long >, bp::bases< SireMol::ChainProp, SireMol::MolViewProperty, SireBase::Property > > ChainIntProperty_exposer_t;
        ChainIntProperty_exposer_t ChainIntProperty_exposer = ChainIntProperty_exposer_t( "ChainIntProperty", "", bp::init< >("") );
        bp::scope ChainIntProperty_scope( ChainIntProperty_exposer );
        ChainIntProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        ChainIntProperty_exposer.def( bp::init< QVector< long long > const & >(( bp::arg("values") ), "") );
        ChainIntProperty_exposer.def( bp::init< SireMol::ChainProperty< long long > const & >(( bp::arg("other") ), "") );
        { //::SireMol::ChainProperty< long long >::array
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::QVector< long long > const & ( ::SireMol::ChainProperty< long long >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::ChainProperty< long long >::array );
            
            ChainIntProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::assertCanConvert
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< long long >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ChainProperty< long long >::assertCanConvert );
            
            ChainIntProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::assignFrom
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef void ( ::SireMol::ChainProperty< long long >::*assignFrom_function_type)( ::SireMol::ChainProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ChainProperty< long long >::assignFrom );
            
            ChainIntProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::at
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long >::*at_function_type)( ::SireMol::ChainIdx const & ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< long long >::at );
            
            ChainIntProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::at
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::ChainProperty< long long >::at );
            
            ChainIntProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::canConvert
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< long long >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ChainProperty< long long >::canConvert );
            
            ChainIntProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::count
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< long long >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::ChainProperty< long long >::count );
            
            ChainIntProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::fromVariant
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::SireMol::ChainProperty< long long > ( *fromVariant_function_type )( ::SireMol::ChainProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::ChainProperty< long long >::fromVariant );
            
            ChainIntProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::get
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long >::*get_function_type)( ::SireMol::ChainIdx const & ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< long long >::get );
            
            ChainIntProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::get
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::ChainProperty< long long >::get );
            
            ChainIntProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::getAsProperty
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::ChainProperty< long long >::*getAsProperty_function_type)( ::SireMol::ChainIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::ChainProperty< long long >::getAsProperty );
            
            ChainIntProperty_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::getAsVariant
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::QVariant ( ::SireMol::ChainProperty< long long >::*getAsVariant_function_type)( ::SireMol::ChainIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::ChainProperty< long long >::getAsVariant );
            
            ChainIntProperty_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::isCompatibleWith
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< long long >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::ChainProperty< long long >::isCompatibleWith );
            
            ChainIntProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::isEmpty
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef bool ( ::SireMol::ChainProperty< long long >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::ChainProperty< long long >::isEmpty );
            
            ChainIntProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::merge
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::SireBase::PropertyList ( ::SireMol::ChainProperty< long long >::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::ChainProperty< long long >::merge );
            
            ChainIntProperty_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::nChains
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< long long >::*nChains_function_type)(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::ChainProperty< long long >::nChains );
            
            ChainIntProperty_exposer.def( 
                "nChains"
                , nChains_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ChainIntProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::ChainProperty< long long >::operator=
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::SireMol::ChainProperty< long long > & ( ::SireMol::ChainProperty< long long >::*assign_function_type)( ::SireMol::ChainProperty< long long > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChainProperty< long long >::operator= );
            
            ChainIntProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ChainIntProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::ChainProperty< long long >::operator[]
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long >::*__getitem___function_type)( ::SireMol::ChainIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< long long >::operator[] );
            
            ChainIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("chainidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::operator[]
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef long long int const & ( ::SireMol::ChainProperty< long long >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< long long >::operator[] );
            
            ChainIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::operator[]
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::QList< long long > ( ::SireMol::ChainProperty< long long >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< long long >::operator[] );
            
            ChainIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::operator[]
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::QList< long long > ( ::SireMol::ChainProperty< long long >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ChainProperty< long long >::operator[] );
            
            ChainIntProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::set
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::SireMol::ChainProperty< long long > & ( ::SireMol::ChainProperty< long long >::*set_function_type)( ::SireMol::ChainIdx,long long int const & ) ;
            set_function_type set_function_value( &::SireMol::ChainProperty< long long >::set );
            
            ChainIntProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("chainidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::size
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef int ( ::SireMol::ChainProperty< long long >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::ChainProperty< long long >::size );
            
            ChainIntProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::toString
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::QString ( ::SireMol::ChainProperty< long long >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChainProperty< long long >::toString );
            
            ChainIntProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::toVariant
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef ::SireMol::ChainProperty< QVariant > ( ::SireMol::ChainProperty< long long >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ChainProperty< long long >::toVariant );
            
            ChainIntProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ChainProperty< long long >::typeName
        
            typedef SireMol::ChainProperty< long long > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChainProperty< long long >::typeName );
            
            ChainIntProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ChainIntProperty_exposer.staticmethod( "fromVariant" );
        ChainIntProperty_exposer.staticmethod( "typeName" );
        ChainIntProperty_exposer.def( "__copy__", &__copy__<SireMol::ChainProperty<long long>>);
        ChainIntProperty_exposer.def( "__deepcopy__", &__copy__<SireMol::ChainProperty<long long>>);
        ChainIntProperty_exposer.def( "clone", &__copy__<SireMol::ChainProperty<long long>>);
        ChainIntProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChainProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainIntProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChainProperty<long long> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainIntProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::ChainProperty<long long> >());
        ChainIntProperty_exposer.def( "__str__", &__str__< ::SireMol::ChainProperty<long long> > );
        ChainIntProperty_exposer.def( "__repr__", &__str__< ::SireMol::ChainProperty<long long> > );
        ChainIntProperty_exposer.def( "__len__", &__len_size< ::SireMol::ChainProperty<long long> > );
    }

}
