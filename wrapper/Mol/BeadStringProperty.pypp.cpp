// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "BeadStringProperty.pypp.hpp"

namespace bp = boost::python;

#include "beadproperty.hpp"

#include "beadproperty.hpp"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

#include "SireMol/atomidxmapping.h"

SireMol::BeadProperty<QString> __copy__(const SireMol::BeadProperty<QString> &other){ return SireMol::BeadProperty<QString>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_BeadStringProperty_class(){

    { //::SireMol::BeadProperty< QString >
        typedef bp::class_< SireMol::BeadProperty< QString >, bp::bases< SireMol::BeadProp, SireMol::MolViewProperty, SireBase::Property > > BeadStringProperty_exposer_t;
        BeadStringProperty_exposer_t BeadStringProperty_exposer = BeadStringProperty_exposer_t( "BeadStringProperty", "", bp::init< >("") );
        bp::scope BeadStringProperty_scope( BeadStringProperty_exposer );
        BeadStringProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const &, SireMol::Beading const & >(( bp::arg("molinfo"), bp::arg("beading") ), "") );
        BeadStringProperty_exposer.def( bp::init< QVector< QString > const &, SireMol::Beading const & >(( bp::arg("values"), bp::arg("beading") ), "") );
        BeadStringProperty_exposer.def( bp::init< SireMol::BeadProperty< QString > const & >(( bp::arg("other") ), "") );
        { //::SireMol::BeadProperty< QString >::array
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QVector< QString > const & ( ::SireMol::BeadProperty< QString >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::BeadProperty< QString >::array );
            
            BeadStringProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::assertCanConvert
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef void ( ::SireMol::BeadProperty< QString >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::BeadProperty< QString >::assertCanConvert );
            
            BeadStringProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::assignFrom
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef void ( ::SireMol::BeadProperty< QString >::*assignFrom_function_type)( ::SireMol::BeadProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::BeadProperty< QString >::assignFrom );
            
            BeadStringProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::at
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::BeadProperty< QString >::*at_function_type)( ::SireMol::BeadIdx const & ) const;
            at_function_type at_function_value( &::SireMol::BeadProperty< QString >::at );
            
            BeadStringProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("beadidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::at
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::BeadProperty< QString >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::BeadProperty< QString >::at );
            
            BeadStringProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::canConvert
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::BeadProperty< QString >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::BeadProperty< QString >::canConvert );
            
            BeadStringProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::count
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef int ( ::SireMol::BeadProperty< QString >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::BeadProperty< QString >::count );
            
            BeadStringProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::fromVariant
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::SireMol::BeadProperty< QString > ( *fromVariant_function_type )( ::SireMol::BeadProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::BeadProperty< QString >::fromVariant );
            
            BeadStringProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::get
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::BeadProperty< QString >::*get_function_type)( ::SireMol::BeadIdx const & ) const;
            get_function_type get_function_value( &::SireMol::BeadProperty< QString >::get );
            
            BeadStringProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("beadidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::get
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::BeadProperty< QString >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::BeadProperty< QString >::get );
            
            BeadStringProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::isCompatibleWith
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::BeadProperty< QString >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::BeadProperty< QString >::isCompatibleWith );
            
            BeadStringProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::isEmpty
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::BeadProperty< QString >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::BeadProperty< QString >::isEmpty );
            
            BeadStringProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::merge
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::SireBase::PropertyList ( ::SireMol::BeadProperty< QString >::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::BeadProperty< QString >::merge );
            
            BeadStringProperty_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::nBeads
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef int ( ::SireMol::BeadProperty< QString >::*nBeads_function_type)(  ) const;
            nBeads_function_type nBeads_function_value( &::SireMol::BeadProperty< QString >::nBeads );
            
            BeadStringProperty_exposer.def( 
                "nBeads"
                , nBeads_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BeadStringProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::BeadProperty< QString >::operator=
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::SireMol::BeadProperty< QString > & ( ::SireMol::BeadProperty< QString >::*assign_function_type)( ::SireMol::BeadProperty< QString > const & ) ;
            assign_function_type assign_function_value( &::SireMol::BeadProperty< QString >::operator= );
            
            BeadStringProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        BeadStringProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::BeadProperty< QString >::operator[]
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::BeadProperty< QString >::*__getitem___function_type)( ::SireMol::BeadIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::BeadProperty< QString >::operator[] );
            
            BeadStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("beadidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::operator[]
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::BeadProperty< QString >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::BeadProperty< QString >::operator[] );
            
            BeadStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::operator[]
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QList< QString > ( ::SireMol::BeadProperty< QString >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::BeadProperty< QString >::operator[] );
            
            BeadStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::operator[]
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QList< QString > ( ::SireMol::BeadProperty< QString >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::BeadProperty< QString >::operator[] );
            
            BeadStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::set
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::SireMol::BeadProperty< QString > & ( ::SireMol::BeadProperty< QString >::*set_function_type)( ::SireMol::BeadIdx,::QString const & ) ;
            set_function_type set_function_value( &::SireMol::BeadProperty< QString >::set );
            
            BeadStringProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("beadidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::size
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef int ( ::SireMol::BeadProperty< QString >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::BeadProperty< QString >::size );
            
            BeadStringProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::toString
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::QString ( ::SireMol::BeadProperty< QString >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::BeadProperty< QString >::toString );
            
            BeadStringProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::toVariant
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef ::SireMol::BeadProperty< QVariant > ( ::SireMol::BeadProperty< QString >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::BeadProperty< QString >::toVariant );
            
            BeadStringProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadProperty< QString >::typeName
        
            typedef SireMol::BeadProperty< QString > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::BeadProperty< QString >::typeName );
            
            BeadStringProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BeadStringProperty_exposer.staticmethod( "fromVariant" );
        BeadStringProperty_exposer.staticmethod( "typeName" );
        BeadStringProperty_exposer.def( "__copy__", &__copy__<SireMol::BeadProperty<QString>>);
        BeadStringProperty_exposer.def( "__deepcopy__", &__copy__<SireMol::BeadProperty<QString>>);
        BeadStringProperty_exposer.def( "clone", &__copy__<SireMol::BeadProperty<QString>>);
        BeadStringProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::BeadProperty<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BeadStringProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::BeadProperty<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BeadStringProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::BeadProperty<QString> >());
        BeadStringProperty_exposer.def( "__str__", &__str__< ::SireMol::BeadProperty<QString> > );
        BeadStringProperty_exposer.def( "__repr__", &__str__< ::SireMol::BeadProperty<QString> > );
        BeadStringProperty_exposer.def( "__len__", &__len_size< ::SireMol::BeadProperty<QString> > );
    }

}
