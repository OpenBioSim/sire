// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "PropertyList.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "propertylist.h"

SireBase::PropertyList __copy__(const SireBase::PropertyList &other){ return SireBase::PropertyList(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_PropertyList_class(){

    { //::SireBase::PropertyList
        typedef bp::class_< SireBase::PropertyList, bp::bases< SireBase::Property > > PropertyList_exposer_t;
        PropertyList_exposer_t PropertyList_exposer = PropertyList_exposer_t( "PropertyList", bp::init< >() );
        bp::scope PropertyList_scope( PropertyList_exposer );
        PropertyList_exposer.def( bp::init< QList< SireBase::PropPtr< SireBase::Property > > const & >(( bp::arg("props") )) );
        PropertyList_exposer.def( bp::init< SireBase::PropertyList const & >(( bp::arg("other") )) );
        { //::SireBase::PropertyList::append
        
            typedef void ( ::SireBase::PropertyList::*append_function_type)( ::SireBase::Property const & ) ;
            append_function_type append_function_value( &::SireBase::PropertyList::append );
            
            PropertyList_exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("property") ) );
        
        }
        { //::SireBase::PropertyList::append
        
            typedef void ( ::SireBase::PropertyList::*append_function_type)( ::QList< SireBase::PropPtr< SireBase::Property > > const & ) ;
            append_function_type append_function_value( &::SireBase::PropertyList::append );
            
            PropertyList_exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("props") ) );
        
        }
        { //::SireBase::PropertyList::array
        
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireBase::PropertyList::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireBase::PropertyList::array );
            
            PropertyList_exposer.def( 
                "array"
                , array_function_value );
        
        }
        { //::SireBase::PropertyList::at
        
            typedef ::SireBase::Property const & ( ::SireBase::PropertyList::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireBase::PropertyList::at );
            
            PropertyList_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireBase::PropertyList::clear
        
            typedef void ( ::SireBase::PropertyList::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireBase::PropertyList::clear );
            
            PropertyList_exposer.def( 
                "clear"
                , clear_function_value );
        
        }
        { //::SireBase::PropertyList::count
        
            typedef int ( ::SireBase::PropertyList::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireBase::PropertyList::count );
            
            PropertyList_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireBase::PropertyList::empty
        
            typedef bool ( ::SireBase::PropertyList::*empty_function_type)(  ) const;
            empty_function_type empty_function_value( &::SireBase::PropertyList::empty );
            
            PropertyList_exposer.def( 
                "empty"
                , empty_function_value );
        
        }
        { //::SireBase::PropertyList::insert
        
            typedef void ( ::SireBase::PropertyList::*insert_function_type)( int,::SireBase::Property const & ) ;
            insert_function_type insert_function_value( &::SireBase::PropertyList::insert );
            
            PropertyList_exposer.def( 
                "insert"
                , insert_function_value
                , ( bp::arg("i"), bp::arg("value") ) );
        
        }
        { //::SireBase::PropertyList::isEmpty
        
            typedef bool ( ::SireBase::PropertyList::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireBase::PropertyList::isEmpty );
            
            PropertyList_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireBase::PropertyList::mid
        
            typedef ::SireBase::PropertyList ( ::SireBase::PropertyList::*mid_function_type)( int,int ) const;
            mid_function_type mid_function_value( &::SireBase::PropertyList::mid );
            
            PropertyList_exposer.def( 
                "mid"
                , mid_function_value
                , ( bp::arg("pos"), bp::arg("length")=(int)(-0x00000000000000001) ) );
        
        }
        { //::SireBase::PropertyList::move
        
            typedef void ( ::SireBase::PropertyList::*move_function_type)( int,int ) ;
            move_function_type move_function_value( &::SireBase::PropertyList::move );
            
            PropertyList_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("from"), bp::arg("to") ) );
        
        }
        PropertyList_exposer.def( bp::self != bp::self );
        PropertyList_exposer.def( bp::self + bp::self );
        { //::SireBase::PropertyList::operator=
        
            typedef ::SireBase::PropertyList & ( ::SireBase::PropertyList::*assign_function_type)( ::SireBase::PropertyList const & ) ;
            assign_function_type assign_function_value( &::SireBase::PropertyList::operator= );
            
            PropertyList_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        PropertyList_exposer.def( bp::self == bp::self );
        { //::SireBase::PropertyList::operator[]
        
            typedef ::SireBase::Property const & ( ::SireBase::PropertyList::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::PropertyList::operator[] );
            
            PropertyList_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireBase::PropertyList::pop_back
        
            typedef void ( ::SireBase::PropertyList::*pop_back_function_type)(  ) ;
            pop_back_function_type pop_back_function_value( &::SireBase::PropertyList::pop_back );
            
            PropertyList_exposer.def( 
                "pop_back"
                , pop_back_function_value );
        
        }
        { //::SireBase::PropertyList::pop_front
        
            typedef void ( ::SireBase::PropertyList::*pop_front_function_type)(  ) ;
            pop_front_function_type pop_front_function_value( &::SireBase::PropertyList::pop_front );
            
            PropertyList_exposer.def( 
                "pop_front"
                , pop_front_function_value );
        
        }
        { //::SireBase::PropertyList::prepend
        
            typedef void ( ::SireBase::PropertyList::*prepend_function_type)( ::SireBase::Property const & ) ;
            prepend_function_type prepend_function_value( &::SireBase::PropertyList::prepend );
            
            PropertyList_exposer.def( 
                "prepend"
                , prepend_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireBase::PropertyList::push_back
        
            typedef void ( ::SireBase::PropertyList::*push_back_function_type)( ::SireBase::Property const & ) ;
            push_back_function_type push_back_function_value( &::SireBase::PropertyList::push_back );
            
            PropertyList_exposer.def( 
                "push_back"
                , push_back_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireBase::PropertyList::push_front
        
            typedef void ( ::SireBase::PropertyList::*push_front_function_type)( ::SireBase::Property const & ) ;
            push_front_function_type push_front_function_value( &::SireBase::PropertyList::push_front );
            
            PropertyList_exposer.def( 
                "push_front"
                , push_front_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireBase::PropertyList::removeAt
        
            typedef void ( ::SireBase::PropertyList::*removeAt_function_type)( int ) ;
            removeAt_function_type removeAt_function_value( &::SireBase::PropertyList::removeAt );
            
            PropertyList_exposer.def( 
                "removeAt"
                , removeAt_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireBase::PropertyList::removeFirst
        
            typedef void ( ::SireBase::PropertyList::*removeFirst_function_type)(  ) ;
            removeFirst_function_type removeFirst_function_value( &::SireBase::PropertyList::removeFirst );
            
            PropertyList_exposer.def( 
                "removeFirst"
                , removeFirst_function_value );
        
        }
        { //::SireBase::PropertyList::removeLast
        
            typedef void ( ::SireBase::PropertyList::*removeLast_function_type)(  ) ;
            removeLast_function_type removeLast_function_value( &::SireBase::PropertyList::removeLast );
            
            PropertyList_exposer.def( 
                "removeLast"
                , removeLast_function_value );
        
        }
        { //::SireBase::PropertyList::replace
        
            typedef void ( ::SireBase::PropertyList::*replace_function_type)( int,::SireBase::Property const & ) ;
            replace_function_type replace_function_value( &::SireBase::PropertyList::replace );
            
            PropertyList_exposer.def( 
                "replace"
                , replace_function_value
                , ( bp::arg("i"), bp::arg("value") ) );
        
        }
        { //::SireBase::PropertyList::size
        
            typedef int ( ::SireBase::PropertyList::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireBase::PropertyList::size );
            
            PropertyList_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireBase::PropertyList::swap
        
            typedef void ( ::SireBase::PropertyList::*swap_function_type)( ::SireBase::PropertyList & ) ;
            swap_function_type swap_function_value( &::SireBase::PropertyList::swap );
            
            PropertyList_exposer.def( 
                "swap"
                , swap_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireBase::PropertyList::swap
        
            typedef void ( ::SireBase::PropertyList::*swap_function_type)( int,int ) ;
            swap_function_type swap_function_value( &::SireBase::PropertyList::swap );
            
            PropertyList_exposer.def( 
                "swap"
                , swap_function_value
                , ( bp::arg("i"), bp::arg("j") ) );
        
        }
        { //::SireBase::PropertyList::takeAt
        
            typedef ::SireBase::PropertyPtr ( ::SireBase::PropertyList::*takeAt_function_type)( int ) ;
            takeAt_function_type takeAt_function_value( &::SireBase::PropertyList::takeAt );
            
            PropertyList_exposer.def( 
                "takeAt"
                , takeAt_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireBase::PropertyList::takeFirst
        
            typedef ::SireBase::PropertyPtr ( ::SireBase::PropertyList::*takeFirst_function_type)(  ) ;
            takeFirst_function_type takeFirst_function_value( &::SireBase::PropertyList::takeFirst );
            
            PropertyList_exposer.def( 
                "takeFirst"
                , takeFirst_function_value );
        
        }
        { //::SireBase::PropertyList::takeLast
        
            typedef ::SireBase::PropertyPtr ( ::SireBase::PropertyList::*takeLast_function_type)(  ) ;
            takeLast_function_type takeLast_function_value( &::SireBase::PropertyList::takeLast );
            
            PropertyList_exposer.def( 
                "takeLast"
                , takeLast_function_value );
        
        }
        { //::SireBase::PropertyList::toList
        
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireBase::PropertyList::*toList_function_type)(  ) const;
            toList_function_type toList_function_value( &::SireBase::PropertyList::toList );
            
            PropertyList_exposer.def( 
                "toList"
                , toList_function_value );
        
        }
        { //::SireBase::PropertyList::toString
        
            typedef ::QString ( ::SireBase::PropertyList::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::PropertyList::toString );
            
            PropertyList_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireBase::PropertyList::toVector
        
            typedef ::QVector< SireBase::PropPtr< SireBase::Property > > ( ::SireBase::PropertyList::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireBase::PropertyList::toVector );
            
            PropertyList_exposer.def( 
                "toVector"
                , toVector_function_value );
        
        }
        { //::SireBase::PropertyList::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::PropertyList::typeName );
            
            PropertyList_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireBase::PropertyList::value
        
            typedef ::SireBase::PropertyPtr ( ::SireBase::PropertyList::*value_function_type)( int ) const;
            value_function_type value_function_value( &::SireBase::PropertyList::value );
            
            PropertyList_exposer.def( 
                "value"
                , value_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireBase::PropertyList::value
        
            typedef ::SireBase::PropertyPtr ( ::SireBase::PropertyList::*value_function_type)( int,::SireBase::Property const & ) const;
            value_function_type value_function_value( &::SireBase::PropertyList::value );
            
            PropertyList_exposer.def( 
                "value"
                , value_function_value
                , ( bp::arg("i"), bp::arg("default_value") ) );
        
        }
        PropertyList_exposer.staticmethod( "typeName" );
        PropertyList_exposer.def( "__copy__", &__copy__);
        PropertyList_exposer.def( "__deepcopy__", &__copy__);
        PropertyList_exposer.def( "clone", &__copy__);
        PropertyList_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::PropertyList >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PropertyList_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::PropertyList >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PropertyList_exposer.def( "__str__", &__str__< ::SireBase::PropertyList > );
        PropertyList_exposer.def( "__repr__", &__str__< ::SireBase::PropertyList > );
        PropertyList_exposer.def( "__len__", &__len_size< ::SireBase::PropertyList > );
    }

}
