// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ChunkedVector_double_.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "chunkedvector.hpp"

#include "chunkedvector.hpp"

SireBase::ChunkedVector<double, 100> __copy__(const SireBase::ChunkedVector<double, 100> &other){ return SireBase::ChunkedVector<double, 100>(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireBase::ChunkedVector<double, 100>&){ return "SireBase::ChunkedVector<double, 100>";}

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_ChunkedVector_double__class(){

    { //::SireBase::ChunkedVector< double, 100 >
        typedef bp::class_< SireBase::ChunkedVector< double, 100 > > ChunkedVector_double__exposer_t;
        ChunkedVector_double__exposer_t ChunkedVector_double__exposer = ChunkedVector_double__exposer_t( "ChunkedVector_double_", "", bp::init< >("") );
        bp::scope ChunkedVector_double__scope( ChunkedVector_double__exposer );
        { //::SireBase::ChunkedVector< double, 100 >::const_iterator
            typedef bp::class_< SireBase::ChunkedVector< double, 100 >::const_iterator > const_iterator_exposer_t;
            const_iterator_exposer_t const_iterator_exposer = const_iterator_exposer_t( "const_iterator", "", bp::init< >("") );
            bp::scope const_iterator_scope( const_iterator_exposer );
            const_iterator_exposer.def( bp::init< SireBase::ChunkedVector< double, 100 >::const_iterator const & >(( bp::arg("other") ), "") );
            const_iterator_exposer.def( bp::self != bp::self );
            { //::SireBase::ChunkedVector< double, 100 >::const_iterator::operator=
            
                typedef ::SireBase::ChunkedVector< double, 100 >::const_iterator & ( ::SireBase::ChunkedVector< double, 100 >::const_iterator::*assign_function_type)( ::SireBase::ChunkedVector< double, 100 >::const_iterator const & ) ;
                assign_function_type assign_function_value( &::SireBase::ChunkedVector< double, 100 >::const_iterator::operator= );
                
                const_iterator_exposer.def( 
                    "assign"
                    , assign_function_value
                    , ( bp::arg("other") )
                    , bp::return_self< >()
                    , "" );
            
            }
            const_iterator_exposer.def( bp::self == bp::self );
        }
        { //::SireBase::ChunkedVector< double, 100 >::iterator
            typedef bp::class_< SireBase::ChunkedVector< double, 100 >::iterator > iterator_exposer_t;
            iterator_exposer_t iterator_exposer = iterator_exposer_t( "iterator", "", bp::init< >("") );
            bp::scope iterator_scope( iterator_exposer );
            iterator_exposer.def( bp::init< SireBase::ChunkedVector< double, 100 >::iterator const & >(( bp::arg("other") ), "") );
            iterator_exposer.def( bp::self != bp::self );
            { //::SireBase::ChunkedVector< double, 100 >::iterator::operator=
            
                typedef ::SireBase::ChunkedVector< double, 100 >::iterator & ( ::SireBase::ChunkedVector< double, 100 >::iterator::*assign_function_type)( ::SireBase::ChunkedVector< double, 100 >::iterator const & ) ;
                assign_function_type assign_function_value( &::SireBase::ChunkedVector< double, 100 >::iterator::operator= );
                
                iterator_exposer.def( 
                    "assign"
                    , assign_function_value
                    , ( bp::arg("other") )
                    , bp::return_self< >()
                    , "" );
            
            }
            iterator_exposer.def( bp::self == bp::self );
        }
        ChunkedVector_double__exposer.def( bp::init< int >(( bp::arg("size") ), "") );
        ChunkedVector_double__exposer.def( bp::init< int, double const & >(( bp::arg("size"), bp::arg("value") ), "") );
        ChunkedVector_double__exposer.def( bp::init< SireBase::ChunkedVector< double, 100 > const & >(( bp::arg("other") ), "") );
        { //::SireBase::ChunkedVector< double, 100 >::append
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef void ( ::SireBase::ChunkedVector< double, 100 >::*append_function_type)( double const & ) ;
            append_function_type append_function_value( &::SireBase::ChunkedVector< double, 100 >::append );
            
            ChunkedVector_double__exposer.def( 
                "append"
                , append_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::at
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef double const & ( ::SireBase::ChunkedVector< double, 100 >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireBase::ChunkedVector< double, 100 >::at );
            
            ChunkedVector_double__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::capacity
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef int ( ::SireBase::ChunkedVector< double, 100 >::*capacity_function_type)(  ) const;
            capacity_function_type capacity_function_value( &::SireBase::ChunkedVector< double, 100 >::capacity );
            
            ChunkedVector_double__exposer.def( 
                "capacity"
                , capacity_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::clear
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef void ( ::SireBase::ChunkedVector< double, 100 >::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireBase::ChunkedVector< double, 100 >::clear );
            
            ChunkedVector_double__exposer.def( 
                "clear"
                , clear_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::count
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef int ( ::SireBase::ChunkedVector< double, 100 >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireBase::ChunkedVector< double, 100 >::count );
            
            ChunkedVector_double__exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::count
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef int ( ::SireBase::ChunkedVector< double, 100 >::*count_function_type)( double const & ) const;
            count_function_type count_function_value( &::SireBase::ChunkedVector< double, 100 >::count );
            
            ChunkedVector_double__exposer.def( 
                "count"
                , count_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::fromList
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef ::SireBase::ChunkedVector< double, 100 > ( *fromList_function_type )( ::QList< double > const & );
            fromList_function_type fromList_function_value( &::SireBase::ChunkedVector< double, 100 >::fromList );
            
            ChunkedVector_double__exposer.def( 
                "fromList"
                , fromList_function_value
                , ( bp::arg("list") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::fromStdVector
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef ::SireBase::ChunkedVector< double, 100 > ( *fromStdVector_function_type )( ::std::vector< double > const & );
            fromStdVector_function_type fromStdVector_function_value( &::SireBase::ChunkedVector< double, 100 >::fromStdVector );
            
            ChunkedVector_double__exposer.def( 
                "fromStdVector"
                , fromStdVector_function_value
                , ( bp::arg("vector") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::fromVector
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef ::SireBase::ChunkedVector< double, 100 > ( *fromVector_function_type )( ::QVector< double > const & );
            fromVector_function_type fromVector_function_value( &::SireBase::ChunkedVector< double, 100 >::fromVector );
            
            ChunkedVector_double__exposer.def( 
                "fromVector"
                , fromVector_function_value
                , ( bp::arg("vector") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::isEmpty
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef bool ( ::SireBase::ChunkedVector< double, 100 >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireBase::ChunkedVector< double, 100 >::isEmpty );
            
            ChunkedVector_double__exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ChunkedVector_double__exposer.def( bp::self != bp::self );
        { //::SireBase::ChunkedVector< double, 100 >::operator=
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef ::SireBase::ChunkedVector< double, 100 > & ( ::SireBase::ChunkedVector< double, 100 >::*assign_function_type)( ::SireBase::ChunkedVector< double, 100 > const & ) ;
            assign_function_type assign_function_value( &::SireBase::ChunkedVector< double, 100 >::operator= );
            
            ChunkedVector_double__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ChunkedVector_double__exposer.def( bp::self == bp::self );
        { //::SireBase::ChunkedVector< double, 100 >::operator[]
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef double & ( ::SireBase::ChunkedVector< double, 100 >::*__getitem___function_type)( int ) ;
            __getitem___function_type __getitem___function_value( &::SireBase::ChunkedVector< double, 100 >::operator[] );
            
            ChunkedVector_double__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_non_const_reference >()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::operator[]
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef double const & ( ::SireBase::ChunkedVector< double, 100 >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::ChunkedVector< double, 100 >::operator[] );
            
            ChunkedVector_double__exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::remove
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef void ( ::SireBase::ChunkedVector< double, 100 >::*remove_function_type)( int ) ;
            remove_function_type remove_function_value( &::SireBase::ChunkedVector< double, 100 >::remove );
            
            ChunkedVector_double__exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::remove
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef void ( ::SireBase::ChunkedVector< double, 100 >::*remove_function_type)( int,int ) ;
            remove_function_type remove_function_value( &::SireBase::ChunkedVector< double, 100 >::remove );
            
            ChunkedVector_double__exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("i"), bp::arg("count") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::reserve
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef void ( ::SireBase::ChunkedVector< double, 100 >::*reserve_function_type)( int ) ;
            reserve_function_type reserve_function_value( &::SireBase::ChunkedVector< double, 100 >::reserve );
            
            ChunkedVector_double__exposer.def( 
                "reserve"
                , reserve_function_value
                , ( bp::arg("count") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::resize
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef void ( ::SireBase::ChunkedVector< double, 100 >::*resize_function_type)( int ) ;
            resize_function_type resize_function_value( &::SireBase::ChunkedVector< double, 100 >::resize );
            
            ChunkedVector_double__exposer.def( 
                "resize"
                , resize_function_value
                , ( bp::arg("count") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::size
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef int ( ::SireBase::ChunkedVector< double, 100 >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireBase::ChunkedVector< double, 100 >::size );
            
            ChunkedVector_double__exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::squeeze
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef void ( ::SireBase::ChunkedVector< double, 100 >::*squeeze_function_type)(  ) ;
            squeeze_function_type squeeze_function_value( &::SireBase::ChunkedVector< double, 100 >::squeeze );
            
            ChunkedVector_double__exposer.def( 
                "squeeze"
                , squeeze_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::toList
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef ::QList< double > ( ::SireBase::ChunkedVector< double, 100 >::*toList_function_type)(  ) const;
            toList_function_type toList_function_value( &::SireBase::ChunkedVector< double, 100 >::toList );
            
            ChunkedVector_double__exposer.def( 
                "toList"
                , toList_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::toStdVector
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef ::std::vector< double > ( ::SireBase::ChunkedVector< double, 100 >::*toStdVector_function_type)(  ) const;
            toStdVector_function_type toStdVector_function_value( &::SireBase::ChunkedVector< double, 100 >::toStdVector );
            
            ChunkedVector_double__exposer.def( 
                "toStdVector"
                , toStdVector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::toVector
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef ::QVector< double > ( ::SireBase::ChunkedVector< double, 100 >::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireBase::ChunkedVector< double, 100 >::toVector );
            
            ChunkedVector_double__exposer.def( 
                "toVector"
                , toVector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::value
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef double ( ::SireBase::ChunkedVector< double, 100 >::*value_function_type)( int ) const;
            value_function_type value_function_value( &::SireBase::ChunkedVector< double, 100 >::value );
            
            ChunkedVector_double__exposer.def( 
                "value"
                , value_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::ChunkedVector< double, 100 >::value
        
            typedef SireBase::ChunkedVector< double, 100 > exported_class_t;
            typedef double ( ::SireBase::ChunkedVector< double, 100 >::*value_function_type)( int,double const & ) const;
            value_function_type value_function_value( &::SireBase::ChunkedVector< double, 100 >::value );
            
            ChunkedVector_double__exposer.def( 
                "value"
                , value_function_value
                , ( bp::arg("i"), bp::arg("default_value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        ChunkedVector_double__exposer.staticmethod( "fromList" );
        ChunkedVector_double__exposer.staticmethod( "fromStdVector" );
        ChunkedVector_double__exposer.staticmethod( "fromVector" );
        ChunkedVector_double__exposer.def( "__copy__", &__copy__<SireBase::ChunkedVector<double, 100>>);
        ChunkedVector_double__exposer.def( "__deepcopy__", &__copy__<SireBase::ChunkedVector<double, 100>>);
        ChunkedVector_double__exposer.def( "clone", &__copy__<SireBase::ChunkedVector<double, 100>>);
        ChunkedVector_double__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::ChunkedVector<double, 100> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChunkedVector_double__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::ChunkedVector<double, 100> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChunkedVector_double__exposer.def_pickle(sire_pickle_suite< ::SireBase::ChunkedVector<double, 100> >());
        ChunkedVector_double__exposer.def( "__str__", &pvt_get_name);
        ChunkedVector_double__exposer.def( "__repr__", &pvt_get_name);
        ChunkedVector_double__exposer.def( "__len__", &__len_size< ::SireBase::ChunkedVector<double, 100> > );
    }

}
