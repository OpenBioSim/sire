// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ExcludedPairs.pypp.hpp"

namespace bp = boost::python;

#include "SireMM/cljnbpairs.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "excludedpairs.h"

#include "excludedpairs.h"

SireMM::ExcludedPairs __copy__(const SireMM::ExcludedPairs &other){ return SireMM::ExcludedPairs(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_ExcludedPairs_class(){

    { //::SireMM::ExcludedPairs
        typedef bp::class_< SireMM::ExcludedPairs, bp::bases< SireMol::MolViewProperty, SireBase::Property > > ExcludedPairs_exposer_t;
        ExcludedPairs_exposer_t ExcludedPairs_exposer = ExcludedPairs_exposer_t( "ExcludedPairs", "This class holds the excluded atoms pair list for a molecule.\nExcluded atom pairs means that the non-bonded coulomb and\nLJ energy between these pair of atoms are not calculated.\n", bp::init< >("") );
        bp::scope ExcludedPairs_scope( ExcludedPairs_exposer );
        ExcludedPairs_exposer.def( bp::init< SireMol::MoleculeView const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molecule"), bp::arg("map")=SireBase::PropertyMap() ), "") );
        ExcludedPairs_exposer.def( bp::init< SireMM::ExcludedPairs const & >(( bp::arg("other") ), "") );
        { //::SireMM::ExcludedPairs::areExcluded
        
            typedef bool ( ::SireMM::ExcludedPairs::*areExcluded_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const & ) const;
            areExcluded_function_type areExcluded_function_value( &::SireMM::ExcludedPairs::areExcluded );
            
            ExcludedPairs_exposer.def( 
                "areExcluded"
                , areExcluded_function_value
                , ( bp::arg("atom0"), bp::arg("atom1") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::count
        
            typedef int ( ::SireMM::ExcludedPairs::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMM::ExcludedPairs::count );
            
            ExcludedPairs_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::info
        
            typedef ::SireMol::MoleculeInfo ( ::SireMM::ExcludedPairs::*info_function_type)(  ) const;
            info_function_type info_function_value( &::SireMM::ExcludedPairs::info );
            
            ExcludedPairs_exposer.def( 
                "info"
                , info_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::isCompatibleWith
        
            typedef bool ( ::SireMM::ExcludedPairs::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMM::ExcludedPairs::isCompatibleWith );
            
            ExcludedPairs_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::merge
        
            typedef ::SireBase::PropertyList ( ::SireMM::ExcludedPairs::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMM::ExcludedPairs::merge );
            
            ExcludedPairs_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "Merge this property with another property" );
        
        }
        { //::SireMM::ExcludedPairs::nExcludedPairs
        
            typedef int ( ::SireMM::ExcludedPairs::*nExcludedPairs_function_type)(  ) const;
            nExcludedPairs_function_type nExcludedPairs_function_value( &::SireMM::ExcludedPairs::nExcludedPairs );
            
            ExcludedPairs_exposer.def( 
                "nExcludedPairs"
                , nExcludedPairs_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ExcludedPairs_exposer.def( bp::self != bp::self );
        { //::SireMM::ExcludedPairs::operator=
        
            typedef ::SireMM::ExcludedPairs & ( ::SireMM::ExcludedPairs::*assign_function_type)( ::SireMM::ExcludedPairs const & ) ;
            assign_function_type assign_function_value( &::SireMM::ExcludedPairs::operator= );
            
            ExcludedPairs_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ExcludedPairs_exposer.def( bp::self == bp::self );
        { //::SireMM::ExcludedPairs::operator[]
        
            typedef ::std::tuple< SireMol::AtomIdx, SireMol::AtomIdx > ( ::SireMM::ExcludedPairs::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMM::ExcludedPairs::operator[] );
            
            ExcludedPairs_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::setExcluded
        
            typedef void ( ::SireMM::ExcludedPairs::*setExcluded_function_type)( ::SireMol::AtomID const &,::SireMol::AtomID const &,bool ) ;
            setExcluded_function_type setExcluded_function_value( &::SireMM::ExcludedPairs::setExcluded );
            
            ExcludedPairs_exposer.def( 
                "setExcluded"
                , setExcluded_function_value
                , ( bp::arg("atom0"), bp::arg("atom1"), bp::arg("are_excluded") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::toString
        
            typedef ::QString ( ::SireMM::ExcludedPairs::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::ExcludedPairs::toString );
            
            ExcludedPairs_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::ExcludedPairs::typeName );
            
            ExcludedPairs_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::updateBondMatrix
        
            typedef void ( ::SireMM::ExcludedPairs::*updateBondMatrix_function_type)( ::QVector< QVector< bool > > & ) const;
            updateBondMatrix_function_type updateBondMatrix_function_value( &::SireMM::ExcludedPairs::updateBondMatrix );
            
            ExcludedPairs_exposer.def( 
                "updateBondMatrix"
                , updateBondMatrix_function_value
                , ( bp::arg("bond_matrix") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::ExcludedPairs::what
        
            typedef char const * ( ::SireMM::ExcludedPairs::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::ExcludedPairs::what );
            
            ExcludedPairs_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ExcludedPairs_exposer.staticmethod( "typeName" );
        ExcludedPairs_exposer.def( "__copy__", &__copy__<SireMM::ExcludedPairs>);
        ExcludedPairs_exposer.def( "__deepcopy__", &__copy__<SireMM::ExcludedPairs>);
        ExcludedPairs_exposer.def( "clone", &__copy__<SireMM::ExcludedPairs>);
        ExcludedPairs_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::ExcludedPairs >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExcludedPairs_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::ExcludedPairs >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExcludedPairs_exposer.def_pickle(sire_pickle_suite< ::SireMM::ExcludedPairs >());
        ExcludedPairs_exposer.def( "__str__", &__str__< ::SireMM::ExcludedPairs > );
        ExcludedPairs_exposer.def( "__repr__", &__str__< ::SireMM::ExcludedPairs > );
        ExcludedPairs_exposer.def( "__len__", &__len_count< ::SireMM::ExcludedPairs > );
    }

}
