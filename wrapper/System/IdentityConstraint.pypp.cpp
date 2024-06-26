// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "IdentityConstraint.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/refcountdata.h"

#include "SireBase/shareddatapointer.hpp"

#include "SireMaths/linearap.h"

#include "SireMaths/nmatrix.h"

#include "SireMaths/nvector.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/molecules.h"

#include "SireMol/moleditor.h"

#include "SireMol/viewsofmol.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/errors.h"

#include "SireVol/space.h"

#include "closemols.h"

#include "delta.h"

#include "identityconstraint.h"

#include "system.h"

#include <QDebug>

#include <QElapsedTimer>

#include <QVarLengthArray>

#include "identityconstraint.h"

SireSystem::IdentityConstraint __copy__(const SireSystem::IdentityConstraint &other){ return SireSystem::IdentityConstraint(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_IdentityConstraint_class(){

    { //::SireSystem::IdentityConstraint
        typedef bp::class_< SireSystem::IdentityConstraint, bp::bases< SireSystem::MoleculeConstraint, SireSystem::Constraint, SireBase::Property > > IdentityConstraint_exposer_t;
        IdentityConstraint_exposer_t IdentityConstraint_exposer = IdentityConstraint_exposer_t( "IdentityConstraint", "An identity constraint provides a method of constraining\nthe identity of molecules based on where they are located.\n\nFor example, it can be useful to be able to identify a\nwater in a binding pocket. However, in a normal simulation,\nwe dont identify waters by location, but by their index\n(e.g. this is the first water, this is the second etc.).\n\nThis means that the identity of the water in the binding\npocket can change, e.g. it can start with the fifth water\nin the pocket, but during the simulation the fifth water\nmay diffuse out of the pocket, and the twentieth water\nwould diffuse in its place. The identity of the water\nin the binding pocket will thus have changed from the\nfifth water to the twentieth water.\n\nAn identity constraint works by constantly monitoring\nthe locations of the waters, and so it can detect when\nthe fifth water is displaced by the twentieth water. When\nit detects that this has occured, the constraint swaps\nthe coordinates of the fifth and twentieth waters,\nthereby ensuring that the fifth water stays in the pocket.\nThis doesnt affect the energy or the statistics of the\nsystem, as waters are indistinguishable (there are N\nequivalent configurations of N waters - we only see them\nas N different configurations as we identify each water,\nwhen really they are indistinguishable).\n\nThe idea of constraining the identity of molecules\nwas first presented by M. Tyka, R. Sessions and A. Clarke in\n\nAbsolute Free-Energy Calculations of Liquids Using\na Harmonic Reference State\n\nJ. Chem. Phys. B,  2007, 111, 9571-9580\n\ndoi:10.1021jp072357w\n\nThey used the method to constrain the identity of all\nmolecules, so that harmonic restraints can be applied\nto them all.\n\nThis identity constraint is more general, and allows\nthe identification of a subset of molecules to be\nconstrained (e.g. just the waters in a binding pocket).\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope IdentityConstraint_scope( IdentityConstraint_exposer );
        IdentityConstraint_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct the constraint that constrains the identities of all\nof the molecules in the passed molecule group. This uses the current\nlocations of the molecules to apply the constraint. The (optionally\nsupplied) property map is used to find the properties required\nof this constraint") );
        IdentityConstraint_exposer.def( bp::init< SireFF::PointRef const &, SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("point"), bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct the constraint that constrains the identity of a single\nmolecule in the passed molecule group - this sets the identity\nof the first molecule to be that of the one closest to the\npassed point. The (optionally supplied) property map is used to\nfind the properties required of this constraint") );
        IdentityConstraint_exposer.def( bp::init< QList< SireBase::PropPtr< SireFF::Point > > const &, SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("points"), bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct the constraint that constrains the identities of the\npoints.count() molecules from the passed molecule group so that\nthe first molecule is identified by the first point, the second\nmolecule is identified by the second point, and the nth molecule\nis identified by the nth point. The (optionally supplied) property\nmap is used to find the properties required of this constraint") );
        IdentityConstraint_exposer.def( bp::init< QVector< SireBase::PropPtr< SireFF::Point > > const &, SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("points"), bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct the constraint that constrains the identities of the\npoints.count() molecules from the passed molecule group so that\nthe first molecule is identified by the first point, the second\nmolecule is identified by the second point, and the nth molecule\nis identified by the nth point. The (optionally supplied) property\nmap is used to find the properties required of this constraint") );
        IdentityConstraint_exposer.def( bp::init< SireSystem::IdentityConstraint const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::IdentityConstraint::constrain
        
            typedef ::SireMol::MolGroupPtr ( *constrain_function_type )( ::SireMol::MoleculeGroup const &,::SireFF::PointRef const &,::SireBase::PropertyMap const & );
            constrain_function_type constrain_function_value( &::SireSystem::IdentityConstraint::constrain );
            
            IdentityConstraint_exposer.def( 
                "constrain"
                , constrain_function_value
                , ( bp::arg("molgroup"), bp::arg("point"), bp::arg("map")=SireBase::PropertyMap() )
                , "Static function used to constrain the identities of the molecules\nin molgroup against the point point. This makes the first molecule\nin the group have the identity that matches this point" );
        
        }
        { //::SireSystem::IdentityConstraint::constrain
        
            typedef ::SireMol::MolGroupPtr ( *constrain_function_type )( ::SireMol::MoleculeGroup const &,::QVector< SireBase::PropPtr< SireFF::Point > > const &,::SireBase::PropertyMap const & );
            constrain_function_type constrain_function_value( &::SireSystem::IdentityConstraint::constrain );
            
            IdentityConstraint_exposer.def( 
                "constrain"
                , constrain_function_value
                , ( bp::arg("molgroup"), bp::arg("points"), bp::arg("map")=SireBase::PropertyMap() )
                , "Static function used to constrain the identities of the molecules\nin molgroup against the identity points in points - the\nfirst npoints molecules in the group are constrained in order\nagainst the points" );
        
        }
        { //::SireSystem::IdentityConstraint::constrain
        
            typedef ::SireMol::MolGroupPtr ( *constrain_function_type )( ::SireMol::MoleculeGroup const &,::QList< SireBase::PropPtr< SireFF::Point > > const &,::SireBase::PropertyMap const & );
            constrain_function_type constrain_function_value( &::SireSystem::IdentityConstraint::constrain );
            
            IdentityConstraint_exposer.def( 
                "constrain"
                , constrain_function_value
                , ( bp::arg("molgroup"), bp::arg("points"), bp::arg("map")=SireBase::PropertyMap() )
                , "Static function used to constrain the identities of the molecules\nin molgroup against the identity points in points - the\nfirst npoints molecules in the group are constrained in order\nagainst the points" );
        
        }
        { //::SireSystem::IdentityConstraint::moleculeGroup
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireSystem::IdentityConstraint::*moleculeGroup_function_type)(  ) const;
            moleculeGroup_function_type moleculeGroup_function_value( &::SireSystem::IdentityConstraint::moleculeGroup );
            
            IdentityConstraint_exposer.def( 
                "moleculeGroup"
                , moleculeGroup_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the molecule group acted on by this constraint" );
        
        }
        IdentityConstraint_exposer.def( bp::self != bp::self );
        { //::SireSystem::IdentityConstraint::operator=
        
            typedef ::SireSystem::IdentityConstraint & ( ::SireSystem::IdentityConstraint::*assign_function_type)( ::SireSystem::IdentityConstraint const & ) ;
            assign_function_type assign_function_value( &::SireSystem::IdentityConstraint::operator= );
            
            IdentityConstraint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IdentityConstraint_exposer.def( bp::self == bp::self );
        { //::SireSystem::IdentityConstraint::points
        
            typedef ::QVector< SireBase::PropPtr< SireFF::Point > > ( ::SireSystem::IdentityConstraint::*points_function_type)(  ) const;
            points_function_type points_function_value( &::SireSystem::IdentityConstraint::points );
            
            IdentityConstraint_exposer.def( 
                "points"
                , points_function_value
                , bp::release_gil_policy()
                , "Return the points used to identify the molecules" );
        
        }
        { //::SireSystem::IdentityConstraint::propertyMap
        
            typedef ::SireBase::PropertyMap const & ( ::SireSystem::IdentityConstraint::*propertyMap_function_type)(  ) const;
            propertyMap_function_type propertyMap_function_value( &::SireSystem::IdentityConstraint::propertyMap );
            
            IdentityConstraint_exposer.def( 
                "propertyMap"
                , propertyMap_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the property map used to find the properties used\nby this constraint" );
        
        }
        { //::SireSystem::IdentityConstraint::toString
        
            typedef ::QString ( ::SireSystem::IdentityConstraint::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireSystem::IdentityConstraint::toString );
            
            IdentityConstraint_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this constraint" );
        
        }
        { //::SireSystem::IdentityConstraint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::IdentityConstraint::typeName );
            
            IdentityConstraint_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireSystem::IdentityConstraint::useFewPointsAlgorithm
        
            typedef void ( ::SireSystem::IdentityConstraint::*useFewPointsAlgorithm_function_type)(  ) ;
            useFewPointsAlgorithm_function_type useFewPointsAlgorithm_function_value( &::SireSystem::IdentityConstraint::useFewPointsAlgorithm );
            
            IdentityConstraint_exposer.def( 
                "useFewPointsAlgorithm"
                , useFewPointsAlgorithm_function_value
                , bp::release_gil_policy()
                , "Function used for debugging that switches this object over\nto using the few points algorithm to apply the constraint" );
        
        }
        { //::SireSystem::IdentityConstraint::useManyPointsAlgorithm
        
            typedef void ( ::SireSystem::IdentityConstraint::*useManyPointsAlgorithm_function_type)(  ) ;
            useManyPointsAlgorithm_function_type useManyPointsAlgorithm_function_value( &::SireSystem::IdentityConstraint::useManyPointsAlgorithm );
            
            IdentityConstraint_exposer.def( 
                "useManyPointsAlgorithm"
                , useManyPointsAlgorithm_function_value
                , bp::release_gil_policy()
                , "Function used for debugging that switches this object over\nto using the many points algorithm to apply the constraint" );
        
        }
        { //::SireSystem::IdentityConstraint::useSinglePointAlgorithm
        
            typedef void ( ::SireSystem::IdentityConstraint::*useSinglePointAlgorithm_function_type)(  ) ;
            useSinglePointAlgorithm_function_type useSinglePointAlgorithm_function_value( &::SireSystem::IdentityConstraint::useSinglePointAlgorithm );
            
            IdentityConstraint_exposer.def( 
                "useSinglePointAlgorithm"
                , useSinglePointAlgorithm_function_value
                , bp::release_gil_policy()
                , "Function used for debugging that switches this object over\nto using the single point algorithm to apply the constraint\nThrow: SireError::invalid_state\n" );
        
        }
        IdentityConstraint_exposer.staticmethod( "constrain" );
        IdentityConstraint_exposer.staticmethod( "typeName" );
        IdentityConstraint_exposer.def( "__copy__", &__copy__<SireSystem::IdentityConstraint>);
        IdentityConstraint_exposer.def( "__deepcopy__", &__copy__<SireSystem::IdentityConstraint>);
        IdentityConstraint_exposer.def( "clone", &__copy__<SireSystem::IdentityConstraint>);
        IdentityConstraint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::IdentityConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IdentityConstraint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::IdentityConstraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IdentityConstraint_exposer.def_pickle(sire_pickle_suite< ::SireSystem::IdentityConstraint >());
        IdentityConstraint_exposer.def( "__str__", &__str__< ::SireSystem::IdentityConstraint > );
        IdentityConstraint_exposer.def( "__repr__", &__str__< ::SireSystem::IdentityConstraint > );
    }

}
