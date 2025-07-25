// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "DihedralPerturbation.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/identities.h"

#include "SireCAS/values.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "core.h"

#include "geometryperturbation.h"

#include "molecule.h"

#include "moleditor.h"

#include "mover.hpp"

#include "geometryperturbation.h"

SireMol::DihedralPerturbation __copy__(const SireMol::DihedralPerturbation &other){ return SireMol::DihedralPerturbation(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_DihedralPerturbation_class(){

    { //::SireMol::DihedralPerturbation
        typedef bp::class_< SireMol::DihedralPerturbation, bp::bases< SireMol::GeometryPerturbation, SireMol::Perturbation, SireBase::Property > > DihedralPerturbation_exposer_t;
        DihedralPerturbation_exposer_t DihedralPerturbation_exposer = DihedralPerturbation_exposer_t( "DihedralPerturbation", "This perturbation moves a dihedral between two sizes.\n\nThis uses the anchors property to anchor parts\nof the molecule, the weight function property\nto weight the motion of the parts of the molecule,\nand the coordinates property to get the coordinates\nto move\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope DihedralPerturbation_scope( DihedralPerturbation_exposer );
        DihedralPerturbation_exposer.def( bp::init< SireMol::DihedralID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("dihedral"), bp::arg("start"), bp::arg("end"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the dihedral dihedral from start to end") );
        DihedralPerturbation_exposer.def( bp::init< SireMol::DihedralID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("dihedral"), bp::arg("start"), bp::arg("end"), bp::arg("mapping_function"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the dihedral dihedral from start to end\nusing the passed mapping function") );
        DihedralPerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("start"), bp::arg("end"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the dihedral between atoms atom0, atom1, atom2 and atom3\nfrom start to end") );
        DihedralPerturbation_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireUnits::Dimension::Angle const &, SireUnits::Dimension::Angle const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3"), bp::arg("start"), bp::arg("end"), bp::arg("mapping_function"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to perturb the dihedral between atoms atom0, atom1, atom2 and atom3\nfrom start to end using the passed mapping function") );
        DihedralPerturbation_exposer.def( bp::init< SireMol::DihedralPerturbation const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::DihedralPerturbation::dihedral
        
            typedef ::SireMol::DihedralID const & ( ::SireMol::DihedralPerturbation::*dihedral_function_type)(  ) const;
            dihedral_function_type dihedral_function_value( &::SireMol::DihedralPerturbation::dihedral );
            
            DihedralPerturbation_exposer.def( 
                "dihedral"
                , dihedral_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ID that identifies that dihedral that will be perturbed" );
        
        }
        { //::SireMol::DihedralPerturbation::end
        
            typedef ::SireUnits::Dimension::Angle const & ( ::SireMol::DihedralPerturbation::*end_function_type)(  ) const;
            end_function_type end_function_value( &::SireMol::DihedralPerturbation::end );
            
            DihedralPerturbation_exposer.def( 
                "end"
                , end_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the end length of the dihedral" );
        
        }
        DihedralPerturbation_exposer.def( bp::self != bp::self );
        { //::SireMol::DihedralPerturbation::operator=
        
            typedef ::SireMol::DihedralPerturbation & ( ::SireMol::DihedralPerturbation::*assign_function_type)( ::SireMol::DihedralPerturbation const & ) ;
            assign_function_type assign_function_value( &::SireMol::DihedralPerturbation::operator= );
            
            DihedralPerturbation_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        DihedralPerturbation_exposer.def( bp::self == bp::self );
        { //::SireMol::DihedralPerturbation::start
        
            typedef ::SireUnits::Dimension::Angle const & ( ::SireMol::DihedralPerturbation::*start_function_type)(  ) const;
            start_function_type start_function_value( &::SireMol::DihedralPerturbation::start );
            
            DihedralPerturbation_exposer.def( 
                "start"
                , start_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the start length of the dihedral" );
        
        }
        { //::SireMol::DihedralPerturbation::toString
        
            typedef ::QString ( ::SireMol::DihedralPerturbation::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::DihedralPerturbation::toString );
            
            DihedralPerturbation_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::DihedralPerturbation::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::DihedralPerturbation::typeName );
            
            DihedralPerturbation_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::DihedralPerturbation::wouldChange
        
            typedef bool ( ::SireMol::DihedralPerturbation::*wouldChange_function_type)( ::SireMol::Molecule const &,::SireCAS::Values const & ) const;
            wouldChange_function_type wouldChange_function_value( &::SireMol::DihedralPerturbation::wouldChange );
            
            DihedralPerturbation_exposer.def( 
                "wouldChange"
                , wouldChange_function_value
                , ( bp::arg("molecule"), bp::arg("values") )
                , bp::release_gil_policy()
                , "Return whether or not this perturbation with the passed values would\nchange the molecule molecule" );
        
        }
        DihedralPerturbation_exposer.staticmethod( "typeName" );
        DihedralPerturbation_exposer.def( "__copy__", &__copy__<SireMol::DihedralPerturbation>);
        DihedralPerturbation_exposer.def( "__deepcopy__", &__copy__<SireMol::DihedralPerturbation>);
        DihedralPerturbation_exposer.def( "clone", &__copy__<SireMol::DihedralPerturbation>);
        DihedralPerturbation_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::DihedralPerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralPerturbation_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::DihedralPerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DihedralPerturbation_exposer.def_pickle(sire_pickle_suite< ::SireMol::DihedralPerturbation >());
        DihedralPerturbation_exposer.def( "__str__", &__str__< ::SireMol::DihedralPerturbation > );
        DihedralPerturbation_exposer.def( "__repr__", &__str__< ::SireMol::DihedralPerturbation > );
    }

}
