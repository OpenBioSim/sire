// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "InterGroupSoftCLJFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intersoftcljff.h"

#include "intersoftcljff.h"

SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > __copy__(const SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > &other){ return SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> >(other); }

#include "Helpers/copy.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_InterGroupSoftCLJFF_class(){

    { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >
        typedef bp::class_< SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >, bp::bases< SireFF::FF3D, SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>, SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential>, SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InterGroupSoftCLJFF_exposer_t;
        InterGroupSoftCLJFF_exposer_t InterGroupSoftCLJFF_exposer = InterGroupSoftCLJFF_exposer_t( "InterGroupSoftCLJFF", "", bp::init< >("") );
        bp::scope InterGroupSoftCLJFF_scope( InterGroupSoftCLJFF_exposer );
        InterGroupSoftCLJFF_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        InterGroupSoftCLJFF_exposer.def( bp::init< SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*energy_function_type)(  ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy );
            
            InterGroupSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*energy_function_type)( ::SireCAS::Symbol const & ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy );
            
            InterGroupSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("component") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy );
            
            InterGroupSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::energy );
            
            InterGroupSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterGroupSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterGroupSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterGroupSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::field );
            
            InterGroupSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force );
            
            InterGroupSoftCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::force );
            
            InterGroupSoftCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        InterGroupSoftCLJFF_exposer.def( bp::self != bp::self );
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::operator=
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > & ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*assign_function_type)( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::operator= );
            
            InterGroupSoftCLJFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        InterGroupSoftCLJFF_exposer.def( bp::self == bp::self );
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterGroupSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterGroupSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterGroupSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::potential );
            
            InterGroupSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::typeName
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::typeName );
            
            InterGroupSoftCLJFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::what
        
            typedef SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Inter2B2G3DFF< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > >::what );
            
            InterGroupSoftCLJFF_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        InterGroupSoftCLJFF_exposer.staticmethod( "typeName" );
        InterGroupSoftCLJFF_exposer.def( "__copy__", &__copy__<SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> >>);
        InterGroupSoftCLJFF_exposer.def( "__deepcopy__", &__copy__<SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> >>);
        InterGroupSoftCLJFF_exposer.def( "clone", &__copy__<SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> >>);
        InterGroupSoftCLJFF_exposer.def( "__str__", &__str__< ::SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
        InterGroupSoftCLJFF_exposer.def( "__repr__", &__str__< ::SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
        InterGroupSoftCLJFF_exposer.def( "__len__", &__len_count< ::SireFF::Inter2B2G3DFF<SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential> > > );
    }

}
