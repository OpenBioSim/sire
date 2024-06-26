// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "IntraSoftCLJFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intrasoftcljff.h"

#include "intrasoftcljff.h"

SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> > __copy__(const SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> > &other){ return SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >(other); }

#include "Helpers/copy.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_IntraSoftCLJFF_class(){

    { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >
        typedef bp::class_< SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >, bp::bases< SireFF::FF3D, SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential>, SireMM::CLJPotentialInterface<SireMM::IntraSoftCLJPotential>, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > IntraSoftCLJFF_exposer_t;
        IntraSoftCLJFF_exposer_t IntraSoftCLJFF_exposer = IntraSoftCLJFF_exposer_t( "IntraSoftCLJFF", "", bp::init< >("") );
        bp::scope IntraSoftCLJFF_scope( IntraSoftCLJFF_exposer );
        IntraSoftCLJFF_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        IntraSoftCLJFF_exposer.def( bp::init< SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*energy_function_type)(  ) ;
            energy_function_type energy_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy );
            
            IntraSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*energy_function_type)( ::SireCAS::Symbol const & ) ;
            energy_function_type energy_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy );
            
            IntraSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("component") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy );
            
            IntraSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::energy );
            
            IntraSoftCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field );
            
            IntraSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field );
            
            IntraSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field );
            
            IntraSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::field );
            
            IntraSoftCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::force
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::force );
            
            IntraSoftCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::force
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::force );
            
            IntraSoftCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        IntraSoftCLJFF_exposer.def( bp::self != bp::self );
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::operator=
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > & ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*assign_function_type)( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::operator= );
            
            IntraSoftCLJFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IntraSoftCLJFF_exposer.def( bp::self == bp::self );
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential );
            
            IntraSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential );
            
            IntraSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential );
            
            IntraSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::potential );
            
            IntraSoftCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::typeName
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::typeName );
            
            IntraSoftCLJFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::what
        
            typedef SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Intra2B3DFF< SireMM::SoftCLJPotentialInterface< SireMM::IntraSoftCLJPotential > >::what );
            
            IntraSoftCLJFF_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IntraSoftCLJFF_exposer.staticmethod( "typeName" );
        IntraSoftCLJFF_exposer.def( "__copy__", &__copy__<SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >>);
        IntraSoftCLJFF_exposer.def( "__deepcopy__", &__copy__<SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >>);
        IntraSoftCLJFF_exposer.def( "clone", &__copy__<SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> >>);
        IntraSoftCLJFF_exposer.def( "__str__", &__str__< ::SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> > > );
        IntraSoftCLJFF_exposer.def( "__repr__", &__str__< ::SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> > > );
        IntraSoftCLJFF_exposer.def( "__len__", &__len_count< ::SireFF::Intra2B3DFF<SireMM::SoftCLJPotentialInterface<SireMM::IntraSoftCLJPotential> > > );
    }

}
