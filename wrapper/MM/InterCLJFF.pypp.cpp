// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "InterCLJFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intercljff.h"

#include "intercljff.h"

SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > __copy__(const SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > &other){ return SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >(other); }

#include "Helpers/copy.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_InterCLJFF_class(){

    { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >
        typedef bp::class_< SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >, bp::bases< SireFF::FF3D, SireMM::CLJPotentialInterface<SireMM::InterCLJPotential>, SireFF::G1FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > InterCLJFF_exposer_t;
        InterCLJFF_exposer_t InterCLJFF_exposer = InterCLJFF_exposer_t( "InterCLJFF", "", bp::init< >("") );
        bp::scope InterCLJFF_scope( InterCLJFF_exposer );
        InterCLJFF_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        InterCLJFF_exposer.def( bp::init< SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*energy_function_type)(  ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy );
            
            InterCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*energy_function_type)( ::SireCAS::Symbol const & ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy );
            
            InterCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("component") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy );
            
            InterCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*energy_function_type)( ::SireFF::EnergyTable &,::SireCAS::Symbol const &,double ) ;
            energy_function_type energy_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::energy );
            
            InterCLJFF_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("energytable"), bp::arg("symbol"), bp::arg("scale_energy")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field );
            
            InterCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field );
            
            InterCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field );
            
            InterCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*field_function_type)( ::SireFF::FieldTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            field_function_type field_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::field );
            
            InterCLJFF_exposer.def( 
                "field"
                , field_function_value
                , ( bp::arg("fieldtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_field")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::force );
            
            InterCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::force
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*force_function_type)( ::SireFF::ForceTable &,::SireCAS::Symbol const &,double ) ;
            force_function_type force_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::force );
            
            InterCLJFF_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("symbol"), bp::arg("scale_force")=1 )
                , "" );
        
        }
        InterCLJFF_exposer.def( bp::self != bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::operator=
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > & ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*assign_function_type)( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::operator= );
            
            InterCLJFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        InterCLJFF_exposer.def( bp::self == bp::self );
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::packCoordinates
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*packCoordinates_function_type)(  ) ;
            packCoordinates_function_type packCoordinates_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::packCoordinates );
            
            InterCLJFF_exposer.def( 
                "packCoordinates"
                , packCoordinates_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential );
            
            InterCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential );
            
            InterCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential );
            
            InterCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*potential_function_type)( ::SireFF::PotentialTable &,::SireCAS::Symbol const &,::SireFF::Probe const &,double ) ;
            potential_function_type potential_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::potential );
            
            InterCLJFF_exposer.def( 
                "potential"
                , potential_function_value
                , ( bp::arg("potentialtable"), bp::arg("component"), bp::arg("probe"), bp::arg("scale_potential")=1 )
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::typeName
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::typeName );
            
            InterCLJFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::what
        
            typedef SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Inter2B3DFF< SireMM::CLJPotentialInterface< SireMM::InterCLJPotential > >::what );
            
            InterCLJFF_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        InterCLJFF_exposer.staticmethod( "typeName" );
        InterCLJFF_exposer.def( "__copy__", &__copy__<SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >>);
        InterCLJFF_exposer.def( "__deepcopy__", &__copy__<SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >>);
        InterCLJFF_exposer.def( "clone", &__copy__<SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> >>);
        InterCLJFF_exposer.def( "__str__", &__str__< ::SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
        InterCLJFF_exposer.def( "__repr__", &__str__< ::SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
        InterCLJFF_exposer.def( "__len__", &__len_count< ::SireFF::Inter2B3DFF<SireMM::CLJPotentialInterface<SireMM::InterCLJPotential> > > );
    }

}
