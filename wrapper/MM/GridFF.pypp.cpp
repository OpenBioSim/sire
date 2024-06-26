// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "GridFF.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/constants.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "atomljs.h"

#include "cljpotential.h"

#include "gridff.h"

#include <QDebug>

#include <QElapsedTimer>

#include <QTime>

#include "gridff.h"

SireMM::GridFF __copy__(const SireMM::GridFF &other){ return SireMM::GridFF(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_GridFF_class(){

    { //::SireMM::GridFF
        typedef bp::class_< SireMM::GridFF, bp::bases< SireFF::FF3D, SireMM::CLJPotentialInterface<SireMM::InterCLJPotential>, SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > GridFF_exposer_t;
        GridFF_exposer_t GridFF_exposer = GridFF_exposer_t( "GridFF", "This class calculates the coulomb and LJ energy between\nall molecules in group 1 and all molecules in group 2.\nThe calculation is optimised, as the molecules in group 2\nare represented using a grid. This is ideal for situations\nwhere the molecules on group 2 move little, or not at all.\n\nAuthor: Christopher Woods\n", bp::init< >("Empty constructor") );
        bp::scope GridFF_scope( GridFF_exposer );
        GridFF_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "Construct a grid forcefield with a specified name") );
        GridFF_exposer.def( bp::init< SireMM::GridFF const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::GridFF::addFixedAtoms
        
            typedef void ( ::SireMM::GridFF::*addFixedAtoms_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            addFixedAtoms_function_type addFixedAtoms_function_value( &::SireMM::GridFF::addFixedAtoms );
            
            GridFF_exposer.def( 
                "addFixedAtoms"
                , addFixedAtoms_function_value
                , ( bp::arg("fixed_atoms"), bp::arg("map")=SireBase::PropertyMap() )
                , "Add fixed atoms to the grid. These are atoms that will never change\nposition or charge during the simulation, and that you wish to be\nincluded in the energy expression. The atoms can be placed here, and\nthen do not need to be added to the simulation System. This is useful\nif you are simulating a small cutout of the system and do not want to\nhave all of the atoms loaded into the system during the simulation" );
        
        }
        { //::SireMM::GridFF::addFixedAtoms
        
            typedef void ( ::SireMM::GridFF::*addFixedAtoms_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) ;
            addFixedAtoms_function_type addFixedAtoms_function_value( &::SireMM::GridFF::addFixedAtoms );
            
            GridFF_exposer.def( 
                "addFixedAtoms"
                , addFixedAtoms_function_value
                , ( bp::arg("fixed_atoms"), bp::arg("map")=SireBase::PropertyMap() )
                , "Add fixed atoms to the grid. These are atoms that will never change\nposition or charge during the simulation, and that you wish to be\nincluded in the energy expression. The atoms can be placed here, and\nthen do not need to be added to the simulation System. This is useful\nif you are simulating a small cutout of the system and do not want to\nhave all of the atoms loaded into the system during the simulation" );
        
        }
        { //::SireMM::GridFF::addFixedAtoms
        
            typedef void ( ::SireMM::GridFF::*addFixedAtoms_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) ;
            addFixedAtoms_function_type addFixedAtoms_function_value( &::SireMM::GridFF::addFixedAtoms );
            
            GridFF_exposer.def( 
                "addFixedAtoms"
                , addFixedAtoms_function_value
                , ( bp::arg("fixed_atoms"), bp::arg("map")=SireBase::PropertyMap() )
                , "Add fixed atoms to the grid. These are atoms that will never change\nposition or charge during the simulation, and that you wish to be\nincluded in the energy expression. The atoms can be placed here, and\nthen do not need to be added to the simulation System. This is useful\nif you are simulating a small cutout of the system and do not want to\nhave all of the atoms loaded into the system during the simulation" );
        
        }
        { //::SireMM::GridFF::addFixedAtoms
        
            typedef void ( ::SireMM::GridFF::*addFixedAtoms_function_type)( ::SireMM::GridFF const & ) ;
            addFixedAtoms_function_type addFixedAtoms_function_value( &::SireMM::GridFF::addFixedAtoms );
            
            GridFF_exposer.def( 
                "addFixedAtoms"
                , addFixedAtoms_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Add fixed atoms to the grid. This will copy the fixed atoms from\nthe passed GridFF. This allows multiple grid forcefields to share the\nmemory cost of a shared set of fixed atoms" );
        
        }
        { //::SireMM::GridFF::buffer
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::GridFF::*buffer_function_type)(  ) const;
            buffer_function_type buffer_function_value( &::SireMM::GridFF::buffer );
            
            GridFF_exposer.def( 
                "buffer"
                , buffer_function_value
                , bp::release_gil_policy()
                , "Return the buffer size used when building grids" );
        
        }
        { //::SireMM::GridFF::coulombCutoff
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::GridFF::*coulombCutoff_function_type)(  ) const;
            coulombCutoff_function_type coulombCutoff_function_value( &::SireMM::GridFF::coulombCutoff );
            
            GridFF_exposer.def( 
                "coulombCutoff"
                , coulombCutoff_function_value
                , bp::release_gil_policy()
                , "Return the cutoff for the coulomb energy" );
        
        }
        { //::SireMM::GridFF::ljCutoff
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::GridFF::*ljCutoff_function_type)(  ) const;
            ljCutoff_function_type ljCutoff_function_value( &::SireMM::GridFF::ljCutoff );
            
            GridFF_exposer.def( 
                "ljCutoff"
                , ljCutoff_function_value
                , bp::release_gil_policy()
                , "Return the cutoff for the LJ energy" );
        
        }
        { //::SireMM::GridFF::mustNowRecalculateFromScratch
        
            typedef void ( ::SireMM::GridFF::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireMM::GridFF::mustNowRecalculateFromScratch );
            
            GridFF_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , bp::release_gil_policy()
                , "Ensure that the next energy evaluation is from scratch" );
        
        }
        GridFF_exposer.def( bp::self != bp::self );
        { //::SireMM::GridFF::operator=
        
            typedef ::SireMM::GridFF & ( ::SireMM::GridFF::*assign_function_type)( ::SireMM::GridFF const & ) ;
            assign_function_type assign_function_value( &::SireMM::GridFF::operator= );
            
            GridFF_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GridFF_exposer.def( bp::self == bp::self );
        { //::SireMM::GridFF::setBuffer
        
            typedef void ( ::SireMM::GridFF::*setBuffer_function_type)( ::SireUnits::Dimension::Length ) ;
            setBuffer_function_type setBuffer_function_value( &::SireMM::GridFF::setBuffer );
            
            GridFF_exposer.def( 
                "setBuffer"
                , setBuffer_function_value
                , ( bp::arg("buffer") )
                , bp::release_gil_policy()
                , "Set the buffer when building the grid. This adds a buffer space\naround the grid when it is built, to try to reduce the number of\ntimes it needs to be rebuilt" );
        
        }
        { //::SireMM::GridFF::setCoulombCutoff
        
            typedef void ( ::SireMM::GridFF::*setCoulombCutoff_function_type)( ::SireUnits::Dimension::Length ) ;
            setCoulombCutoff_function_type setCoulombCutoff_function_value( &::SireMM::GridFF::setCoulombCutoff );
            
            GridFF_exposer.def( 
                "setCoulombCutoff"
                , setCoulombCutoff_function_value
                , ( bp::arg("cutoff") )
                , bp::release_gil_policy()
                , "Set the cutoff for the coulomb energy - this can be greater\nthan the box size as multiple periodic images can be used" );
        
        }
        { //::SireMM::GridFF::setGridSpacing
        
            typedef void ( ::SireMM::GridFF::*setGridSpacing_function_type)( ::SireUnits::Dimension::Length ) ;
            setGridSpacing_function_type setGridSpacing_function_value( &::SireMM::GridFF::setGridSpacing );
            
            GridFF_exposer.def( 
                "setGridSpacing"
                , setGridSpacing_function_value
                , ( bp::arg("spacing") )
                , bp::release_gil_policy()
                , "Set the grid spacing (the distance between grid points). The\nsmaller the spacing the more memory is required to hold the grid,\nbut the more accurate the energy" );
        
        }
        { //::SireMM::GridFF::setLJCutoff
        
            typedef void ( ::SireMM::GridFF::*setLJCutoff_function_type)( ::SireUnits::Dimension::Length ) ;
            setLJCutoff_function_type setLJCutoff_function_value( &::SireMM::GridFF::setLJCutoff );
            
            GridFF_exposer.def( 
                "setLJCutoff"
                , setLJCutoff_function_value
                , ( bp::arg("cutoff") )
                , bp::release_gil_policy()
                , "Set the cutoff for the LJ energy - this can be greater than\nthe box size as multiple periodic images can be used" );
        
        }
        { //::SireMM::GridFF::setReactionFieldDielectric
        
            typedef bool ( ::SireMM::GridFF::*setReactionFieldDielectric_function_type)( double ) ;
            setReactionFieldDielectric_function_type setReactionFieldDielectric_function_value( &::SireMM::GridFF::setReactionFieldDielectric );
            
            GridFF_exposer.def( 
                "setReactionFieldDielectric"
                , setReactionFieldDielectric_function_value
                , ( bp::arg("dielectric") )
                , bp::release_gil_policy()
                , "Set the dielectric constant to use with the reaction field potential" );
        
        }
        { //::SireMM::GridFF::setShiftElectrostatics
        
            typedef bool ( ::SireMM::GridFF::*setShiftElectrostatics_function_type)( bool ) ;
            setShiftElectrostatics_function_type setShiftElectrostatics_function_value( &::SireMM::GridFF::setShiftElectrostatics );
            
            GridFF_exposer.def( 
                "setShiftElectrostatics"
                , setShiftElectrostatics_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Turn on or off use of the force shifted potential" );
        
        }
        { //::SireMM::GridFF::setUseReactionField
        
            typedef bool ( ::SireMM::GridFF::*setUseReactionField_function_type)( bool ) ;
            setUseReactionField_function_type setUseReactionField_function_value( &::SireMM::GridFF::setUseReactionField );
            
            GridFF_exposer.def( 
                "setUseReactionField"
                , setUseReactionField_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Turn on or off the use of the reaction field" );
        
        }
        { //::SireMM::GridFF::spacing
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::GridFF::*spacing_function_type)(  ) const;
            spacing_function_type spacing_function_value( &::SireMM::GridFF::spacing );
            
            GridFF_exposer.def( 
                "spacing"
                , spacing_function_value
                , bp::release_gil_policy()
                , "Return the spacing between grid points" );
        
        }
        { //::SireMM::GridFF::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::GridFF::typeName );
            
            GridFF_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::GridFF::what
        
            typedef char const * ( ::SireMM::GridFF::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::GridFF::what );
            
            GridFF_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        GridFF_exposer.staticmethod( "typeName" );
        GridFF_exposer.def( "__copy__", &__copy__<SireMM::GridFF>);
        GridFF_exposer.def( "__deepcopy__", &__copy__<SireMM::GridFF>);
        GridFF_exposer.def( "clone", &__copy__<SireMM::GridFF>);
        GridFF_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::GridFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GridFF_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::GridFF >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GridFF_exposer.def_pickle(sire_pickle_suite< ::SireMM::GridFF >());
        GridFF_exposer.def( "__str__", &__str__< ::SireMM::GridFF > );
        GridFF_exposer.def( "__repr__", &__str__< ::SireMM::GridFF > );
        GridFF_exposer.def( "__len__", &__len_count< ::SireMM::GridFF > );
    }

}
