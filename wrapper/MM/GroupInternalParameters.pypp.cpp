// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "GroupInternalParameters.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireFF/errors.h"

#include "SireMol/cgidx.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireVol/coordgroup.h"

#include "internalparameters.h"

#include "sireglobal.h"

#include "tostring.h"

#include <algorithm>

#include "internalparameters.h"

SireMM::GroupInternalParameters __copy__(const SireMM::GroupInternalParameters &other){ return SireMM::GroupInternalParameters(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMM::GroupInternalParameters&){ return "SireMM::GroupInternalParameters";}

#include "Helpers/release_gil_policy.hpp"

void register_GroupInternalParameters_class(){

    { //::SireMM::GroupInternalParameters
        typedef bp::class_< SireMM::GroupInternalParameters > GroupInternalParameters_exposer_t;
        GroupInternalParameters_exposer_t GroupInternalParameters_exposer = GroupInternalParameters_exposer_t( "GroupInternalParameters", "This class holds all of the internal parameters for one group\ncombination within a molecule\n\nThere are several types of internal parameters, defined\nby the type of internal used to provide the coordinates,\nand the quantities calculated from those coordinates that\ncan be used in the function\n\nBond            : Input Bond - function uses interatomic distance (1-2), r\nAngle           : Input Angle - function uses angle (1-2-3), theta\nDihedral        : Input Dihedral - function uses torsion (1-2-3-4), phi\n\nImproper        : Input Improper - function uses either torsion angle\n(1-3-4-2), phi, or out of plane angle, theta\n\nUrey-Bradley    : Input Angle - function uses distance (1-3), r\n\nStretch-Stretch : Input Angle - function uses distances (1-2), r12, and\n(3-2), r32\n\nStretch-Bend    : Input Angle - function uses distances angle (1-2-3), theta,\nand distances (1-2), r12, and (3-2) r32\n\nBend-Bend       : Input Improper - function uses angles (1-2-3), (3-2-4), (4-2-1),\ntheta123, theta324, theta421\n\nStretch-Bend    : Input Dihedral - function uses torsion (1-2-3-4), phi,\n-Torsion                         distances (1-2), (2-3), (3-4), (1-4)\nr12, r23, r34, r14 and angles\n(1-2-3) and (2-3-4), theta123, theta234\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope GroupInternalParameters_scope( GroupInternalParameters_exposer );
        GroupInternalParameters_exposer.def( bp::init< SireMM::GroupInternalParameters const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::GroupInternalParameters::angleForces
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*angleForces_function_type)(  ) const;
            angleForces_function_type angleForces_function_value( &::SireMM::GroupInternalParameters::angleForces );
            
            GroupInternalParameters_exposer.def( 
                "angleForces"
                , angleForces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the angle force ( -dEdtheta )" );
        
        }
        { //::SireMM::GroupInternalParameters::anglePotential
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*anglePotential_function_type)(  ) const;
            anglePotential_function_type anglePotential_function_value( &::SireMM::GroupInternalParameters::anglePotential );
            
            GroupInternalParameters_exposer.def( 
                "anglePotential"
                , anglePotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the angle potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::bendBendPotential
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*bendBendPotential_function_type)(  ) const;
            bendBendPotential_function_type bendBendPotential_function_value( &::SireMM::GroupInternalParameters::bendBendPotential );
            
            GroupInternalParameters_exposer.def( 
                "bendBendPotential"
                , bendBendPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the bend-bend potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::bendBend_Theta012_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*bendBend_Theta012_Forces_function_type)(  ) const;
            bendBend_Theta012_Forces_function_type bendBend_Theta012_Forces_function_value( &::SireMM::GroupInternalParameters::bendBend_Theta012_Forces );
            
            GroupInternalParameters_exposer.def( 
                "bendBend_Theta012_Forces"
                , bendBend_Theta012_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return bend-bend force ( -dEdtheta_012 )" );
        
        }
        { //::SireMM::GroupInternalParameters::bendBend_Theta213_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*bendBend_Theta213_Forces_function_type)(  ) const;
            bendBend_Theta213_Forces_function_type bendBend_Theta213_Forces_function_value( &::SireMM::GroupInternalParameters::bendBend_Theta213_Forces );
            
            GroupInternalParameters_exposer.def( 
                "bendBend_Theta213_Forces"
                , bendBend_Theta213_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return bend-bend force ( -dEdtheta_012 )" );
        
        }
        { //::SireMM::GroupInternalParameters::bendBend_Theta310_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*bendBend_Theta310_Forces_function_type)(  ) const;
            bendBend_Theta310_Forces_function_type bendBend_Theta310_Forces_function_value( &::SireMM::GroupInternalParameters::bendBend_Theta310_Forces );
            
            GroupInternalParameters_exposer.def( 
                "bendBend_Theta310_Forces"
                , bendBend_Theta310_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return bend-bend force ( -dEdtheta_012 )" );
        
        }
        { //::SireMM::GroupInternalParameters::bondForces
        
            typedef ::QVector< SireMM::TwoAtomFunction > const & ( ::SireMM::GroupInternalParameters::*bondForces_function_type)(  ) const;
            bondForces_function_type bondForces_function_value( &::SireMM::GroupInternalParameters::bondForces );
            
            GroupInternalParameters_exposer.def( 
                "bondForces"
                , bondForces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the bond force ( -dEdr )" );
        
        }
        { //::SireMM::GroupInternalParameters::bondPotential
        
            typedef ::QVector< SireMM::TwoAtomFunction > const & ( ::SireMM::GroupInternalParameters::*bondPotential_function_type)(  ) const;
            bondPotential_function_type bondPotential_function_value( &::SireMM::GroupInternalParameters::bondPotential );
            
            GroupInternalParameters_exposer.def( 
                "bondPotential"
                , bondPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return all of the bond potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::cgIdx0
        
            typedef ::SireMol::CGIdx ( ::SireMM::GroupInternalParameters::*cgIdx0_function_type)(  ) const;
            cgIdx0_function_type cgIdx0_function_value( &::SireMM::GroupInternalParameters::cgIdx0 );
            
            GroupInternalParameters_exposer.def( 
                "cgIdx0"
                , cgIdx0_function_value
                , bp::release_gil_policy()
                , "Return the index of the first group" );
        
        }
        { //::SireMM::GroupInternalParameters::cgIdx1
        
            typedef ::SireMol::CGIdx ( ::SireMM::GroupInternalParameters::*cgIdx1_function_type)(  ) const;
            cgIdx1_function_type cgIdx1_function_value( &::SireMM::GroupInternalParameters::cgIdx1 );
            
            GroupInternalParameters_exposer.def( 
                "cgIdx1"
                , cgIdx1_function_value
                , bp::release_gil_policy()
                , "Return the index of the first group" );
        
        }
        { //::SireMM::GroupInternalParameters::cgIdx2
        
            typedef ::SireMol::CGIdx ( ::SireMM::GroupInternalParameters::*cgIdx2_function_type)(  ) const;
            cgIdx2_function_type cgIdx2_function_value( &::SireMM::GroupInternalParameters::cgIdx2 );
            
            GroupInternalParameters_exposer.def( 
                "cgIdx2"
                , cgIdx2_function_value
                , bp::release_gil_policy()
                , "Return the index of the first group" );
        
        }
        { //::SireMM::GroupInternalParameters::cgIdx3
        
            typedef ::SireMol::CGIdx ( ::SireMM::GroupInternalParameters::*cgIdx3_function_type)(  ) const;
            cgIdx3_function_type cgIdx3_function_value( &::SireMM::GroupInternalParameters::cgIdx3 );
            
            GroupInternalParameters_exposer.def( 
                "cgIdx3"
                , cgIdx3_function_value
                , bp::release_gil_policy()
                , "Return the index of the first group" );
        
        }
        { //::SireMM::GroupInternalParameters::dihedralForces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*dihedralForces_function_type)(  ) const;
            dihedralForces_function_type dihedralForces_function_value( &::SireMM::GroupInternalParameters::dihedralForces );
            
            GroupInternalParameters_exposer.def( 
                "dihedralForces"
                , dihedralForces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the dihedral force ( -dEdphi )" );
        
        }
        { //::SireMM::GroupInternalParameters::dihedralPotential
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*dihedralPotential_function_type)(  ) const;
            dihedralPotential_function_type dihedralPotential_function_value( &::SireMM::GroupInternalParameters::dihedralPotential );
            
            GroupInternalParameters_exposer.def( 
                "dihedralPotential"
                , dihedralPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the dihedral potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::hasCrossTerms
        
            typedef bool ( ::SireMM::GroupInternalParameters::*hasCrossTerms_function_type)(  ) const;
            hasCrossTerms_function_type hasCrossTerms_function_value( &::SireMM::GroupInternalParameters::hasCrossTerms );
            
            GroupInternalParameters_exposer.def( 
                "hasCrossTerms"
                , hasCrossTerms_function_value
                , bp::release_gil_policy()
                , "Return whether or not this has any cross terms\n(stretch-stretch, stretch-bend, bend-bend, stretch-bend-torsion)" );
        
        }
        { //::SireMM::GroupInternalParameters::hasNonPhysicalParameters
        
            typedef bool ( ::SireMM::GroupInternalParameters::*hasNonPhysicalParameters_function_type)(  ) const;
            hasNonPhysicalParameters_function_type hasNonPhysicalParameters_function_value( &::SireMM::GroupInternalParameters::hasNonPhysicalParameters );
            
            GroupInternalParameters_exposer.def( 
                "hasNonPhysicalParameters"
                , hasNonPhysicalParameters_function_value
                , bp::release_gil_policy()
                , "Return whether or not this has any non-physical parameters\n(Urey-Bradley or improper terms)" );
        
        }
        { //::SireMM::GroupInternalParameters::hasPhysicalParameters
        
            typedef bool ( ::SireMM::GroupInternalParameters::*hasPhysicalParameters_function_type)(  ) const;
            hasPhysicalParameters_function_type hasPhysicalParameters_function_value( &::SireMM::GroupInternalParameters::hasPhysicalParameters );
            
            GroupInternalParameters_exposer.def( 
                "hasPhysicalParameters"
                , hasPhysicalParameters_function_value
                , bp::release_gil_policy()
                , "Return whether or not this has any physical parameters\n(bond, angle or dihedral)" );
        
        }
        { //::SireMM::GroupInternalParameters::improperPotential
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*improperPotential_function_type)(  ) const;
            improperPotential_function_type improperPotential_function_value( &::SireMM::GroupInternalParameters::improperPotential );
            
            GroupInternalParameters_exposer.def( 
                "improperPotential"
                , improperPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the improper potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::improper_Phi_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*improper_Phi_Forces_function_type)(  ) const;
            improper_Phi_Forces_function_type improper_Phi_Forces_function_value( &::SireMM::GroupInternalParameters::improper_Phi_Forces );
            
            GroupInternalParameters_exposer.def( 
                "improper_Phi_Forces"
                , improper_Phi_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the improper force ( -dEdphi )" );
        
        }
        { //::SireMM::GroupInternalParameters::improper_Theta_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*improper_Theta_Forces_function_type)(  ) const;
            improper_Theta_Forces_function_type improper_Theta_Forces_function_value( &::SireMM::GroupInternalParameters::improper_Theta_Forces );
            
            GroupInternalParameters_exposer.def( 
                "improper_Theta_Forces"
                , improper_Theta_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the improper force ( -dEdtheta )" );
        
        }
        { //::SireMM::GroupInternalParameters::isDoubleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isDoubleCutGroup_function_type)(  ) const;
            isDoubleCutGroup_function_type isDoubleCutGroup_function_value( &::SireMM::GroupInternalParameters::isDoubleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isDoubleCutGroup"
                , isDoubleCutGroup_function_value
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from only\ntwo CutGroups" );
        
        }
        { //::SireMM::GroupInternalParameters::isDoubleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isDoubleCutGroup_function_type)( ::SireMol::CGIdx,::SireMol::CGIdx ) const;
            isDoubleCutGroup_function_type isDoubleCutGroup_function_value( &::SireMM::GroupInternalParameters::isDoubleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isDoubleCutGroup"
                , isDoubleCutGroup_function_value
                , ( bp::arg("cgidx0"), bp::arg("cgidx1") )
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from\nonly two CutGroups, with indicies cgidx0 and cgidx1" );
        
        }
        { //::SireMM::GroupInternalParameters::isEmpty
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMM::GroupInternalParameters::isEmpty );
            
            GroupInternalParameters_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether this group is empty (contains no parameters)" );
        
        }
        { //::SireMM::GroupInternalParameters::isQuadrupleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isQuadrupleCutGroup_function_type)(  ) const;
            isQuadrupleCutGroup_function_type isQuadrupleCutGroup_function_value( &::SireMM::GroupInternalParameters::isQuadrupleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isQuadrupleCutGroup"
                , isQuadrupleCutGroup_function_value
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from\nfour CutGroups" );
        
        }
        { //::SireMM::GroupInternalParameters::isQuadrupleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isQuadrupleCutGroup_function_type)( ::SireMol::CGIdx,::SireMol::CGIdx,::SireMol::CGIdx,::SireMol::CGIdx ) const;
            isQuadrupleCutGroup_function_type isQuadrupleCutGroup_function_value( &::SireMM::GroupInternalParameters::isQuadrupleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isQuadrupleCutGroup"
                , isQuadrupleCutGroup_function_value
                , ( bp::arg("cgidx0"), bp::arg("cgidx1"), bp::arg("cgidx2"), bp::arg("cgidx3") )
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from\nfour CutGroups, with indicies cgidx0, cgidx1, cgidx2 and cgidx3" );
        
        }
        { //::SireMM::GroupInternalParameters::isSingleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isSingleCutGroup_function_type)(  ) const;
            isSingleCutGroup_function_type isSingleCutGroup_function_value( &::SireMM::GroupInternalParameters::isSingleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isSingleCutGroup"
                , isSingleCutGroup_function_value
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from only\na single CutGroup" );
        
        }
        { //::SireMM::GroupInternalParameters::isSingleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isSingleCutGroup_function_type)( ::SireMol::CGIdx ) const;
            isSingleCutGroup_function_type isSingleCutGroup_function_value( &::SireMM::GroupInternalParameters::isSingleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isSingleCutGroup"
                , isSingleCutGroup_function_value
                , ( bp::arg("cgidx0") )
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from\nonly a single CutGroup, with index cgidx0" );
        
        }
        { //::SireMM::GroupInternalParameters::isTripleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isTripleCutGroup_function_type)(  ) const;
            isTripleCutGroup_function_type isTripleCutGroup_function_value( &::SireMM::GroupInternalParameters::isTripleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isTripleCutGroup"
                , isTripleCutGroup_function_value
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from only\nthree CutGroups" );
        
        }
        { //::SireMM::GroupInternalParameters::isTripleCutGroup
        
            typedef bool ( ::SireMM::GroupInternalParameters::*isTripleCutGroup_function_type)( ::SireMol::CGIdx,::SireMol::CGIdx,::SireMol::CGIdx ) const;
            isTripleCutGroup_function_type isTripleCutGroup_function_value( &::SireMM::GroupInternalParameters::isTripleCutGroup );
            
            GroupInternalParameters_exposer.def( 
                "isTripleCutGroup"
                , isTripleCutGroup_function_value
                , ( bp::arg("cgidx0"), bp::arg("cgidx1"), bp::arg("cgidx2") )
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters from\nthree CutGroups, with indicies cgidx0, cgidx1 and cgidx2" );
        
        }
        GroupInternalParameters_exposer.def( bp::self != bp::self );
        { //::SireMM::GroupInternalParameters::operator=
        
            typedef ::SireMM::GroupInternalParameters & ( ::SireMM::GroupInternalParameters::*assign_function_type)( ::SireMM::GroupInternalParameters const & ) ;
            assign_function_type assign_function_value( &::SireMM::GroupInternalParameters::operator= );
            
            GroupInternalParameters_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GroupInternalParameters_exposer.def( bp::self == bp::self );
        { //::SireMM::GroupInternalParameters::refersTo
        
            typedef bool ( ::SireMM::GroupInternalParameters::*refersTo_function_type)( ::SireMol::CGIdx ) const;
            refersTo_function_type refersTo_function_value( &::SireMM::GroupInternalParameters::refersTo );
            
            GroupInternalParameters_exposer.def( 
                "refersTo"
                , refersTo_function_value
                , ( bp::arg("cgidx") )
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters\nthat involve atoms in the CutGroup at index cgidx" );
        
        }
        { //::SireMM::GroupInternalParameters::refersTo
        
            typedef bool ( ::SireMM::GroupInternalParameters::*refersTo_function_type)( ::QSet< SireMol::CGIdx > const & ) const;
            refersTo_function_type refersTo_function_value( &::SireMM::GroupInternalParameters::refersTo );
            
            GroupInternalParameters_exposer.def( 
                "refersTo"
                , refersTo_function_value
                , ( bp::arg("cgidxs") )
                , bp::release_gil_policy()
                , "Return whether or not this group contains parameters\nthat involve atoms in any of the CutGroups whose indicies are\nin cgidxs" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendPotential
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendPotential_function_type)(  ) const;
            stretchBendPotential_function_type stretchBendPotential_function_value( &::SireMM::GroupInternalParameters::stretchBendPotential );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendPotential"
                , stretchBendPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsionPotential
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsionPotential_function_type)(  ) const;
            stretchBendTorsionPotential_function_type stretchBendTorsionPotential_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsionPotential );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsionPotential"
                , stretchBendTorsionPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsion_Phi_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsion_Phi_Forces_function_type)(  ) const;
            stretchBendTorsion_Phi_Forces_function_type stretchBendTorsion_Phi_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsion_Phi_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsion_Phi_Forces"
                , stretchBendTorsion_Phi_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion force ( -dEdphi )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsion_R01_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsion_R01_Forces_function_type)(  ) const;
            stretchBendTorsion_R01_Forces_function_type stretchBendTorsion_R01_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsion_R01_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsion_R01_Forces"
                , stretchBendTorsion_R01_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion force ( -dEdr_01 )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsion_R03_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsion_R03_Forces_function_type)(  ) const;
            stretchBendTorsion_R03_Forces_function_type stretchBendTorsion_R03_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsion_R03_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsion_R03_Forces"
                , stretchBendTorsion_R03_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion force ( -dEdr_03 )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsion_R12_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsion_R12_Forces_function_type)(  ) const;
            stretchBendTorsion_R12_Forces_function_type stretchBendTorsion_R12_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsion_R12_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsion_R12_Forces"
                , stretchBendTorsion_R12_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion force ( -dEdr_12 )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsion_R32_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsion_R32_Forces_function_type)(  ) const;
            stretchBendTorsion_R32_Forces_function_type stretchBendTorsion_R32_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsion_R32_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsion_R32_Forces"
                , stretchBendTorsion_R32_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion force ( -dEdr_32 )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsion_Theta012_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsion_Theta012_Forces_function_type)(  ) const;
            stretchBendTorsion_Theta012_Forces_function_type stretchBendTorsion_Theta012_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsion_Theta012_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsion_Theta012_Forces"
                , stretchBendTorsion_Theta012_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion force ( -dEdtheta_012 )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBendTorsion_Theta321_Forces
        
            typedef ::QVector< SireMM::FourAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBendTorsion_Theta321_Forces_function_type)(  ) const;
            stretchBendTorsion_Theta321_Forces_function_type stretchBendTorsion_Theta321_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBendTorsion_Theta321_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBendTorsion_Theta321_Forces"
                , stretchBendTorsion_Theta321_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend-torsion force ( -dEdtheta_321 )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBend_R01_Forces
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBend_R01_Forces_function_type)(  ) const;
            stretchBend_R01_Forces_function_type stretchBend_R01_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBend_R01_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBend_R01_Forces"
                , stretchBend_R01_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend force ( -dEdtheta )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBend_R21_Forces
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBend_R21_Forces_function_type)(  ) const;
            stretchBend_R21_Forces_function_type stretchBend_R21_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBend_R21_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBend_R21_Forces"
                , stretchBend_R21_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend force ( -dEdtheta )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchBend_Theta_Forces
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchBend_Theta_Forces_function_type)(  ) const;
            stretchBend_Theta_Forces_function_type stretchBend_Theta_Forces_function_value( &::SireMM::GroupInternalParameters::stretchBend_Theta_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchBend_Theta_Forces"
                , stretchBend_Theta_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-bend force ( -dEdtheta )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchStretchPotential
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchStretchPotential_function_type)(  ) const;
            stretchStretchPotential_function_type stretchStretchPotential_function_value( &::SireMM::GroupInternalParameters::stretchStretchPotential );
            
            GroupInternalParameters_exposer.def( 
                "stretchStretchPotential"
                , stretchStretchPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-stretch potentials for this group" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchStretch_R01_Forces
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchStretch_R01_Forces_function_type)(  ) const;
            stretchStretch_R01_Forces_function_type stretchStretch_R01_Forces_function_value( &::SireMM::GroupInternalParameters::stretchStretch_R01_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchStretch_R01_Forces"
                , stretchStretch_R01_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-stretch force ( -dEdr_01 )" );
        
        }
        { //::SireMM::GroupInternalParameters::stretchStretch_R21_Forces
        
            typedef ::QVector< SireMM::ThreeAtomFunction > const & ( ::SireMM::GroupInternalParameters::*stretchStretch_R21_Forces_function_type)(  ) const;
            stretchStretch_R21_Forces_function_type stretchStretch_R21_Forces_function_value( &::SireMM::GroupInternalParameters::stretchStretch_R21_Forces );
            
            GroupInternalParameters_exposer.def( 
                "stretchStretch_R21_Forces"
                , stretchStretch_R21_Forces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the stretch-stretch force ( -dEdr_01 )" );
        
        }
        { //::SireMM::GroupInternalParameters::ureyBradleyForces
        
            typedef ::QVector< SireMM::TwoAtomFunction > const & ( ::SireMM::GroupInternalParameters::*ureyBradleyForces_function_type)(  ) const;
            ureyBradleyForces_function_type ureyBradleyForces_function_value( &::SireMM::GroupInternalParameters::ureyBradleyForces );
            
            GroupInternalParameters_exposer.def( 
                "ureyBradleyForces"
                , ureyBradleyForces_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the Urey-Bradley force ( -dEdr )" );
        
        }
        { //::SireMM::GroupInternalParameters::ureyBradleyPotential
        
            typedef ::QVector< SireMM::TwoAtomFunction > const & ( ::SireMM::GroupInternalParameters::*ureyBradleyPotential_function_type)(  ) const;
            ureyBradleyPotential_function_type ureyBradleyPotential_function_value( &::SireMM::GroupInternalParameters::ureyBradleyPotential );
            
            GroupInternalParameters_exposer.def( 
                "ureyBradleyPotential"
                , ureyBradleyPotential_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the Urey-Bradley potentials for this group" );
        
        }
        GroupInternalParameters_exposer.def( "__copy__", &__copy__<SireMM::GroupInternalParameters>);
        GroupInternalParameters_exposer.def( "__deepcopy__", &__copy__<SireMM::GroupInternalParameters>);
        GroupInternalParameters_exposer.def( "clone", &__copy__<SireMM::GroupInternalParameters>);
        GroupInternalParameters_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::GroupInternalParameters >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GroupInternalParameters_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::GroupInternalParameters >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GroupInternalParameters_exposer.def_pickle(sire_pickle_suite< ::SireMM::GroupInternalParameters >());
        GroupInternalParameters_exposer.def( "__str__", &pvt_get_name);
        GroupInternalParameters_exposer.def( "__repr__", &pvt_get_name);
    }

}
