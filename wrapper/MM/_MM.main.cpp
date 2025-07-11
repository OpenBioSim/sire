// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License


#include "boost/python.hpp"

#include "Helpers/clone_const_reference.hpp"

#include "AmberAngle.pypp.hpp"

#include "AmberBond.pypp.hpp"

#include "AmberDihPart.pypp.hpp"

#include "AmberDihedral.pypp.hpp"

#include "AmberNB14.pypp.hpp"

#include "AmberParams.pypp.hpp"

#include "Angle.pypp.hpp"

#include "AngleComponent.pypp.hpp"

#include "AngleParameterName.pypp.hpp"

#include "AngleRestraint.pypp.hpp"

#include "AngleRestraints.pypp.hpp"

#include "AngleSymbols.pypp.hpp"

#include "AtomFunction.pypp.hpp"

#include "AtomFunctions.pypp.hpp"

#include "AtomLJs.pypp.hpp"

#include "AtomPairs_CLJScaleFactor_.pypp.hpp"

#include "AtomPairs_CoulombScaleFactor_.pypp.hpp"

#include "AtomPairs_LJScaleFactor_.pypp.hpp"

#include "BendBendComponent.pypp.hpp"

#include "BendBendParameterName.pypp.hpp"

#include "BendBendSymbols.pypp.hpp"

#include "Bond.pypp.hpp"

#include "BondComponent.pypp.hpp"

#include "BondParameterName.pypp.hpp"

#include "BondRestraint.pypp.hpp"

#include "BondRestraints.pypp.hpp"

#include "BondSymbols.pypp.hpp"

#include "BoreschRestraint.pypp.hpp"

#include "BoreschRestraints.pypp.hpp"

#include "CHARMMSwitchingFunction.pypp.hpp"

#include "CLJ14Group.pypp.hpp"

#include "CLJAtom.pypp.hpp"

#include "CLJAtoms.pypp.hpp"

#include "CLJBox.pypp.hpp"

#include "CLJBoxDistance.pypp.hpp"

#include "CLJBoxIndex.pypp.hpp"

#include "CLJBoxes.pypp.hpp"

#include "CLJCalculator.pypp.hpp"

#include "CLJComponent.pypp.hpp"

#include "CLJCutoffFunction.pypp.hpp"

#include "CLJDelta.pypp.hpp"

#include "CLJExtractor.pypp.hpp"

#include "CLJFunction.pypp.hpp"

#include "CLJGrid.pypp.hpp"

#include "CLJGroup.pypp.hpp"

#include "CLJIntraFunction.pypp.hpp"

#include "CLJIntraRFFunction.pypp.hpp"

#include "CLJIntraShiftFunction.pypp.hpp"

#include "CLJNBPairs.pypp.hpp"

#include "CLJParameterNames.pypp.hpp"

#include "CLJParameterNames3D.pypp.hpp"

#include "CLJPotentialInterface_InterCLJPotential_.pypp.hpp"

#include "CLJPotentialInterface_InterSoftCLJPotential_.pypp.hpp"

#include "CLJPotentialInterface_IntraCLJPotential_.pypp.hpp"

#include "CLJPotentialInterface_IntraSoftCLJPotential_.pypp.hpp"

#include "CLJProbe.pypp.hpp"

#include "CLJRFFunction.pypp.hpp"

#include "CLJScaleFactor.pypp.hpp"

#include "CLJShiftFunction.pypp.hpp"

#include "CLJSoftFunction.pypp.hpp"

#include "CLJSoftIntraFunction.pypp.hpp"

#include "CLJSoftIntraRFFunction.pypp.hpp"

#include "CLJSoftIntraShiftFunction.pypp.hpp"

#include "CLJSoftRFFunction.pypp.hpp"

#include "CLJSoftShiftFunction.pypp.hpp"

#include "CLJWorkspace.pypp.hpp"

#include "CMAPFunction.pypp.hpp"

#include "CMAPFunctions.pypp.hpp"

#include "CMAPParameter.pypp.hpp"

#include "ChargeParameterName.pypp.hpp"

#include "ChargeParameterName3D.pypp.hpp"

#include "CoulombComponent.pypp.hpp"

#include "CoulombNBPairs.pypp.hpp"

#include "CoulombPotentialInterface_InterCoulombPotential_.pypp.hpp"

#include "CoulombPotentialInterface_IntraCoulombPotential_.pypp.hpp"

#include "CoulombProbe.pypp.hpp"

#include "CoulombScaleFactor.pypp.hpp"

#include "Dihedral.pypp.hpp"

#include "DihedralComponent.pypp.hpp"

#include "DihedralParameterName.pypp.hpp"

#include "DihedralRestraint.pypp.hpp"

#include "DihedralRestraints.pypp.hpp"

#include "DihedralSymbols.pypp.hpp"

#include "DistanceRestraint.pypp.hpp"

#include "DoubleDistanceRestraint.pypp.hpp"

#include "ExcludedPairs.pypp.hpp"

#include "FourAtomFunction.pypp.hpp"

#include "FourAtomFunctions.pypp.hpp"

#include "FourAtomPerturbation.pypp.hpp"

#include "GridFF.pypp.hpp"

#include "GridFF2.pypp.hpp"

#include "GromacsAngle.pypp.hpp"

#include "GromacsAtomType.pypp.hpp"

#include "GromacsBond.pypp.hpp"

#include "GromacsDihedral.pypp.hpp"

#include "GroupInternalParameters.pypp.hpp"

#include "HarmonicSwitchingFunction.pypp.hpp"

#include "Improper.pypp.hpp"

#include "ImproperComponent.pypp.hpp"

#include "ImproperParameterName.pypp.hpp"

#include "ImproperSymbols.pypp.hpp"

#include "InterCLJFF.pypp.hpp"

#include "InterCLJFFBase.pypp.hpp"

#include "InterCoulombFF.pypp.hpp"

#include "InterCoulombFFBase.pypp.hpp"

#include "InterFF.pypp.hpp"

#include "InterGroupCLJFF.pypp.hpp"

#include "InterGroupCLJFFBase.pypp.hpp"

#include "InterGroupCoulombFF.pypp.hpp"

#include "InterGroupCoulombFFBase.pypp.hpp"

#include "InterGroupFF.pypp.hpp"

#include "InterGroupLJFF.pypp.hpp"

#include "InterGroupLJFFBase.pypp.hpp"

#include "InterGroupSoftCLJFF.pypp.hpp"

#include "InterGroupSoftCLJFFBase.pypp.hpp"

#include "InterLJFF.pypp.hpp"

#include "InterLJFFBase.pypp.hpp"

#include "InterSoftCLJFF.pypp.hpp"

#include "InterSoftCLJFFBase.pypp.hpp"

#include "InternalComponent.pypp.hpp"

#include "InternalFF.pypp.hpp"

#include "InternalGroupFF.pypp.hpp"

#include "InternalParameterNames.pypp.hpp"

#include "InternalParameterNames3D.pypp.hpp"

#include "InternalParameters.pypp.hpp"

#include "InternalParameters3D.pypp.hpp"

#include "InternalPerturbation.pypp.hpp"

#include "InternalSymbols.pypp.hpp"

#include "InternalSymbolsBase.pypp.hpp"

#include "Intra14Component.pypp.hpp"

#include "Intra14CoulombComponent.pypp.hpp"

#include "Intra14LJComponent.pypp.hpp"

#include "IntraCLJFF.pypp.hpp"

#include "IntraCLJFFBase.pypp.hpp"

#include "IntraCoulombFF.pypp.hpp"

#include "IntraCoulombFFBase.pypp.hpp"

#include "IntraFF.pypp.hpp"

#include "IntraGroupCLJFF.pypp.hpp"

#include "IntraGroupCLJFFBase.pypp.hpp"

#include "IntraGroupCoulombFF.pypp.hpp"

#include "IntraGroupCoulombFFBase.pypp.hpp"

#include "IntraGroupFF.pypp.hpp"

#include "IntraGroupLJFF.pypp.hpp"

#include "IntraGroupLJFFBase.pypp.hpp"

#include "IntraGroupSoftCLJFF.pypp.hpp"

#include "IntraGroupSoftCLJFFBase.pypp.hpp"

#include "IntraLJFF.pypp.hpp"

#include "IntraLJFFBase.pypp.hpp"

#include "IntraSoftCLJFF.pypp.hpp"

#include "IntraSoftCLJFFBase.pypp.hpp"

#include "LJ1264Parameter.pypp.hpp"

#include "LJComponent.pypp.hpp"

#include "LJException.pypp.hpp"

#include "LJExceptionID.pypp.hpp"

#include "LJNBPairs.pypp.hpp"

#include "LJParameter.pypp.hpp"

#include "LJParameterName.pypp.hpp"

#include "LJParameterName3D.pypp.hpp"

#include "LJPerturbation.pypp.hpp"

#include "LJPotentialInterface_InterLJPotential_.pypp.hpp"

#include "LJPotentialInterface_IntraLJPotential_.pypp.hpp"

#include "LJProbe.pypp.hpp"

#include "LJScaleFactor.pypp.hpp"

#include "MMDetail.pypp.hpp"

#include "Mover_Angle_.pypp.hpp"

#include "Mover_Bond_.pypp.hpp"

#include "Mover_Dihedral_.pypp.hpp"

#include "Mover_Improper_.pypp.hpp"

#include "Mover_SelectorAngle_.pypp.hpp"

#include "Mover_SelectorBond_.pypp.hpp"

#include "Mover_SelectorDihedral_.pypp.hpp"

#include "Mover_SelectorImproper_.pypp.hpp"

#include "MultiCLJComponent.pypp.hpp"

#include "NoCutoff.pypp.hpp"

#include "NullCLJFunction.pypp.hpp"

#include "NullRestraint.pypp.hpp"

#include "PositionalRestraint.pypp.hpp"

#include "PositionalRestraints.pypp.hpp"

#include "Restraint.pypp.hpp"

#include "Restraint3D.pypp.hpp"

#include "RestraintComponent.pypp.hpp"

#include "RestraintFF.pypp.hpp"

#include "Restraints.pypp.hpp"

#include "ScaledCLJParameterNames3D.pypp.hpp"

#include "ScaledChargeParameterNames3D.pypp.hpp"

#include "ScaledLJParameterNames3D.pypp.hpp"

#include "SelectorAngle.pypp.hpp"

#include "SelectorBond.pypp.hpp"

#include "SelectorDihedral.pypp.hpp"

#include "SelectorImproper.pypp.hpp"

#include "SelectorMAngle.pypp.hpp"

#include "SelectorMBond.pypp.hpp"

#include "SelectorMDihedral.pypp.hpp"

#include "SelectorMImproper.pypp.hpp"

#include "SoftCLJComponent.pypp.hpp"

#include "SoftCLJPotentialInterface_InterSoftCLJPotential_.pypp.hpp"

#include "SoftCLJPotentialInterface_IntraSoftCLJPotential_.pypp.hpp"

#include "StretchBendComponent.pypp.hpp"

#include "StretchBendParameterName.pypp.hpp"

#include "StretchBendSymbols.pypp.hpp"

#include "StretchBendTorsionComponent.pypp.hpp"

#include "StretchBendTorsionParameterName.pypp.hpp"

#include "StretchBendTorsionSymbols.pypp.hpp"

#include "StretchStretchComponent.pypp.hpp"

#include "StretchStretchParameterName.pypp.hpp"

#include "StretchStretchSymbols.pypp.hpp"

#include "SwitchingFunction.pypp.hpp"

#include "TestFF.pypp.hpp"

#include "ThreeAtomFunction.pypp.hpp"

#include "ThreeAtomFunctions.pypp.hpp"

#include "ThreeAtomPerturbation.pypp.hpp"

#include "TripleDistanceRestraint.pypp.hpp"

#include "TwoAtomFunction.pypp.hpp"

#include "TwoAtomFunctions.pypp.hpp"

#include "TwoAtomPerturbation.pypp.hpp"

#include "UreyBradleyComponent.pypp.hpp"

#include "UreyBradleyParameterName.pypp.hpp"

namespace bp = boost::python;

#include "SireMM_containers.h"

#include "SireMM_registrars.h"

#include "SireMM_properties.h"

BOOST_PYTHON_MODULE(_MM){
    register_SireMM_objects();

    register_SireMM_containers();

    register_CLJPotentialInterface_InterCLJPotential__class();

    register_InterGroupCLJFFBase_class();

    register_InterGroupCLJFF_class();

    register_CoulombPotentialInterface_InterCoulombPotential__class();

    register_InterGroupCoulombFFBase_class();

    register_InterGroupCoulombFF_class();

    register_LJPotentialInterface_InterLJPotential__class();

    register_InterGroupLJFFBase_class();

    register_InterGroupLJFF_class();

    register_CLJPotentialInterface_InterSoftCLJPotential__class();

    register_SoftCLJPotentialInterface_InterSoftCLJPotential__class();

    register_InterGroupSoftCLJFFBase_class();

    register_InterGroupSoftCLJFF_class();

    register_InterCLJFFBase_class();

    register_InterCLJFF_class();

    register_InterCoulombFFBase_class();

    register_InterCoulombFF_class();

    register_InterLJFFBase_class();

    register_InterLJFF_class();

    register_InterSoftCLJFFBase_class();

    register_InterSoftCLJFF_class();

    register_CLJPotentialInterface_IntraCLJPotential__class();

    register_IntraGroupCLJFFBase_class();

    register_IntraGroupCLJFF_class();

    register_CoulombPotentialInterface_IntraCoulombPotential__class();

    register_IntraGroupCoulombFFBase_class();

    register_IntraGroupCoulombFF_class();

    register_LJPotentialInterface_IntraLJPotential__class();

    register_IntraGroupLJFFBase_class();

    register_IntraGroupLJFF_class();

    register_CLJPotentialInterface_IntraSoftCLJPotential__class();

    register_SoftCLJPotentialInterface_IntraSoftCLJPotential__class();

    register_IntraGroupSoftCLJFFBase_class();

    register_IntraGroupSoftCLJFF_class();

    register_IntraCLJFFBase_class();

    register_IntraCLJFF_class();

    register_IntraCoulombFFBase_class();

    register_IntraCoulombFF_class();

    register_IntraLJFFBase_class();

    register_IntraLJFF_class();

    register_IntraSoftCLJFFBase_class();

    register_IntraSoftCLJFF_class();

    register_AmberAngle_class();

    register_AmberBond_class();

    register_AmberDihPart_class();

    register_AmberDihedral_class();

    register_AmberNB14_class();

    register_AmberParams_class();

    register_Angle_class();

    register_AngleComponent_class();

    register_AngleParameterName_class();

    register_AngleRestraint_class();

    register_Restraints_class();

    register_AngleRestraints_class();

    register_InternalSymbolsBase_class();

    register_AngleSymbols_class();

    register_AtomFunction_class();

    register_AtomFunctions_class();

    register_CoulombScaleFactor_class();

    register_LJScaleFactor_class();

    register_CLJScaleFactor_class();

    register_AtomPairs_CLJScaleFactor__class();

    register_AtomPairs_CoulombScaleFactor__class();

    register_AtomPairs_LJScaleFactor__class();

    register_BendBendComponent_class();

    register_BendBendParameterName_class();

    register_BendBendSymbols_class();

    register_Bond_class();

    register_BondComponent_class();

    register_BondParameterName_class();

    register_BondRestraint_class();

    register_BondRestraints_class();

    register_BondSymbols_class();

    register_BoreschRestraint_class();

    register_BoreschRestraints_class();

    register_SwitchingFunction_class();

    register_CHARMMSwitchingFunction_class();

    register_CLJFunction_class();

    register_CLJ14Group_class();

    register_CLJAtom_class();

    register_CLJAtoms_class();

    register_CLJBox_class();

    register_CLJBoxDistance_class();

    register_CLJBoxIndex_class();

    register_CLJBoxes_class();

    register_CLJCalculator_class();

    register_CLJComponent_class();

    register_CLJCutoffFunction_class();

    register_CLJDelta_class();

    register_CLJExtractor_class();

    register_CLJGrid_class();

    register_CLJGroup_class();

    register_CLJIntraFunction_class();

    register_CLJIntraRFFunction_class();

    register_CLJIntraShiftFunction_class();

    register_CLJNBPairs_class();

    register_ChargeParameterName_class();

    register_LJParameterName_class();

    register_CLJParameterNames_class();

    register_CLJParameterNames3D_class();

    register_CLJProbe_class();

    register_CLJRFFunction_class();

    register_CLJShiftFunction_class();

    register_CLJSoftFunction_class();

    register_CLJSoftIntraFunction_class();

    register_CLJSoftIntraRFFunction_class();

    register_CLJSoftIntraShiftFunction_class();

    register_CLJSoftRFFunction_class();

    register_CLJSoftShiftFunction_class();

    register_CLJWorkspace_class();

    register_CMAPFunction_class();

    register_CMAPFunctions_class();

    register_CMAPParameter_class();

    register_ChargeParameterName3D_class();

    register_CoulombComponent_class();

    register_CoulombNBPairs_class();

    register_CoulombProbe_class();

    register_Dihedral_class();

    register_DihedralComponent_class();

    register_DihedralParameterName_class();

    register_DihedralRestraint_class();

    register_DihedralRestraints_class();

    register_DihedralSymbols_class();

    register_Restraint_class();

    register_Restraint3D_class();

    register_DistanceRestraint_class();

    register_DoubleDistanceRestraint_class();

    register_ExcludedPairs_class();

    register_FourAtomFunction_class();

    register_FourAtomFunctions_class();

    register_InternalPerturbation_class();

    register_FourAtomPerturbation_class();

    register_GridFF_class();

    register_GridFF2_class();

    register_GromacsAngle_class();

    register_GromacsAtomType_class();

    register_GromacsBond_class();

    register_GromacsDihedral_class();

    register_GroupInternalParameters_class();

    register_HarmonicSwitchingFunction_class();

    register_Improper_class();

    register_ImproperComponent_class();

    register_ImproperParameterName_class();

    register_ImproperSymbols_class();

    register_InterFF_class();

    register_InterGroupFF_class();

    register_InternalComponent_class();

    register_InternalFF_class();

    register_InternalGroupFF_class();

    register_StretchBendParameterName_class();

    register_StretchBendTorsionParameterName_class();

    register_StretchStretchParameterName_class();

    register_UreyBradleyParameterName_class();

    register_InternalParameterNames_class();

    register_InternalParameterNames3D_class();

    register_InternalParameters_class();

    register_InternalParameters3D_class();

    register_InternalSymbols_class();

    register_Intra14Component_class();

    register_Intra14CoulombComponent_class();

    register_Intra14LJComponent_class();

    register_IntraFF_class();

    register_IntraGroupFF_class();

    register_LJ1264Parameter_class();

    register_LJComponent_class();

    register_LJException_class();

    register_LJExceptionID_class();

    register_LJNBPairs_class();

    register_LJParameter_class();

    register_LJParameterName3D_class();

    register_LJPerturbation_class();

    register_LJProbe_class();

    register_MMDetail_class();

    register_MultiCLJComponent_class();

    register_NoCutoff_class();

    register_NullCLJFunction_class();

    register_NullRestraint_class();

    register_PositionalRestraint_class();

    register_PositionalRestraints_class();

    register_RestraintComponent_class();

    register_RestraintFF_class();

    register_ScaledCLJParameterNames3D_class();

    register_ScaledChargeParameterNames3D_class();

    register_ScaledLJParameterNames3D_class();

    register_SelectorAngle_class();

    register_SelectorBond_class();

    register_SelectorDihedral_class();

    register_SelectorImproper_class();

    register_SelectorMAngle_class();

    register_SelectorMBond_class();

    register_SelectorMDihedral_class();

    register_SelectorMImproper_class();

    register_SoftCLJComponent_class();

    register_StretchBendComponent_class();

    register_StretchBendSymbols_class();

    register_StretchBendTorsionComponent_class();

    register_StretchBendTorsionSymbols_class();

    register_StretchStretchComponent_class();

    register_StretchStretchSymbols_class();

    register_TestFF_class();

    register_ThreeAtomFunction_class();

    register_ThreeAtomFunctions_class();

    register_ThreeAtomPerturbation_class();

    register_TripleDistanceRestraint_class();

    register_TwoAtomFunction_class();

    register_TwoAtomFunctions_class();

    register_TwoAtomPerturbation_class();

    register_UreyBradleyComponent_class();

    register_AtomLJs_class();

    register_Mover_Angle__class();

    register_Mover_Bond__class();

    register_Mover_Dihedral__class();

    register_Mover_Improper__class();

    register_Mover_SelectorAngle__class();

    register_Mover_SelectorBond__class();

    register_Mover_SelectorDihedral__class();

    register_SireMM_properties();

    register_Mover_SelectorImproper__class();
}

