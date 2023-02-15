//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireMM_registrars.h"

#include "amberparams.h"
#include "angle.h"
#include "anglerestraint.h"
#include "atomljs.h"
#include "bond.h"
#include "clj14group.h"
#include "cljatoms.h"
#include "cljboxes.h"
#include "cljcalculator.h"
#include "cljcomponent.h"
#include "cljdelta.h"
#include "cljextractor.h"
#include "cljfunction.h"
#include "cljgrid.h"
#include "cljgroup.h"
#include "cljnbpairs.h"
#include "cljparam.h"
#include "cljprobe.h"
#include "cljrffunction.h"
#include "cljshiftfunction.h"
#include "cljworkspace.h"
#include "dihedral.h"
#include "dihedralrestraint.h"
#include "distancerestraint.h"
#include "excludedpairs.h"
#include "fouratomfunctions.h"
#include "gridff.h"
#include "gridff2.h"
#include "gromacsparams.h"
#include "improper.h"
#include "intercljff.h"
#include "intercoulombff.h"
#include "interff.h"
#include "intergroupff.h"
#include "interljff.h"
#include "internalcomponent.h"
#include "internalff.h"
#include "internalgroupff.h"
#include "internalparameters.h"
#include "internalperturbation.h"
#include "intersoftcljff.h"
#include "intracljff.h"
#include "intracoulombff.h"
#include "intraff.h"
#include "intragroupff.h"
#include "intraljff.h"
#include "intrasoftcljff.h"
#include "ljpair.h"
#include "ljparameter.h"
#include "ljperturbation.h"
#include "mmdetail.h"
#include "multicljcomponent.h"
#include "restraint.h"
#include "restraintcomponent.h"
#include "restraintff.h"
#include "selectorangle.h"
#include "selectorbond.h"
#include "selectordihedral.h"
#include "selectorimproper.h"
#include "selectormangle.h"
#include "selectormbond.h"
#include "selectormdihedral.h"
#include "selectormimproper.h"
#include "softcljcomponent.h"
#include "switchingfunction.h"
#include "threeatomfunctions.h"
#include "twoatomfunctions.h"

#include "Helpers/objectregistry.hpp"

void register_SireMM_objects()
{

    ObjectRegistry::registerConverterFor< SireMM::AmberParams >();
    ObjectRegistry::registerConverterFor< SireMM::AmberBond >();
    ObjectRegistry::registerConverterFor< SireMM::AmberAngle >();
    ObjectRegistry::registerConverterFor< SireMM::AmberDihPart >();
    ObjectRegistry::registerConverterFor< SireMM::AmberDihedral >();
    ObjectRegistry::registerConverterFor< SireMM::AmberNB14 >();
    ObjectRegistry::registerConverterFor< SireMM::AmberNBDihPart >();
    ObjectRegistry::registerConverterFor< SireMM::Angle >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::Angle> >();
    ObjectRegistry::registerConverterFor< SireMM::AngleRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::AtomLJs >();
    ObjectRegistry::registerConverterFor< SireMM::Bond >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::Bond> >();
    ObjectRegistry::registerConverterFor< SireMM::CLJ14Group >();
    ObjectRegistry::registerConverterFor< SireMM::CLJAtom >();
    ObjectRegistry::registerConverterFor< SireMM::CLJAtoms >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBox >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxIndex >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxes >();
    ObjectRegistry::registerConverterFor< SireMM::CLJBoxDistance >();
    ObjectRegistry::registerConverterFor< SireMM::CLJCalculator >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombComponent >();
    ObjectRegistry::registerConverterFor< SireMM::LJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::CLJDelta >();
    ObjectRegistry::registerConverterFor< SireMM::CLJExtractor >();
    ObjectRegistry::registerConverterFor< SireMM::NullCLJFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJGrid >();
    ObjectRegistry::registerConverterFor< SireMM::CLJGroup >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::LJScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::CLJScaleFactor >();
    ObjectRegistry::registerConverterFor< SireMM::CLJNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::LJNBPairs >();
    ObjectRegistry::registerConverterFor< SireMM::CLJParams >();
    ObjectRegistry::registerConverterFor< SireMM::CLJParamsArray >();
    ObjectRegistry::registerConverterFor< SireMM::CoulombProbe >();
    ObjectRegistry::registerConverterFor< SireMM::LJProbe >();
    ObjectRegistry::registerConverterFor< SireMM::CLJProbe >();
    ObjectRegistry::registerConverterFor< SireMM::CLJRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJIntraRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftIntraRFFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJIntraShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJSoftIntraShiftFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CLJWorkspace >();
    ObjectRegistry::registerConverterFor< SireMM::Dihedral >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::Dihedral> >();
    ObjectRegistry::registerConverterFor< SireMM::DihedralRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::DistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::DoubleDistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::TripleDistanceRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::HarmonicDistanceForceConstant >();
    ObjectRegistry::registerConverterFor< SireMM::ExcludedPairs >();
    ObjectRegistry::registerConverterFor< SireMM::FourAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::GridFF >();
    ObjectRegistry::registerConverterFor< SireMM::GridFF2 >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsAtomType >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsBond >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsAngle >();
    ObjectRegistry::registerConverterFor< SireMM::GromacsDihedral >();
    ObjectRegistry::registerConverterFor< SireMM::Improper >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::Improper> >();
    ObjectRegistry::registerConverterFor< SireMM::InterCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::BondComponent >();
    ObjectRegistry::registerConverterFor< SireMM::AngleComponent >();
    ObjectRegistry::registerConverterFor< SireMM::DihedralComponent >();
    ObjectRegistry::registerConverterFor< SireMM::ImproperComponent >();
    ObjectRegistry::registerConverterFor< SireMM::UreyBradleyComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchStretchComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchBendComponent >();
    ObjectRegistry::registerConverterFor< SireMM::BendBendComponent >();
    ObjectRegistry::registerConverterFor< SireMM::StretchBendTorsionComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14CoulombComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14LJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::Intra14Component >();
    ObjectRegistry::registerConverterFor< SireMM::InternalComponent >();
    ObjectRegistry::registerConverterFor< SireMM::InternalFF >();
    ObjectRegistry::registerConverterFor< SireMM::InternalGroupFF >();
    ObjectRegistry::registerConverterFor< SireMM::InternalParameters >();
    ObjectRegistry::registerConverterFor< SireMM::InternalParameters3D >();
    ObjectRegistry::registerConverterFor< SireMM::TwoAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::ThreeAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::FourAtomPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::InterSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::InterGroupSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCoulombFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupCoulombFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupSoftCLJFFBase >();
    ObjectRegistry::registerConverterFor< SireMM::IntraGroupSoftCLJFF >();
    ObjectRegistry::registerConverterFor< SireMM::LJPair >();
    ObjectRegistry::registerConverterFor< SireMM::LJParameter >();
    ObjectRegistry::registerConverterFor< SireMM::LJPerturbation >();
    ObjectRegistry::registerConverterFor< SireMM::MMDetail >();
    ObjectRegistry::registerConverterFor< SireMM::MultiCLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::NullRestraint >();
    ObjectRegistry::registerConverterFor< SireMM::RestraintComponent >();
    ObjectRegistry::registerConverterFor< SireMM::RestraintFF >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorAngle >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::SelectorAngle> >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorBond >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::SelectorBond> >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorDihedral >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::SelectorDihedral> >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorImproper >();
    ObjectRegistry::registerConverterFor< SireMol::Mover<SireMM::SelectorImproper> >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorMAngle >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorMBond >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorMDihedral >();
    ObjectRegistry::registerConverterFor< SireMM::SelectorMImproper >();
    ObjectRegistry::registerConverterFor< SireMM::SoftCLJComponent >();
    ObjectRegistry::registerConverterFor< SireMM::NoCutoff >();
    ObjectRegistry::registerConverterFor< SireMM::HarmonicSwitchingFunction >();
    ObjectRegistry::registerConverterFor< SireMM::CHARMMSwitchingFunction >();
    ObjectRegistry::registerConverterFor< SireMM::ThreeAtomFunctions >();
    ObjectRegistry::registerConverterFor< SireMM::TwoAtomFunctions >();

}

