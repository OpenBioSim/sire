//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!
#include <Python.h>

#include "SireMove_registrars.h"

#include "supramove.h"
#include "volumechanger.h"
#include "uniformsampler.h"
#include "moves.h"
#include "rbworkspace.h"
#include "suprasubmoves.h"
#include "repexmove.h"
#include "velocityverlet.h"
#include "move.h"
#include "hybridmc.h"
#include "zmatrix.h"
#include "prefsampler.h"
#include "rigidbodymc.h"
#include "weightedmoves.h"
#include "replica.h"
#include "replicas.h"
#include "zmatmove.h"
#include "moleculardynamics.h"
#include "suprasubmove.h"
#include "supramoves.h"
#include "titrationmove.h"
#include "openmmfrenergyst.h"
#include "rbworkspacejm.h"
#include "suprasubsystem.h"
#include "openmmmdintegrator.h"
#include "suprasystem.h"
#include "integrator.h"
#include "mtsmc.h"
#include "suprasimpacket.h"
#include "internalmove.h"
#include "titrator.h"
#include "openmmfrenergydt.h"
#include "velocitygenerator.h"
#include "integratorworkspacejm.h"
#include "repexmove2.h"
#include "integratorworkspace.h"
#include "getpoint.h"
#include "simpacket.h"
#include "moldeleter.h"
#include "molinserter.h"
#include "flexibility.h"
#include "volumemove.h"
#include "suprasubsimpacket.h"
#include "internalmovesingle.h"
#include "simstore.h"
#include "dlmrigidbody.h"
#include "ensemble.h"

#include "Helpers/objectregistry.hpp"

void register_SireMove_objects()
{

    ObjectRegistry::registerConverterFor< SireMove::NullSupraMove >();
    ObjectRegistry::registerConverterFor< SireMove::NullVolumeChanger >();
    ObjectRegistry::registerConverterFor< SireMove::ScaleVolumeFromCenter >();
    ObjectRegistry::registerConverterFor< SireMove::UniformSampler >();
    ObjectRegistry::registerConverterFor< SireMove::SameMoves >();
    ObjectRegistry::registerConverterFor< SireMove::RBWorkspace >();
    ObjectRegistry::registerConverterFor< SireMove::SameSupraSubMoves >();
    ObjectRegistry::registerConverterFor< SireMove::RepExMove >();
    ObjectRegistry::registerConverterFor< SireMove::RepExSubMove >();
    ObjectRegistry::registerConverterFor< SireMove::VelocityVerlet >();
    ObjectRegistry::registerConverterFor< SireMove::NullMove >();
    ObjectRegistry::registerConverterFor< SireMove::HybridMC >();
    ObjectRegistry::registerConverterFor< SireMove::HMCGenerator >();
    ObjectRegistry::registerConverterFor< SireMove::ZMatrix >();
    ObjectRegistry::registerConverterFor< SireMove::ZMatrixLine >();
    ObjectRegistry::registerConverterFor< SireMove::ZMatrixCoords >();
    ObjectRegistry::registerConverterFor< SireMove::ZMatrixCoordsLine >();
    ObjectRegistry::registerConverterFor< SireMove::PrefSampler >();
    ObjectRegistry::registerConverterFor< SireMove::RigidBodyMC >();
    ObjectRegistry::registerConverterFor< SireMove::WeightedMoves >();
    ObjectRegistry::registerConverterFor< SireMove::Replica >();
    ObjectRegistry::registerConverterFor< SireMove::Replicas >();
    ObjectRegistry::registerConverterFor< SireMove::ZMatMove >();
    ObjectRegistry::registerConverterFor< SireMove::MolecularDynamics >();
    ObjectRegistry::registerConverterFor< SireMove::NullSupraSubMove >();
    ObjectRegistry::registerConverterFor< SireMove::SameSupraMoves >();
    ObjectRegistry::registerConverterFor< SireMove::TitrationMove >();
    ObjectRegistry::registerConverterFor< SireMove::OpenMMFrEnergyST >();
    ObjectRegistry::registerConverterFor< SireMove::OpenMMFrEnergyST >();
    ObjectRegistry::registerConverterFor< SireMove::RBWorkspaceJM >();
    ObjectRegistry::registerConverterFor< SireMove::SupraSubSystem >();
    ObjectRegistry::registerConverterFor< SireMove::OpenMMMDIntegrator >();
    ObjectRegistry::registerConverterFor< SireMove::OpenMMMDIntegrator >();
    ObjectRegistry::registerConverterFor< SireMove::SupraSystem >();
    ObjectRegistry::registerConverterFor< SireMove::NullIntegrator >();
    ObjectRegistry::registerConverterFor< SireMove::MTSMC >();
    ObjectRegistry::registerConverterFor< SireMove::SupraSimPacket >();
    ObjectRegistry::registerConverterFor< SireMove::InternalMove >();
    ObjectRegistry::registerConverterFor< SireMove::Titrator >();
    ObjectRegistry::registerConverterFor< SireMove::OpenMMFrEnergyDT >();
    ObjectRegistry::registerConverterFor< SireMove::OpenMMFrEnergyDT >();
    ObjectRegistry::registerConverterFor< SireMove::NullVelocityGenerator >();
    ObjectRegistry::registerConverterFor< SireMove::VelocitiesFromProperty >();
    ObjectRegistry::registerConverterFor< SireMove::MaxwellBoltzmann >();
    ObjectRegistry::registerConverterFor< SireMove::NullIntegratorWorkspaceJM >();
    ObjectRegistry::registerConverterFor< SireMove::AtomicVelocityWorkspaceJM >();
    ObjectRegistry::registerConverterFor< SireMove::RepExMove2 >();
    ObjectRegistry::registerConverterFor< SireMove::NullIntegratorWorkspace >();
    ObjectRegistry::registerConverterFor< SireMove::AtomicVelocityWorkspace >();
    ObjectRegistry::registerConverterFor< SireMove::NullGetPoint >();
    ObjectRegistry::registerConverterFor< SireMove::GetCOMPoint >();
    ObjectRegistry::registerConverterFor< SireMove::GetCOGPoint >();
    ObjectRegistry::registerConverterFor< SireMove::GetCentroidPoint >();
    ObjectRegistry::registerConverterFor< SireMove::SimPacket >();
    ObjectRegistry::registerConverterFor< SireMove::NullDeleter >();
    ObjectRegistry::registerConverterFor< SireMove::SpecifiedGroupsDeleter >();
    ObjectRegistry::registerConverterFor< SireMove::SystemWideDeleter >();
    ObjectRegistry::registerConverterFor< SireMove::NullInserter >();
    ObjectRegistry::registerConverterFor< SireMove::UniformInserter >();
    ObjectRegistry::registerConverterFor< SireMove::DofID >();
    ObjectRegistry::registerConverterFor< SireMove::Flexibility >();
    ObjectRegistry::registerConverterFor< SireMove::VolumeMove >();
    ObjectRegistry::registerConverterFor< SireMove::SupraSubSimPacket >();
    ObjectRegistry::registerConverterFor< SireMove::InternalMoveSingle >();
    ObjectRegistry::registerConverterFor< SireMove::SimStore >();
    ObjectRegistry::registerConverterFor< SireMove::DLMRigidBody >();
    ObjectRegistry::registerConverterFor< SireMove::Ensemble >();

}

