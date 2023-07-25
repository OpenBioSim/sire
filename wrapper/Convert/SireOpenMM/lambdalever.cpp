/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2023  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the website
  *  at http://sire.openbiosim.org
  *
\*********************************************/

#include "lambdalever.h"

#include "SireCAS/values.h"

using namespace SireOpenMM;
using namespace SireCAS;

LambdaLever::LambdaLever() : SireBase::ConcreteProperty<LambdaLever, SireBase::Property>()
{
}

LambdaLever::LambdaLever(const LambdaLever &other)
    : SireBase::ConcreteProperty<LambdaLever, SireBase::Property>(other),
      name_to_ffidx(other.name_to_ffidx),
      lambda_schedule(other.lambda_schedule),
      perturbable_mols(other.perturbable_mols),
      start_indicies(other.start_indicies)
{
}

LambdaLever::~LambdaLever()
{
}

LambdaLever &LambdaLever::operator=(const LambdaLever &other)
{
    if (this != &other)
    {
        name_to_ffidx = other.name_to_ffidx;
        lambda_schedule = other.lambda_schedule;
        perturbable_mols = other.perturbable_mols;
        start_indicies = other.start_indicies;
        Property::operator=(other);
    }

    return *this;
}

bool LambdaLever::operator==(const LambdaLever &other) const
{
    return name_to_ffidx == other.name_to_ffidx and
           lambda_schedule == other.lambda_schedule and
           perturbable_mols == other.perturbable_mols and
           start_indicies == other.start_indicies;
}

bool LambdaLever::operator!=(const LambdaLever &other) const
{
    return not this->operator==(other);
}

LambdaLever *LambdaLever::clone() const
{
    return new LambdaLever(*this);
}

const char *LambdaLever::what() const
{
    return LambdaLever::typeName();
}

const char *LambdaLever::typeName()
{
    return QMetaType::typeName(qMetaTypeId<LambdaLever>());
}

template <class T>
T *get_force(const QString &name, OpenMM::Context &context,
             const QHash<QString, int> &name_to_index,
             const QString &force_type)
{
    auto it = name_to_index.constFind(name);

    if (it == name_to_index.constEnd())
    {
        return 0;
    }

    const int idx = it.value();

    OpenMM::System &system = const_cast<OpenMM::System &>(context.getSystem());

    const int num_forces = system.getNumForces();

    if (idx < 0 or idx >= num_forces)
    {
        throw SireError::invalid_key(QObject::tr(
                                         "The index for the Force called '%1', %2, is invalid for an "
                                         "OpenMM System which has %3 forces.")
                                         .arg(name)
                                         .arg(idx)
                                         .arg(num_forces),
                                     CODELOC);
    }

    OpenMM::Force &force = system.getForce(idx);

    T *t_force = dynamic_cast<T *>(&force);

    if (t_force == 0)
    {
        throw SireError::invalid_cast(QObject::tr(
                                          "Cannot case the force called '%1' to a %2.")
                                          .arg(name)
                                          .arg(force_type),
                                      CODELOC);
    }

    return t_force;
}

void LambdaLever::set_lambda(OpenMM::Context &context,
                             double lambda_value) const
{
    // go over each forcefield and update the parameters in the forcefield,
    // and then pass those updated parameters to the context

    // get copies of the forcefields in which the parameters will be changed
    auto cljff = get_force<OpenMM::NonbondedForce>("clj", context, name_to_ffidx, "NonbondedForce");
    auto bondff = get_force<OpenMM::HarmonicBondForce>("bond", context, name_to_ffidx, "HarmonicBondForce");
    auto angff = get_force<OpenMM::HarmonicAngleForce>("angle", context, name_to_ffidx, "HarmonicAngleForce");
    auto dihff = get_force<OpenMM::PeriodicTorsionForce>("torsion", context, name_to_ffidx, "PeriodicTorsionForce");

    // change the parameters for all of the perturbable molecules
    for (int i = 0; i < this->perturbable_mols.count(); ++i)
    {
        const auto &perturbable_mol = this->perturbable_mols[i];
        const auto &start_idxs = this->start_indicies[i];

        // calculate the new parameters for this lambda value
        const auto morphed_charges = this->lambda_schedule.morph(
            "charge",
            perturbable_mol.getCharges(),
            perturbable_mol.perturbed->getCharges(),
            lambda_value);

        const auto morphed_sigmas = this->lambda_schedule.morph(
            "sigma",
            perturbable_mol.getSigmas(),
            perturbable_mol.perturbed->getSigmas(),
            lambda_value);

        const auto morphed_epsilons = this->lambda_schedule.morph(
            "epsilon",
            perturbable_mol.getEpsilons(),
            perturbable_mol.perturbed->getEpsilons(),
            lambda_value);

        const auto morphed_bond_k = this->lambda_schedule.morph(
            "bond_k",
            perturbable_mol.getBondKs(),
            perturbable_mol.perturbed->getBondKs(),
            lambda_value);

        const auto morphed_bond_length = this->lambda_schedule.morph(
            "bond_length",
            perturbable_mol.getBondLengths(),
            perturbable_mol.perturbed->getBondLengths(),
            lambda_value);

        const auto morphed_angle_k = this->lambda_schedule.morph(
            "angle_k",
            perturbable_mol.getAngleKs(),
            perturbable_mol.perturbed->getAngleKs(),
            lambda_value);

        const auto morphed_angle_size = this->lambda_schedule.morph(
            "angle_size",
            perturbable_mol.getAngleSizes(),
            perturbable_mol.perturbed->getAngleSizes(),
            lambda_value);

        const auto morphed_torsion_phase = this->lambda_schedule.morph(
            "torsion_phase",
            perturbable_mol.getTorsionPhases(),
            perturbable_mol.perturbed->getTorsionPhases(),
            lambda_value);

        const auto morphed_torsion_periodicity = this->lambda_schedule.morph(
            "torsion_periodicity",
            perturbable_mol.getTorsionPeriodicities(),
            perturbable_mol.perturbed->getTorsionPeriodicities(),
            lambda_value);

        const auto morphed_torsion_k = this->lambda_schedule.morph(
            "torsion_k",
            perturbable_mol.getTorsionKs(),
            perturbable_mol.perturbed->getTorsionKs(),
            lambda_value);

        // now update the forcefields
        int start_index = start_idxs.value("clj", -1);

        if (start_index != -1)
        {
            const int nparams = morphed_charges.count();

            if (cljff == 0)
                throw SireError::incompatible_error(QObject::tr(
                                                        "There is no NonbondedForce called 'clj' in this context!"),
                                                    CODELOC);

            for (int j = 0; j < nparams; ++j)
            {
                cljff->setParticleParameters(start_index + j, morphed_charges[j], morphed_sigmas[j], morphed_epsilons[j]);
            }
        }

        start_index = start_idxs.value("bond", -1);

        if (start_index != -1)
        {
            const int nparams = morphed_bond_k.count();

            if (bondff == 0)
                throw SireError::incompatible_error(QObject::tr(
                                                        "There is no HarmonicBondForce called 'bond' in this context!"),
                                                    CODELOC);

            for (int j = 0; j < nparams; ++j)
            {
                const int index = start_index + j;

                int particle1, particle2;
                double length, k;

                bondff->getBondParameters(index, particle1, particle2,
                                          length, k);

                bondff->setBondParameters(index, particle1, particle2,
                                          morphed_bond_length[j],
                                          morphed_bond_k[j]);
            }
        }

        start_index = start_idxs.value("angle", -1);

        if (start_index != -1)
        {
            const int nparams = morphed_angle_k.count();

            if (angff == 0)
                throw SireError::incompatible_error(QObject::tr(
                                                        "There is no HarmonicAngleForce called 'angle' in this context!"),
                                                    CODELOC);

            for (int j = 0; j < nparams; ++j)
            {
                const int index = start_index + j;

                int particle1, particle2, particle3;
                double size, k;

                angff->getAngleParameters(index,
                                          particle1, particle2, particle3,
                                          size, k);

                angff->setAngleParameters(index,
                                          particle1, particle2, particle3,
                                          morphed_angle_size[j],
                                          morphed_angle_k[j]);
            }
        }

        start_index = start_idxs.value("torsion", -1);

        if (start_index != -1)
        {
            const int nparams = morphed_torsion_k.count();

            if (dihff == 0)
                throw SireError::incompatible_error(QObject::tr(
                                                        "There is no PeriodicTorsionForce called 'torsion' in this context!"),
                                                    CODELOC);

            for (int j = 0; j < nparams; ++j)
            {
                const int index = start_index + j;

                int particle1, particle2, particle3, particle4;
                double phase, k;
                int periodicity;

                dihff->getTorsionParameters(index,
                                            particle1, particle2,
                                            particle3, particle4,
                                            periodicity, phase, k);

                dihff->setTorsionParameters(index,
                                            particle1, particle2,
                                            particle3, particle4,
                                            morphed_torsion_periodicity[j],
                                            morphed_torsion_phase[j],
                                            morphed_torsion_k[j]);
            }
        }
    }

    if (cljff)
        cljff->updateParametersInContext(context);

    if (bondff)
        bondff->updateParametersInContext(context);

    if (angff)
        angff->updateParametersInContext(context);

    if (dihff)
        dihff->updateParametersInContext(context);
}

void LambdaLever::set_force_index(const QString &force, int index)
{
    if (index < 0)
        throw SireError::invalid_index(QObject::tr(
                                           "This looks like an invalid index for force '%1' : %2")
                                           .arg(force)
                                           .arg(index),
                                       CODELOC);

    this->name_to_ffidx.insert(force, index);
}

void LambdaLever::add_perturbable_molecule(const OpenMMMolecule &molecule,
                                           const QHash<QString, qint32> &starts)
{
    // should add in some sanity checks for these inputs
    this->perturbable_mols.append(molecule);
    this->start_indicies.append(starts);
}

LambdaSchedule LambdaLever::schedule() const
{
    return lambda_schedule;
}

void LambdaLever::set_schedule(const LambdaSchedule &sched)
{
    lambda_schedule = sched;
}
