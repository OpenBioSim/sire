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
T &get_force(const QString &name, OpenMM::Context &context,
             const QHash<QString, int> &name_to_index,
             const QString &force_type)
{
    auto it = name_to_index.constFind(name);

    if (it == name_to_index.constEnd())
    {
        throw SireError::invalid_key(QObject::tr(
                                         "There is no Force called '%1' in the passed OpenMM Context. "
                                         "Available forces are [ %2 ]")
                                         .arg(name)
                                         .arg(name_to_index.keys().join(", ")),
                                     CODELOC);
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

    return *t_force;
}

void LambdaLever::set_lambda(OpenMM::Context &context,
                             double lambda_value) const
{
    // go over each forcefield and update the parameters in the forcefield,
    // and then pass those updated parameters to the context

    // get copies of the forcefields in which the parameters will be changed
    auto &cljff = get_force<OpenMM::NonbondedForce>("clj", context, name_to_ffidx, "NonbondedForce");

    // change the parameters for all of the perturbable molecules
    for (int i = 0; i < this->perturbable_mols.count(); ++i)
    {
        const auto &perturbable_mol = this->perturbable_mols[i];
        const auto &start_idxs = this->start_indicies[i];

        const auto charge0 = perturbable_mol.getCharges();
        const auto charge1 = perturbable_mol.perturbed->getCharges();

        const auto sigma0 = perturbable_mol.getSigmas();
        const auto sigma1 = perturbable_mol.perturbed->getSigmas();

        const auto epsilon0 = perturbable_mol.getEpsilons();
        const auto epsilon1 = perturbable_mol.getEpsilons();

        // calculate the new parameters for this lambda value
        const auto morphed_charges = this->lambda_schedule.morph("charge", charge0, charge1, lambda_value);
        const auto morphed_sigmas = this->lambda_schedule.morph("sigma", sigma0, sigma1, lambda_value);
        const auto morphed_epsilons = this->lambda_schedule.morph("epsilon", epsilon0, epsilon1, lambda_value);

        // now update the forcefield
        const auto start_index = start_idxs.value("clj", -1);

        if (start_index == -1)
            throw SireError::program_bug(QObject::tr(
                                             "No start index for force %1, despite this force existing?")
                                             .arg("clj"),
                                         CODELOC);

        const int nparams = morphed_charges.count();

        for (int j = 0; j < nparams; ++j)
        {
            int idx = start_index + j;

            // qDebug() << idx << j << "I" << charge0[j] << sigma0[j] << epsilon0[j];
            // qDebug() << idx << j << "M" << morphed_charges[j] << morphed_sigmas[j] << morphed_epsilons[j];
            // qDebug() << idx << j << "F" << charge1[j] << sigma1[j] << epsilon1[j];

            cljff.setParticleParameters(idx, morphed_charges[j], morphed_sigmas[j], morphed_epsilons[j]);
        }
    }

    cljff.updateParametersInContext(context);
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
