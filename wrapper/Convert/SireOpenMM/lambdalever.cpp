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
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#include "emle.h"
#include "lambdalever.h"

#include "SireCAS/values.h"

#include "tostring.h"

using namespace SireOpenMM;
using namespace SireCAS;

LambdaLever::LambdaLever() : SireBase::ConcreteProperty<LambdaLever, SireBase::Property>()
{
}

LambdaLever::LambdaLever(const LambdaLever &other)
    : SireBase::ConcreteProperty<LambdaLever, SireBase::Property>(other),
      name_to_ffidx(other.name_to_ffidx),
      name_to_restraintidx(other.name_to_restraintidx),
      lambda_schedule(other.lambda_schedule),
      perturbable_mols(other.perturbable_mols),
      start_indices_pert(other.start_indices_pert),
      perturbable_maps(other.perturbable_maps)
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
        name_to_restraintidx = other.name_to_restraintidx;
        lambda_schedule = other.lambda_schedule;
        perturbable_mols = other.perturbable_mols;
        start_indices_pert = other.start_indices_pert;
        perturbable_maps = other.perturbable_maps;
        Property::operator=(other);
    }

    return *this;
}

bool LambdaLever::operator==(const LambdaLever &other) const
{
    return name_to_ffidx == other.name_to_ffidx and
           name_to_restraintidx == other.name_to_restraintidx and
           lambda_schedule == other.lambda_schedule and
           perturbable_mols == other.perturbable_mols and
           start_indices_pert == other.start_indices_pert and
           perturbable_maps == other.perturbable_maps;
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

bool LambdaLever::hasLever(const QString &lever_name)
{
    return this->lambda_schedule.getLevers().contains(lever_name);
}

void LambdaLever::addLever(const QString &lever_name)
{
    this->lambda_schedule.addLever(lever_name);
}

/** Get the index of the force called 'name'. Returns -1 if
 *  there is no force with this name
 */
int LambdaLever::getForceIndex(const QString &name) const
{
    auto it = name_to_ffidx.constFind(name);

    if (it == name_to_ffidx.constEnd())
    {
        return -1;
    }

    return it.value();
}

/** Get the C++ type of the force called 'name'. Returns an
 *  empty string if there is no such force
 */
QString LambdaLever::getForceType(const QString &name,
                                  const OpenMM::System &system) const
{
    int idx = this->getForceIndex(name);

    if (idx == -1)
        return QString();

    if (idx < 0 or idx >= system.getNumForces())
    {
        throw SireError::invalid_key(QObject::tr(
                                         "The index for the Force called '%1', %2, is invalid for an "
                                         "OpenMM System which has %3 forces.")
                                         .arg(name)
                                         .arg(idx)
                                         .arg(system.getNumForces()),
                                     CODELOC);
    }

    const OpenMM::Force &force = system.getForce(idx);

    return QString::fromStdString(force.getName());
}

std::tuple<int, int, double, double, double, double>
get_exception(int atom0, int atom1, int start_index,
              double coul_14_scl, double lj_14_scl,
              const QVector<double> &morphed_charges,
              const QVector<double> &morphed_sigmas,
              const QVector<double> &morphed_epsilons,
              const QVector<double> &morphed_alphas)
{
    double charge = 0.0;
    double sigma = 0.0;
    double epsilon = 0.0;
    double alpha = 0.0;

    if (coul_14_scl != 0 or lj_14_scl != 0)
    {
        const int n = morphed_charges.count();

        if (morphed_sigmas.count() != n or morphed_epsilons.count() != n)
        {
            throw SireError::program_bug(QObject::tr(
                                             "Expect same number of parameters! %1 vs %2 vs %3")
                                             .arg(n)
                                             .arg(morphed_sigmas.count())
                                             .arg(morphed_epsilons.count()),
                                         CODELOC);
        }

        if (atom0 < 0 or atom0 >= n or atom1 < 0 or atom1 >= n)
        {
            throw SireError::invalid_index(QObject::tr(
                                               "Wrong number of parameters (%1) for these atoms (%2, %3)")
                                               .arg(n)
                                               .arg(atom0)
                                               .arg(atom1),
                                           CODELOC);
        }

        charge = coul_14_scl * morphed_charges.constData()[atom0] * morphed_charges.constData()[atom1];
        sigma = 0.5 * (morphed_sigmas.constData()[atom0] + morphed_sigmas.constData()[atom1]);
        epsilon = lj_14_scl * std::sqrt(morphed_epsilons.constData()[atom0] * morphed_epsilons.constData()[atom1]);

        if (not morphed_alphas.isEmpty())
        {
            alpha = std::max(morphed_alphas[atom0], morphed_alphas[atom1]);
        }

        // clamp alpha between 0 and 1
        if (alpha < 0)
            alpha = 0;
        else if (alpha > 1)
            alpha = 1;
    }

    if (charge == 0 and epsilon == 0)
    {
        // openmm tries to optimise away zero parameters - this is an issue
        // as perturbation requires that we don't remove them!
        // If we don't do this, then we get a
        // "updateParametersInContext: The set of non-excluded exceptions has changed"
        /// exception when we update parameters in context
        sigma = 1e-9;
        epsilon = 1e-9;
    }

    return std::make_tuple(atom0 + start_index,
                           atom1 + start_index,
                           charge, sigma, epsilon,
                           alpha);
}

/** Set the value of lambda in the passed context. Returns the
 *  actual value of lambda set.
 */
double LambdaLever::setLambda(OpenMM::Context &context,
                              double lambda_value) const
{
    // go over each forcefield and update the parameters in the forcefield,
    // and then pass those updated parameters to the context
    if (this->lambda_schedule.isNull())
        return 0.0;

    lambda_value = this->lambda_schedule.clamp(lambda_value);

    // we need an editable reference to the system to get editable
    // pointers to the forces...
    OpenMM::System &system = const_cast<OpenMM::System &>(context.getSystem());

    // get copies of the forcefields in which the parameters will be changed
    auto qmff = this->getForce<QMMMForce>("qmff", system);
    auto cljff = this->getForce<OpenMM::NonbondedForce>("clj", system);
    auto ghost_ghostff = this->getForce<OpenMM::CustomNonbondedForce>("ghost/ghost", system);
    auto ghost_nonghostff = this->getForce<OpenMM::CustomNonbondedForce>("ghost/non-ghost", system);
    auto ghost_14ff = this->getForce<OpenMM::CustomBondForce>("ghost-14", system);
    auto bondff = this->getForce<OpenMM::HarmonicBondForce>("bond", system);
    auto angff = this->getForce<OpenMM::HarmonicAngleForce>("angle", system);
    auto dihff = this->getForce<OpenMM::PeriodicTorsionForce>("torsion", system);

    // we know if we have peturbable ghost atoms if we have the ghost forcefields
    const bool have_ghost_atoms = (ghost_ghostff != 0 or ghost_nonghostff != 0);

    std::vector<double> custom_params = {0.0, 0.0, 0.0, 0.0};

    if (qmff != 0)
    {
        double lam = this->lambda_schedule.morph("qmff", 1.0, 0.0, lambda_value);
        qmff->setLambda(lam);
    }

    // change the parameters for all of the perturbable molecules
    for (int i = 0; i < this->perturbable_mols.count(); ++i)
    {
        const auto &perturbable_mol = this->perturbable_mols[i];
        const auto &start_idxs = this->start_indices_pert[i];

        // calculate the new parameters for this lambda value
        const auto morphed_charges = this->lambda_schedule.morph(
            "charge",
            perturbable_mol.getCharges0(),
            perturbable_mol.getCharges1(),
            lambda_value);

        const auto morphed_sigmas = this->lambda_schedule.morph(
            "sigma",
            perturbable_mol.getSigmas0(),
            perturbable_mol.getSigmas1(),
            lambda_value);

        const auto morphed_epsilons = this->lambda_schedule.morph(
            "epsilon",
            perturbable_mol.getEpsilons0(),
            perturbable_mol.getEpsilons1(),
            lambda_value);

        const auto morphed_alphas = this->lambda_schedule.morph(
            "alpha",
            perturbable_mol.getAlphas0(),
            perturbable_mol.getAlphas1(),
            lambda_value);

        const auto morphed_bond_k = this->lambda_schedule.morph(
            "bond_k",
            perturbable_mol.getBondKs0(),
            perturbable_mol.getBondKs1(),
            lambda_value);

        const auto morphed_bond_length = this->lambda_schedule.morph(
            "bond_length",
            perturbable_mol.getBondLengths0(),
            perturbable_mol.getBondLengths1(),
            lambda_value);

        const auto morphed_angle_k = this->lambda_schedule.morph(
            "angle_k",
            perturbable_mol.getAngleKs0(),
            perturbable_mol.getAngleKs1(),
            lambda_value);

        const auto morphed_angle_size = this->lambda_schedule.morph(
            "angle_size",
            perturbable_mol.getAngleSizes0(),
            perturbable_mol.getAngleSizes1(),
            lambda_value);

        const auto morphed_torsion_phase = this->lambda_schedule.morph(
            "torsion_phase",
            perturbable_mol.getTorsionPhases0(),
            perturbable_mol.getTorsionPhases1(),
            lambda_value);

        const auto morphed_torsion_k = this->lambda_schedule.morph(
            "torsion_k",
            perturbable_mol.getTorsionKs0(),
            perturbable_mol.getTorsionKs1(),
            lambda_value);

        const auto morphed_charge_scale = this->lambda_schedule.morph(
            "charge_scale",
            perturbable_mol.getChargeScales0(),
            perturbable_mol.getChargeScales1(),
            lambda_value);

        const auto morphed_lj_scale = this->lambda_schedule.morph(
            "lj_scale",
            perturbable_mol.getLJScales0(),
            perturbable_mol.getLJScales1(),
            lambda_value);

        // now update the forcefields
        int start_index = start_idxs.value("clj", -1);

        if (start_index != -1 and cljff != 0)
        {
            const int nparams = morphed_charges.count();

            if (have_ghost_atoms)
            {
                for (int j = 0; j < nparams; ++j)
                {
                    const bool is_from_ghost = perturbable_mol.getFromGhostIdxs().contains(j);
                    const bool is_to_ghost = perturbable_mol.getToGhostIdxs().contains(j);

                    // reduced charge
                    custom_params[0] = morphed_charges[j];
                    // half_sigma
                    custom_params[1] = 0.5 * morphed_sigmas[j];
                    // two_sqrt_epsilon
                    custom_params[2] = 2.0 * std::sqrt(morphed_epsilons[j]);
                    // alpha
                    custom_params[3] = morphed_alphas[j];

                    // clamp alpha between 0 and 1
                    if (custom_params[3] < 0)
                        custom_params[3] = 0;
                    else if (custom_params[3] > 1)
                        custom_params[3] = 1;

                    ghost_ghostff->setParticleParameters(start_index + j, custom_params);
                    ghost_nonghostff->setParticleParameters(start_index + j, custom_params);

                    if (is_from_ghost or is_to_ghost)
                    {
                        // don't set the LJ parameters in the cljff
                        cljff->setParticleParameters(start_index + j, morphed_charges[j], 0.0, 0.0);
                    }
                    else
                    {
                        cljff->setParticleParameters(start_index + j, morphed_charges[j], morphed_sigmas[j], morphed_epsilons[j]);
                    }
                }
            }
            else
            {
                for (int j = 0; j < nparams; ++j)
                {
                    cljff->setParticleParameters(start_index + j, morphed_charges[j], morphed_sigmas[j], morphed_epsilons[j]);
                }
            }

            const auto idxs = perturbable_mol.getExceptionIndices("clj");

            if (not idxs.isEmpty())
            {
                const auto exception_atoms = perturbable_mol.getExceptionAtoms();

                for (int j = 0; j < exception_atoms.count(); ++j)
                {
                    const auto &atoms = exception_atoms[j];

                    const auto atom0 = std::get<0>(atoms);
                    const auto atom1 = std::get<1>(atoms);

                    const auto coul_14_scale = morphed_charge_scale[j];
                    const auto lj_14_scale = morphed_lj_scale[j];

                    const bool atom0_is_ghost = perturbable_mol.isGhostAtom(atom0);
                    const bool atom1_is_ghost = perturbable_mol.isGhostAtom(atom1);

                    const auto p = get_exception(atom0, atom1,
                                                start_index, coul_14_scale, lj_14_scale,
                                                morphed_charges, morphed_sigmas, morphed_epsilons,
                                                morphed_alphas);

                    // don't set LJ terms for ghost atoms
                    if (atom0_is_ghost or atom1_is_ghost)
                    {
                        cljff->setExceptionParameters(
                            std::get<0>(idxs[j]),
                            std::get<0>(p), std::get<1>(p),
                            std::get<2>(p), 1e-9, 1e-9);

                        if (coul_14_scale != 0 or lj_14_scale != 0)
                        {
                            // this is a 1-4 parameter - need to update
                            // the ghost 1-4 forcefield
                            int nbidx = std::get<1>(idxs[j]);

                            if (nbidx < 0)
                                throw SireError::program_bug(QObject::tr(
                                                                "Unset NB14 index for a ghost atom?"),
                                                            CODELOC);

                            if (ghost_14ff != 0)
                            {
                                // parameters are q, sigma, four_epsilon and alpha
                                std::vector<double> params14 =
                                    {std::get<2>(p), std::get<3>(p),
                                    4.0 * std::get<4>(p), std::get<5>(p)};

                                ghost_14ff->setBondParameters(nbidx,
                                                            std::get<0>(p),
                                                            std::get<1>(p),
                                                            params14);
                            }
                        }
                    }
                    else
                    {
                        cljff->setExceptionParameters(
                            std::get<0>(idxs[j]),
                            std::get<0>(p), std::get<1>(p),
                            std::get<2>(p), std::get<3>(p),
                            std::get<4>(p));
                    }
                }
            }
        }

        start_index = start_idxs.value("bond", -1);

        if (start_index != -1 and bondff != 0)
        {
            const int nparams = morphed_bond_k.count();

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

        if (start_index != -1 and angff != 0)
        {
            const int nparams = morphed_angle_k.count();

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

        if (start_index != -1 and dihff != 0)
        {
            const int nparams = morphed_torsion_k.count();

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
                                            periodicity,
                                            morphed_torsion_phase[j],
                                            morphed_torsion_k[j]);
            }
        }
    }

    // update the parameters in the context
    if (cljff)
        cljff->updateParametersInContext(context);

    if (ghost_ghostff)
        ghost_ghostff->updateParametersInContext(context);

    if (ghost_nonghostff)
        ghost_nonghostff->updateParametersInContext(context);

    if (ghost_14ff)
        ghost_14ff->updateParametersInContext(context);

    // in OpenMM 8.1beta updating the bond parameters past lambda=0.25
    // causes a "All Forces must have identical exclusions" error,
    // when running minimisation without h-bond constraints...
    if (bondff)
        bondff->updateParametersInContext(context);

    if (angff)
        angff->updateParametersInContext(context);

    if (dihff)
        dihff->updateParametersInContext(context);

    // now update any restraints that are scaled
    for (const auto &restraint : this->name_to_restraintidx.keys())
    {
        // restraints always morph between 1 and 1 (i.e. they fully
        // follow whatever is set by lambda, e.g. 'initial*lambda'
        // to switch them on, or `final*lambda` to switch them off)
        const double rho = lambda_schedule.morph(restraint,
                                                1.0, 1.0,
                                                lambda_value);

        for (auto &ff : this->getRestraints(restraint, system))
        {
            if (ff != 0)
            {
                this->updateRestraintInContext(*ff, rho, context);
            }
        }
    }

    return lambda_value;
}

/** Update the parameters for a CustomCompoundBondForce for scale factor 'rho'
 *  in the passed context */
void _update_restraint_in_context(OpenMM::CustomCompoundBondForce *ff, double rho,
                                  OpenMM::Context &context)
{
    if (ff == 0)
        throw SireError::invalid_cast(QObject::tr(
                                          "Unable to cast the restraint force to an OpenMM::CustomCompoundBondForce, "
                                          "despite it reporting that is was an object of this type..."),
                                      CODELOC);

    const int nbonds = ff->getNumBonds();

    if (nbonds == 0)
        // nothing to update
        return;

    const int nparams = ff->getNumPerBondParameters();

    if (nparams == 0)
        throw SireError::incompatible_error(QObject::tr(
                                                "Unable to set 'rho' for this restraint as it has no custom parameters!"),
                                            CODELOC);

    // we set the first parameter - we can see what the current value
    // is from the first restraint. This is because rho should be the
    // first parameter and have the same value for all restraints
    std::vector<double> custom_params;
    std::vector<int> particles;

    custom_params.resize(nparams);
    particles.resize(ff->getNumParticlesPerBond());

    ff->getBondParameters(0, particles, custom_params);

    if (custom_params[0] == rho)
        // nothing to do - it is already equal to this value
        return;

    for (int i = 0; i < nbonds; ++i)
    {
        ff->getBondParameters(i, particles, custom_params);
        custom_params[0] = rho;
        ff->setBondParameters(i, particles, custom_params);
    }

    ff->updateParametersInContext(context);
}

/** Update the parameters for a CustomBondForce for scale factor 'rho'
 *  in the passed context */
void _update_restraint_in_context(OpenMM::CustomBondForce *ff, double rho,
                                  OpenMM::Context &context)
{
    if (ff == 0)
        throw SireError::invalid_cast(QObject::tr(
                                          "Unable to cast the restraint force to an OpenMM::CustomBondForce, "
                                          "despite it reporting that is was an object of this type..."),
                                      CODELOC);

    const int nbonds = ff->getNumBonds();

    if (nbonds == 0)
        // nothing to update
        return;

    const int nparams = ff->getNumPerBondParameters();

    if (nparams == 0)
        throw SireError::incompatible_error(QObject::tr(
                                                "Unable to set 'rho' for this restraint as it has no custom parameters!"),
                                            CODELOC);

    // we set the first parameter - we can see what the current value
    // is from the first restraint. This is because rho should be the
    // first parameter and have the same value for all restraints
    std::vector<double> custom_params;
    custom_params.resize(nparams);
    int atom0, atom1;

    ff->getBondParameters(0, atom0, atom1, custom_params);

    if (custom_params[0] == rho)
        // nothing to do - it is already equal to this value
        return;

    for (int i = 0; i < nbonds; ++i)
    {
        ff->getBondParameters(i, atom0, atom1, custom_params);
        custom_params[0] = rho;
        ff->setBondParameters(i, atom0, atom1, custom_params);
    }

    ff->updateParametersInContext(context);
}

/** Internal function used to update the restraint force to use the
 *  supplied value of rho in the passed context */
void LambdaLever::updateRestraintInContext(OpenMM::Force &ff, double rho,
                                           OpenMM::Context &context) const
{
    // what is the type of this force...?
    const auto ff_type = ff.getName();

    if (ff_type == "CustomBondForce")
    {
        _update_restraint_in_context(
            dynamic_cast<OpenMM::CustomBondForce *>(&ff),
            rho, context);
    }
    else if (ff_type == "CustomCompoundBondForce")
    {
        _update_restraint_in_context(
            dynamic_cast<OpenMM::CustomCompoundBondForce *>(&ff),
            rho, context);
    }
    else
    {
        throw SireError::unknown_type(QObject::tr(
                                          "Unable to update the restraints for the passed force as it has "
                                          "an unknown type (%1). We currently only support a limited number "
                                          "of force types, e.g. CustomBondForce etc")
                                          .arg(QString::fromStdString(ff_type)),
                                      CODELOC);
    }
}

/** Set the index of the force called 'force' in the OpenMM System.
 *  There can only be one force with this name. Attempts to add
 *  a duplicate will cause an error to be raised.
 */
void LambdaLever::setForceIndex(const QString &force,
                                int index)
{
    if (index < 0)
        throw SireError::invalid_index(QObject::tr(
                                           "This looks like an invalid index for force '%1' : %2")
                                           .arg(force)
                                           .arg(index),
                                       CODELOC);

    this->name_to_ffidx.insert(force, index);
}

/** Add the index of a restraint force called 'restraint' in the
 *  OpenMM System. There can be multiple restraint forces with
 *  the same name
 */
void LambdaLever::addRestraintIndex(const QString &restraint,
                                    int index)
{
    this->name_to_restraintidx.insertMulti(restraint, index);
}

/** Return the pointers to all of the forces from the passed System
 *  are restraints called 'restraint'. This returns an empty
 *  list if there are no restraints with this name */
QList<OpenMM::Force *> LambdaLever::getRestraints(const QString &name,
                                                  OpenMM::System &system) const
{
    QList<OpenMM::Force *> forces;

    const int num_forces = system.getNumForces();

    for (const auto &idx : this->name_to_restraintidx.values(name))
    {
        if (idx < 0 or idx >= num_forces)
        {
            throw SireError::invalid_key(QObject::tr(
                                             "The index for the Restraint Force called '%1', %2, is invalid for an "
                                             "OpenMM System which has %3 forces.")
                                             .arg(name)
                                             .arg(idx)
                                             .arg(num_forces),
                                         CODELOC);
        }

        forces.append(&(system.getForce(idx)));
    }

    return forces;
}

/** Add info for the passed perturbable OpenMMMolecule, returning
 *  its index in the list of perturbable molecules
 */
int LambdaLever::addPerturbableMolecule(const OpenMMMolecule &molecule,
                                        const QHash<QString, qint32> &starts)
{
    // should add in some sanity checks for these inputs
    this->perturbable_mols.append(PerturbableOpenMMMolecule(molecule));
    this->start_indices_pert.append(starts);
    this->perturbable_maps.insert(molecule.number, molecule.perturtable_map);
    return this->perturbable_mols.count() - 1;
}


/** Set the exception indices for the perturbable molecule at
 *  index 'mol_idx'
 */
void LambdaLever::setExceptionIndices(int mol_idx,
                                      const QString &name,
                                      const QVector<std::pair<int, int>> &exception_idxs)
{
    mol_idx = SireID::Index(mol_idx).map(this->perturbable_mols.count());
    this->perturbable_mols[mol_idx].setExceptionIndices(name, exception_idxs);
}

/** Return all of the property maps used to find the perturbable properties
 *  of the perturbable molecules. This is indexed by molecule number
 */
QHash<SireMol::MolNum, SireBase::PropertyMap> LambdaLever::getPerturbableMoleculeMaps() const
{
    return perturbable_maps;
}

LambdaSchedule LambdaLever::getSchedule() const
{
    return lambda_schedule;
}

void LambdaLever::setSchedule(const LambdaSchedule &sched)
{
    QStringList levers = lambda_schedule.getLevers();

    lambda_schedule = sched;

    for (const auto &lever : levers)
    {
        lambda_schedule.addLever(lever);
    }
}
