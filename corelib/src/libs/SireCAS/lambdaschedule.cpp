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

#include "lambdaschedule.h"

#include "SireCAS/values.h"

#include "SireBase/console.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;

QString _get_lever_name(QString force, QString lever)
{
    force = force.trimmed().simplified().replace(" ", "_").replace(":", ".");
    lever = lever.trimmed().simplified().replace(" ", "_").replace(":", ".");

    return force + "::" + lever;
}

QString _fix_lever_name(const QString &lever)
{
    if (lever.contains("::"))
    {
        return lever;
    }
    else
    {
        return "*::" + lever;
    }
}

static RegisterMetaType<LambdaSchedule> r_schedule;

QDataStream &operator<<(QDataStream &ds, const LambdaSchedule &schedule)
{
    writeHeader(ds, r_schedule, 3);

    SharedDataStream sds(ds);

    sds << schedule.constant_values
        << schedule.force_names
        << schedule.lever_names << schedule.stage_names
        << schedule.default_equations
        << schedule.stage_equations
        << schedule.mol_schedules
        << static_cast<const Property &>(schedule);

    return ds;
}

Symbol LambdaSchedule::lambda_symbol("λ");
Symbol LambdaSchedule::initial_symbol("initial");
Symbol LambdaSchedule::final_symbol("final");

Expression LambdaSchedule::default_morph_equation = (1.0 - LambdaSchedule::lam()) * LambdaSchedule::initial() +
                                                    LambdaSchedule::lam() * LambdaSchedule::final();

QDataStream &operator>>(QDataStream &ds, LambdaSchedule &schedule)
{
    VersionID v = readHeader(ds, r_schedule);

    if (v == 1 or v == 2 or v == 3)
    {
        SharedDataStream sds(ds);

        sds >> schedule.constant_values;

        if (v == 3)
            sds >> schedule.force_names;

        sds >> schedule.lever_names >> schedule.stage_names >>
            schedule.default_equations >> schedule.stage_equations;

        if (v == 2 or v == 3)
            sds >> schedule.mol_schedules;

        if (v < 3)
        {
            // need to make sure that the lever names are namespaced
            auto fixed_lever_names = QStringList();

            for (auto &lever : schedule.lever_names)
            {
                fixed_lever_names.append(_fix_lever_name(lever));

                for (auto &stage_equations : schedule.stage_equations)
                {
                    if (stage_equations.contains(lever))
                    {
                        auto fixed_lever = _fix_lever_name(lever);
                        stage_equations[fixed_lever] = stage_equations.take(lever);
                    }
                }
            }

            schedule.lever_names = fixed_lever_names;
        }

        sds >> static_cast<Property &>(schedule);

        for (auto &expression : schedule.default_equations)
        {
            if (expression == LambdaSchedule::default_morph_equation)
                expression = LambdaSchedule::default_morph_equation;
        }

        for (auto &stage_equations : schedule.stage_equations)
        {
            for (auto &expression : stage_equations)
            {
                if (expression == LambdaSchedule::default_morph_equation)
                    expression = LambdaSchedule::default_morph_equation;
            }
        }
    }
    else
        throw version_error(v, "1, 2, 3", r_schedule, CODELOC);

    return ds;
}

LambdaSchedule::LambdaSchedule() : ConcreteProperty<LambdaSchedule, Property>()
{
}

LambdaSchedule::LambdaSchedule(const LambdaSchedule &other)
    : ConcreteProperty<LambdaSchedule, Property>(other),
      mol_schedules(other.mol_schedules),
      constant_values(other.constant_values),
      force_names(other.force_names),
      lever_names(other.lever_names), stage_names(other.stage_names),
      default_equations(other.default_equations),
      stage_equations(other.stage_equations)
{
}

LambdaSchedule::~LambdaSchedule()
{
}

LambdaSchedule &LambdaSchedule::operator=(const LambdaSchedule &other)
{
    if (this != &other)
    {
        mol_schedules = other.mol_schedules;
        constant_values = other.constant_values;
        force_names = other.force_names;
        lever_names = other.lever_names;
        stage_names = other.stage_names;
        default_equations = other.default_equations;
        stage_equations = other.stage_equations;
        Property::operator=(other);
    }

    return *this;
}

bool LambdaSchedule::operator==(const LambdaSchedule &other) const
{
    return mol_schedules == other.mol_schedules and
           force_names == other.force_names and
           constant_values == other.constant_values and
           lever_names == other.lever_names and
           stage_names == other.stage_names and
           default_equations == other.default_equations and
           stage_equations == other.stage_equations;
}

bool LambdaSchedule::operator!=(const LambdaSchedule &other) const
{
    return not this->operator==(other);
}

LambdaSchedule *LambdaSchedule::clone() const
{
    return new LambdaSchedule(*this);
}

const char *LambdaSchedule::what() const
{
    return LambdaSchedule::typeName();
}

const char *LambdaSchedule::typeName()
{
    return QMetaType::typeName(qMetaTypeId<LambdaSchedule>());
}

bool LambdaSchedule::isNull() const
{
    return this->stage_names.isEmpty() and this->lever_names.isEmpty();
}

QString LambdaSchedule::toString() const
{
    if (this->isNull())
        return QObject::tr("LambdaSchedule::null");

    QStringList lines;

    for (int i = 0; i < this->stage_names.count(); ++i)
    {

        lines.append(QString("  %1: %2")
                         .arg(this->stage_names[i])
                         .arg(this->default_equations[i].toOpenMMString()));

        auto keys = this->stage_equations[i].keys();
        std::sort(keys.begin(), keys.end());

        for (const auto &lever : keys)
        {
            auto output_name = lever;
            output_name.replace("*::", "");
            lines.append(QString("    %1: %2")
                             .arg(output_name)
                             .arg(this->stage_equations[i][lever].toOpenMMString()));
        }
    }

    for (const auto &constant : this->constant_values.keys())
    {
        lines.append(QString("  %1 == %2")
                         .arg(constant.toString())
                         .arg(this->constant_values[constant]));
    }

    if (not this->mol_schedules.isEmpty())
    {
        lines.append("  Molecule schedules:");

        for (const auto &mol_id : this->mol_schedules.keys())
        {
            lines.append(QString("    %1: %2")
                             .arg(mol_id)
                             .arg(this->mol_schedules[mol_id].toString()));
        }
    }

    return QObject::tr("LambdaSchedule(\n%1\n)")
        .arg(lines.join("\n"));
}

/** Return a LambdaSchedule that represents a standard morph,
 *  where every forcefield parameter is scaled by
 *  (1-:lambda:).initial + :lambda:.final
 */
LambdaSchedule LambdaSchedule::standard_morph()
{
    LambdaSchedule l;
    l.addMorphStage();
    return l;
}

/** Return a LambdaSchedule that represents a central morph
 *  stage that is sandwiched between a charge descaling,
 *  and a charge rescaling stage. The first stage scales
 *  the "charge" lever down from 1.0 to `scale`. This
 *  is followed by a standard morph stage using the
 *  descaled charges. This the finished with a recharging
 *  stage that restores the charges back to their
 *  original values.
 */
LambdaSchedule LambdaSchedule::charge_scaled_morph(double scale)
{
    LambdaSchedule l;
    l.addMorphStage();
    l.addChargeScaleStages(scale);

    return l;
}

/** Return a schedule that can be used for a standard double-decoupling
 *  free energy perturbation. If `perturbed_is_decoupled` is true, then
 *  the perturbed state is decoupled, otherwise the reference state is
 *  decoupled.
 */
LambdaSchedule LambdaSchedule::standard_decouple(bool perturbed_is_decoupled)
{
    LambdaSchedule l;
    l.addDecoupleStage(perturbed_is_decoupled);

    return l;
}

/** Return a schedule that can be used for a standard double-decoupling
 *  free energy perturbation. If `perturbed_is_decoupled` is true, then
 *  the perturbed state is decoupled, otherwise the reference state is
 *  decoupled. In this case also add states to decharge and recharge
 *  the molecule either side of the decoupling stage, where the charges
 *  are scaled to 'scale' times their original value.
 */
LambdaSchedule LambdaSchedule::charge_scaled_decouple(double scale, bool perturbed_is_decoupled)
{
    LambdaSchedule l;
    l.addDecoupleStage(perturbed_is_decoupled);
    l.addChargeScaleStages(scale);

    return l;
}

/** Return a schedule that can be used for a standard double-annihilation
 *  free energy perturbation. If `perturbed_is_annihilated` is true, then
 *  the perturbed state is annihilated, otherwise the reference state is
 *  annihilated.
 */
LambdaSchedule LambdaSchedule::standard_annihilate(bool perturbed_is_annihilated)
{
    LambdaSchedule l;
    l.addAnnihilateStage(perturbed_is_annihilated);

    return l;
}

/** Return a schedule that can be used for a standard double-annihilation
 *  free energy perturbation. If `perturbed_is_annihilated` is true, then
 *  the perturbed state is annihilated, otherwise the reference state is
 *  annihilated. In this case also add states to decharge and recharge
 *  the molecule either side of the annihilation stage, where the charges
 *  are scaled to 'scale' times their original value.
 */
LambdaSchedule LambdaSchedule::charge_scaled_annihilate(double scale, bool perturbed_is_annihilated)
{
    LambdaSchedule l;
    l.addAnnihilateStage(perturbed_is_annihilated);
    l.addChargeScaleStages(scale);

    return l;
}

/** Return the symbol used to represent the :lambda: coordinate.
 *  This symbol is used to represent the per-stage :lambda:
 *  variable that goes from 0.0-1.0 within that stage.
 */
Symbol LambdaSchedule::lam()
{
    return lambda_symbol;
}

/** Return the symbol used to represent the initial
 *  (:lambda:=0) value of the forcefield parameter
 */
Symbol LambdaSchedule::initial()
{
    return initial_symbol;
}

/** Return the symbol used to represent the final
 *  (:lambda:=1) value of the forcefield parameter
 */
Symbol LambdaSchedule::final()
{
    return final_symbol;
}

/** Set the value of a constant that may be used in any
 *  of the stage equations.
 */
SireCAS::Symbol LambdaSchedule::setConstant(const QString &constant, double value)
{
    return this->setConstant(this->getConstantSymbol(constant), value);
}

/** Set the value of a constant that may be used in any
 *  of the stage equations.
 */
SireCAS::Symbol LambdaSchedule::setConstant(const SireCAS::Symbol &constant,
                                            double value)
{
    this->constant_values.set(constant, value);
    return constant;
}

/** Return the value of the passed constant that may be
 *  used in any of the stage equations
 */
double LambdaSchedule::getConstant(const QString &constant)
{
    return this->getConstant(this->getConstantSymbol(constant));
}

/** Return the value of the passed constant that may be
 *  used in any of the stage equations
 */
double LambdaSchedule::getConstant(const SireCAS::Symbol &constant) const
{
    return this->constant_values.value(constant);
}

/** Get the Symbol used to represent the named constant 'constant' */
SireCAS::Symbol LambdaSchedule::getConstantSymbol(const QString &constant) const
{
    return SireCAS::Symbol(constant);
}

/** Add a lever to the schedule. This is only useful if you want to
 *  plot how the equations would affect the lever. Levers will be
 *  automatically added by any perturbation run that needs them,
 *  so you don't need to add them manually yourself.
 */
void LambdaSchedule::addLever(const QString &lever)
{
    if (lever == "*" or this->lever_names.contains(lever))
        return;

    this->lever_names.append(lever);
}

/** Add some levers to the schedule. This is only useful if you want to
 *  plot how the equations would affect the lever. Levers will be
 *  automatically added by any perturbation run that needs them,
 *  so you don't need to add them manually yourself.
 */
void LambdaSchedule::addLevers(const QStringList &levers)
{
    for (const auto &lever : levers)
    {
        if (not(lever == "*" or this->lever_names.contains(lever)))
            this->lever_names.append(lever);
    }
}

/** Remove a lever from the schedule. This will not impact any
 *  perturbation runs that use this schedule, as any missing
 *  levers will be re-added.
 */
void LambdaSchedule::removeLever(const QString &lever)
{
    if (not this->lever_names.contains(lever))
        return;

    int idx = this->lever_names.indexOf(lever);

    this->default_equations.removeAt(idx);
    this->stage_equations.removeAt(idx);
    this->stage_names.removeAt(idx);
}

/** Remove some levers from the schedule. This will not impact any
 *  perturbation runs that use this schedule, as any missing
 *  levers will be re-added.
 */
void LambdaSchedule::removeLevers(const QStringList &levers)
{
    for (const auto &lever : levers)
    {
        this->removeLever(lever);
    }
}

/** Return the number of levers that have been explicitly added
 *  to the schedule. Note that levers will be automatically added
 *  by any perturbation run that needs them, so you don't normally
 *  need to manage them manually yourself.
 */
int LambdaSchedule::nLevers() const
{
    return this->lever_names.count();
}

/** Return all of the levers that have been explicitly added
 *  to the schedule. Note that levers will be automatically added
 *  by any perturbation run that needs them, so you don't normally
 *  need to manage them manually yourself.
 */
QStringList LambdaSchedule::getLevers() const
{
    return this->lever_names;
}

/** Add a force to a schedule. This is only useful if you want to
 *  plot how the equations would affect the lever. Forces will be
 *  automatically added by any perturbation run that needs them,
 *  so you don't need to add them manually yourself.
 */
void LambdaSchedule::addForce(const QString &force)
{
    if (force == "*" or this->force_names.contains(force))
        return;

    this->force_names.append(force);
}

/** Add some forces to a schedule. This is only useful if you want to
 *  plot how the equations would affect the lever. Forces will be
 *  automatically added by any perturbation run that needs them,
 *  so you don't need to add them manually yourself.
 */
void LambdaSchedule::addForces(const QStringList &forces)
{
    for (const auto &force : forces)
    {
        if (not(force == "*" or this->force_names.contains(force)))
            this->force_names.append(force);
    }
}

/** Remove a force from a schedule. This will not impact any
 *  perturbation runs that use this schedule, as any missing
 *  forces will be re-added.
 */
void LambdaSchedule::removeForce(const QString &force)
{
    if (not this->force_names.contains(force))
        return;

    int idx = this->force_names.indexOf(force);

    this->force_names.removeAt(idx);
}

/** Remove some forces from a schedule. This will not impact any
 *  perturbation runs that use this schedule, as any missing
 *  forces will be re-added.
 */
void LambdaSchedule::removeForces(const QStringList &forces)
{
    for (const auto &force : forces)
    {
        this->removeForce(force);
    }
}

/** Return the number of forces that have been explicitly added
 *  to the schedule. Note that forces will be automatically added
 *  by any perturbation run that needs them, so you don't normally
 *  need to manage them manually yourself.
 */
int LambdaSchedule::nForces() const
{
    return this->force_names.count();
}

/** Return all of the forces that have been explicitly added
 *  to the schedule. Note that forces will be automatically added
 *  by any perturbation run that needs them, so you don't normally
 *  need to manage them manually yourself.
 */
QStringList LambdaSchedule::getForces() const
{
    return this->force_names;
}

/** Return the number of stages in this schedule */
int LambdaSchedule::nStages() const
{
    return this->stage_names.count();
}

/** Return the names of all of the stages in this schedule, in
 *  the order they will be performed
 */
QStringList LambdaSchedule::getStages() const
{
    return this->stage_names;
}

/** Clamp and return the passed lambda value so that it is between a valid
 *  range for this schedule (typically between [0.0-1.0] inclusive).
 */
double LambdaSchedule::clamp(double lambda_value) const
{
    return std::max(0.0, std::min(lambda_value, 1.0));
}

/** Internal function used to resolve the global value of lambda down
 *  to a stage number and stage-lambda value
 */
std::tuple<int, double> LambdaSchedule::resolve_lambda(double lambda_value) const
{
    if (this->nStages() == 0)
        throw SireError::invalid_key(QObject::tr(
                                         "Cannot resolve λ to a stage as no stages have been added. "
                                         "Add a stage via the 'add_stage' function and try again."),
                                     CODELOC);

    lambda_value = this->clamp(lambda_value);

    if (lambda_value == 0.0)
    {
        return std::tuple<int, double>(0, 0.0);
    }
    else if (lambda_value == 1.0)
    {
        return std::tuple<int, double>(this->nStages() - 1, 1.0);
    }

    double stage_width = 1.0 / this->nStages();

    double resolved = lambda_value / stage_width;

    double stage = std::floor(resolved);

    return std::tuple<int, double>(int(stage), resolved - stage);
}

/** Return the name of the stage that controls the forcefield parameters
 *  at the global value of :lambda: equal to `lambda_value`
 */
QString LambdaSchedule::getStage(double lambda_value) const
{
    auto resolved = this->resolve_lambda(lambda_value);

    return this->stage_names[std::get<0>(resolved)];
}

/** Return the stage-local value of :lambda: that corresponds to the
 *  global value of :lambda: at `lambda_value`
 */
double LambdaSchedule::getLambdaInStage(double lambda_value) const
{
    if (this->nStages() == 0)
        return this->clamp(lambda_value);

    auto resolved = this->resolve_lambda(lambda_value);

    return std::get<1>(resolved);
}

/** Completely clear all stages and levers */
void LambdaSchedule::clear()
{
    this->stage_names.clear();
    this->stage_equations.clear();
    this->default_equations.clear();
    this->constant_values = Values();
}

/** Append a morph stage onto this schedule. The morph stage is a
 *  standard stage that scales each forcefield parameter by
 *  (1-:lambda:).initial + :lambda:.final
 */
void LambdaSchedule::addMorphStage(const QString &name)
{
    this->addStage(name, default_morph_equation);
}

/** Append a morph stage onto this schedule. The morph stage is a
 *  standard stage that scales each forcefield parameter by
 *  (1-:lambda:).initial + :lambda:.final
 */
void LambdaSchedule::addMorphStage()
{
    this->addMorphStage("morph");
}

/** Add a stage to the schedule that will decouple the perturbed
 *  state if `perturbed_is_decoupled` is true, otherwise the
 *  reference state is decoupled. The stage will be called 'decouple'.
 */
void LambdaSchedule::addDecoupleStage(bool perturbed_is_decoupled)
{
    this->addDecoupleStage("decouple", perturbed_is_decoupled);
}

/** Add a named stage to the schedule that will decouple the perturbed
 *  state if `perturbed_is_decoupled` is true, otherwise the
 *  reference state is decoupled.
 */
void LambdaSchedule::addDecoupleStage(const QString &name, bool perturbed_is_decoupled)
{
    this->addStage(name, default_morph_equation);

    // we now need to ensure that the ghost/ghost parameters are not
    // perturbed
    if (perturbed_is_decoupled)
    {
        this->setEquation(name, "ghost/ghost", "*", this->initial());

        // we also need to scale down kappa as the decoupled state is
        // not evaluated in the NonbondedForce, so must not be cancelled
        this->setEquation(name, "ghost/ghost", "kappa", 1.0 - this->lam());
    }
    else
    {
        this->setEquation(name, "ghost/ghost", "*", this->final());

        // we also need to scale up kappa as the decoupled state is
        // not evaluated in the NonbondedForce, so must not be cancelled
        this->setEquation(name, "ghost/ghost", "kappa", this->lam());
    }
}

/** Add a stage to the schedule that will annihilate the perturbed
 *  state if `perturbed_is_annihilated` is true, otherwise the
 *  reference state is annihilated. The stage will be called 'annihilate'.
 */
void LambdaSchedule::addAnnihilateStage(bool perturbed_is_annihilated)
{
    this->addAnnihilateStage("annihilate", perturbed_is_annihilated);
}

/** Add a named stage to the schedule that will annihilate the perturbed
 *  state if `perturbed_is_annihilated` is true, otherwise the
 *  reference state is annihilated.
 */
void LambdaSchedule::addAnnihilateStage(const QString &name, bool perturbed_is_annihilated)
{
    this->addStage(name, default_morph_equation);
}

/** Sandwich the current set of stages with a charge-descaling and
 *  a charge-scaling stage. This prepends a charge-descaling stage
 *  that scales the charge parameter down from `initial` to
 *  :gamma:.initial (where :gamma:=`scale`). The charge parameter in all of
 *  the exising stages in this schedule are then multiplied
 *  by :gamma:. A final charge-rescaling stage is then appended that
 *  scales the charge parameter from :gamma:.final to final.
 */
void LambdaSchedule::addChargeScaleStages(const QString &decharge_name,
                                          const QString &recharge_name,
                                          double scale)
{
    auto scl = this->setConstant("γ", scale);

    // make sure all of the existing stages for the charge lever are scaled
    for (int i = 0; i < this->stage_names.count(); ++i)
    {
        this->setEquation(this->stage_names[i], "*", "charge",
                          scale * this->stage_equations[i].value("charge", this->default_equations[i]));
    }

    // now prepend the decharging stage, and append the recharging stage
    this->prependStage(decharge_name, this->initial());
    this->appendStage(recharge_name, this->final());

    this->setEquation(decharge_name, "*", "charge", (1.0 - ((1.0 - scl) * this->lam())) * this->initial());
    this->setEquation(recharge_name, "*", "charge", (1.0 - ((1.0 - scl) * (1.0 - this->lam()))) * this->final());
}

/** Sandwich the current set of stages with a charge-descaling and
 *  a charge-scaling stage. This prepends a charge-descaling stage
 *  that scales the charge parameter down from `initial` to
 *  :gamma:.initial (where :gamma:=`scale`). The charge parameter in all of
 *  the exising stages in this schedule are then multiplied
 *  by :gamma:. A final charge-rescaling stage is then appended that
 *  scales the charge parameter from :gamma:.final to final.
 */
void LambdaSchedule::addChargeScaleStages(double scale)
{
    this->addChargeScaleStages("decharge", "recharge", scale);
}

/** Prepend a stage called 'name' which uses the passed 'equation'
 *  to the start of this schedule. The equation will be the default
 *  equation that scales all parameters (levers) that don't have
 *  a custom lever for this stage.
 */
void LambdaSchedule::prependStage(const QString &name,
                                  const SireCAS::Expression &equation)
{
    if (name == "*")
        throw SireError::invalid_key(QObject::tr(
                                         "The stage name '*' is reserved and cannot be used."),
                                     CODELOC);

    auto e = equation;

    if (e == default_morph_equation)
        e = default_morph_equation;

    if (this->nStages() == 0)
    {
        this->appendStage(name, e);
        return;
    }

    if (this->stage_names.contains(name))
        throw SireError::invalid_key(QObject::tr(
                                         "Cannot prepend the stage %1 as it already exists.")
                                         .arg(name),
                                     CODELOC);

    this->stage_names.prepend(name);
    this->default_equations.prepend(e);
    this->stage_equations.prepend(QHash<QString, Expression>());
}

/** Append a stage called 'name' which uses the passed 'equation'
 *  to the end of this schedule. The equation will be the default
 *  equation that scales all parameters (levers) that don't have
 *  a custom lever for this stage.
 */
void LambdaSchedule::appendStage(const QString &name,
                                 const SireCAS::Expression &equation)
{
    if (name == "*")
        throw SireError::invalid_key(QObject::tr(
                                         "The stage name '*' is reserved and cannot be used."),
                                     CODELOC);

    else if (this->stage_names.contains(name))
        throw SireError::invalid_key(QObject::tr(
                                         "Cannot append the stage %1 as it already exists.")
                                         .arg(name),
                                     CODELOC);

    auto e = equation;

    if (e == default_morph_equation)
        e = default_morph_equation;

    this->stage_names.append(name);
    this->default_equations.append(e);
    this->stage_equations.append(QHash<QString, Expression>());
}

/** Insert a stage called 'name' at position `i` which uses the passed
 *  'equation'. The equation will be the default
 *  equation that scales all parameters (levers) that don't have
 *  a custom lever for this stage.
 */
void LambdaSchedule::insertStage(int i,
                                 const QString &name,
                                 const SireCAS::Expression &equation)
{
    if (name == "*")
        throw SireError::invalid_key(QObject::tr(
                                         "The stage name '*' is reserved and cannot be used."),
                                     CODELOC);

    auto e = equation;

    if (e == default_morph_equation)
        e = default_morph_equation;

    if (i == 0)
    {
        this->prependStage(name, e);
        return;
    }
    else if (i >= this->nStages())
    {
        this->appendStage(name, e);
        return;
    }

    if (this->stage_names.contains(name))
        throw SireError::invalid_key(QObject::tr(
                                         "Cannot append the stage %1 as it already exists.")
                                         .arg(name),
                                     CODELOC);

    this->stage_names.insert(i, name);
    this->default_equations.insert(i, e);
    this->stage_equations.insert(i, QHash<QString, Expression>());
}

/** Remove the stage 'stage' */
void LambdaSchedule::removeStage(const QString &stage)
{
    if (not this->stage_names.contains(stage))
        return;

    int idx = this->stage_names.indexOf(stage);

    this->stage_names.removeAt(idx);
    this->default_equations.removeAt(idx);
    this->stage_equations.removeAt(idx);
}

/** Append a stage called 'name' which uses the passed 'equation'
 *  to the end of this schedule. The equation will be the default
 *  equation that scales all parameters (levers) that don't have
 *  a custom lever for this stage.
 */
void LambdaSchedule::addStage(const QString &name,
                              const Expression &equation)
{
    if (name == "*")
        throw SireError::invalid_key(QObject::tr(
                                         "The stage name '*' is reserved and cannot be used."),
                                     CODELOC);

    this->appendStage(name, equation);
}

/** Find the index of the stage called 'stage'. This returns
 *  the order in which this stage will take place in the
 *  schedule.
 */
int LambdaSchedule::find_stage(const QString &stage) const
{
    int idx = this->stage_names.indexOf(stage);

    if (idx < 0)
        throw SireError::invalid_key(QObject::tr(
                                         "There is no stage name called '%1'. Valid stages are %2.")
                                         .arg(stage)
                                         .arg(this->stage_names.join(", ")),
                                     CODELOC);

    return idx;
}

/** Set the default equation used to control levers for the
 *  stage 'stage' to 'equation'. This equation will be used
 *  to control any levers in this stage that don't have
 *  their own custom equation.
 */
void LambdaSchedule::setDefaultStageEquation(const QString &stage,
                                             const Expression &equation)
{
    auto e = equation;

    if (e == default_morph_equation)
        e = default_morph_equation;

    this->default_equations[this->find_stage(stage)] = e;
}

/** Set the custom equation used to control the specified 'lever'
 *  for the specified 'force' at the stage 'stage' to 'equation'.
 *  This equation will only be used to control the parameters for the
 *  specified lever in the specified force at the specified stage
 */
void LambdaSchedule::setEquation(const QString &stage,
                                 const QString &force,
                                 const QString &lever,
                                 const SireCAS::Expression &equation)
{
    if (stage == "*")
    {
        // we do this for all stages
        for (int i = 0; i < this->nStages(); ++i)
        {
            this->setEquation(this->stage_names[i], force, lever, equation);
        }

        return;
    }

    auto e = equation;

    if (e == default_morph_equation)
        e = default_morph_equation;

    auto &lever_expressions = this->stage_equations[this->find_stage(stage)];

    if (lever != "*" and not this->lever_names.contains(lever))
        this->addLever(lever);

    if (force != "*" and not this->force_names.contains(force))
        this->addForce(force);

    lever_expressions[_get_lever_name(force, lever)] = e;
}

/** Remove the custom equation for the specified `lever` in the
 *  specified 'force' at the specified `stage`.
 *  The lever will now use the equation specified for this
 *  lever for this stage, or the default lever for the stage
 *  if this isn't set
 */
void LambdaSchedule::removeEquation(const QString &stage,
                                    const QString &force,
                                    const QString &lever)
{
    if (stage == "*")
    {
        // remove from all stages
        for (int i = 0; i < this->nStages(); ++i)
        {
            this->removeEquation(this->stage_names[i], force, lever);
        }

        return;
    }

    int idx = this->stage_names.indexOf(stage);

    this->stage_equations[idx].remove(_get_lever_name(force, lever));
}

/** Return whether or not the specified 'lever' in the specified 'force'
 *  at the specified 'stage' has a custom equation set for it
 */
bool LambdaSchedule::hasForceSpecificEquation(const QString &stage,
                                              const QString &force,
                                              const QString &lever) const
{
    if (stage == "*")
        throw SireError::invalid_key(QObject::tr(
                                         "The stage name '*' is reserved and cannot be used "
                                         "when querying for force-specific equations."),
                                     CODELOC);

    int idx = this->stage_names.indexOf(stage);

    if (idx < 0)
        throw SireError::invalid_key(QObject::tr(
                                         "There is no stage name called '%1'. Valid stages are %2.")
                                         .arg(stage)
                                         .arg(this->stage_names.join(", ")),
                                     CODELOC);

    if (force == "*")
        return false;
    else
        return this->stage_equations[idx].contains(_get_lever_name(force, lever));
}

SireCAS::Expression LambdaSchedule::_getEquation(int stage,
                                                 const QString &force,
                                                 const QString &lever) const
{
    if (stage < 0 or stage >= this->nStages())
        throw SireError::invalid_key(QObject::tr(
                                         "There is no stage number %1. Valid stages are 0-%2.")
                                         .arg(stage)
                                         .arg(this->nStages() - 1),
                                     CODELOC);

    const auto default_lever = _get_lever_name("*", lever);
    const auto default_force = _get_lever_name(force, "*");
    const auto lever_name = _get_lever_name(force, lever);

    const auto equations = this->stage_equations[stage];

    // search from most specific to least specific
    auto it = equations.find(lever_name);

    if (it != equations.end())
    {
        return it.value();
    }

    it = equations.find(default_force);

    if (it != equations.end())
    {
        return it.value();
    }

    it = equations.find(default_lever);

    if (it != equations.end())
    {
        return it.value();
    }

    // we don't have any match, so return the default equation for this stage
    return this->default_equations[stage];
}

/** Return the equation used to control the specified 'lever'
 *  in the specified 'force' at the specified 'stage'. This will
 *  be a custom equation if that has been set for this lever in this
 *  force, or else it will be a custom equation set for this lever,
 *  else it will be the default equation for this stage
 */
Expression LambdaSchedule::getEquation(const QString &stage,
                                       const QString &force,
                                       const QString &lever) const
{
    if (stage == "*")
        throw SireError::invalid_key(QObject::tr(
                                         "The stage name '*' is reserved and cannot be used "
                                         "when getting individual equations."),
                                     CODELOC);

    int idx = this->stage_names.indexOf(stage);

    if (idx < 0)
        throw SireError::invalid_key(QObject::tr(
                                         "There is no stage name called '%1'. Valid stages are %2.")
                                         .arg(stage)
                                         .arg(this->stage_names.join(", ")),
                                     CODELOC);

    return _getEquation(idx, force, lever);
}

/** Set 'schedule' as the molecule-specific schedule for the
 *  perturbable molecule (or part of molecule) that is identified by the
 *  passed 'pert_mol_id'. This schedule will be used to control
 *  all of the levers for this molecule (or part of molecule),
 *  and replaces any levers provided by this schedule
 */
void LambdaSchedule::setMoleculeSchedule(int pert_mol_id,
                                         const LambdaSchedule &schedule)
{
    this->mol_schedules.insert(pert_mol_id, schedule);
    this->mol_schedules[pert_mol_id].mol_schedules.clear();
}

/** Return whether or not the perturbable molecule (or part of molecule)
 *  that is identified by passed 'pert_mol_id' has its own schedule */
bool LambdaSchedule::hasMoleculeSchedule(int pert_mol_id) const
{
    return this->mol_schedules.contains(pert_mol_id);
}

/** Remove the perturbable molecule-specific schedule associated
 *  with the perturbable molecule (or part of molecule) that is
 *  identified by the passed 'pert_mol_id'.
 */
void LambdaSchedule::removeMoleculeSchedule(int pert_mol_id)
{
    this->mol_schedules.remove(pert_mol_id);
}

/** Remove the perturbable molecule-specific schedule associated
 *  with the perturbable molecule (or part of molecule) that is
 *  identified by the passed 'pert_mol_id'. This returns the
 *  schedule that was removed. If no such schedule exists, then
 *  a copy of this schedule is returned.
 */
LambdaSchedule LambdaSchedule::takeMoleculeSchedule(int pert_mol_id)
{
    if (this->mol_schedules.contains(pert_mol_id))
    {
        return this->mol_schedules.take(pert_mol_id);
    }
    else
    {
        auto ret = *this;
        ret.mol_schedules.clear();
        return ret;
    }
}

/** Return the schedule used to control perturbations for the
 *  perturbable molecule (or part of molecule) that is identified by the
 *  passed 'pert_mol_id'. This schedule will be used to control
 *  all of the levers for this molecule (or part of molecule).
 *
 *  This returns this schedule if there is no specified schedule
 *  for this molecule
 */
const LambdaSchedule &LambdaSchedule::getMoleculeSchedule(int pert_mol_id) const
{
    auto it = this->mol_schedules.constFind(pert_mol_id);

    if (it == this->mol_schedules.constEnd())
        return *this;
    else
        return it.value();
}

QVector<double> generate_lambdas(int num_values)
{
    if (num_values < 1)
    {
        return QVector<double>(1, 0.0);
    }
    else if (num_values > 10000)
    {
        num_values = 10000;
    }

    QVector<double> lambda_values(num_values);

    double width = 1.0 / (num_values - 1);

    lambda_values[0] = 0.0;
    lambda_values[num_values - 1] = 1.0;

    for (int i = 1; i < num_values - 1; ++i)
    {
        lambda_values[i] = i * width;
    }

    return lambda_values;
}

/** Return the list of stages that are used for the passed list
 *  of lambda values. The stage names will be returned in the matching
 *  order of the lambda values.
 */
QStringList LambdaSchedule::getLeverStages(const QVector<double> &lambda_values) const
{
    QStringList stages;

    if (this->nStages() == 0)
    {
        for (int i = 0; i < lambda_values.count(); ++i)
        {
            stages.append(QString());
        }

        return stages;
    }

    for (const auto &lambda_value : lambda_values)
    {
        stages.append(this->stage_names[std::get<0>(this->resolve_lambda(lambda_value))]);
    }

    return stages;
}

/** Return the stages used for the list of `nvalue` lambda values
 *  generated for the global lambda value between 0 and 1 inclusive.
 */
QStringList LambdaSchedule::getLeverStages(int nvalues) const
{
    return this->getLeverStages(generate_lambdas(nvalues));
}

/** Return the stage name and parameter values for that lever
 *  for the specified list of lambda values, assuming that a
 *  parameter for that stage has an initial value of
 *  `initial_value` and a final value of `final_value`. This
 *  is mostly useful for testing and graphing how this
 *  schedule would change some hyperthetical forcefield
 *  parameters for the specified lambda values.
 */
QHash<QString, QVector<double>> LambdaSchedule::getLeverValues(
    const QVector<double> &lambda_values,
    double initial_value, double final_value) const
{
    QVector<double> values(lambda_values.count(), NAN);

    QHash<QString, QVector<double>> lever_values;

    // get all of the lever / force combinations in use
    QSet<QString> all_levers;

    for (const auto &equations : this->stage_equations)
    {
        for (const auto &lever : equations.keys())
        {
            all_levers.insert(lever);
        }
    }

    QStringList levers = all_levers.values();
    std::sort(levers.begin(), levers.end());

    lever_values.reserve(levers.count() + 2);

    lever_values.insert("λ", values);

    lever_values.insert("default", values);

    for (const auto &lever : levers)
    {
        if (lever.startsWith("*::"))
            lever_values.insert(lever.mid(3), values);
        else
            lever_values.insert(lever, values);
    }

    if (this->nStages() == 0)
        return lever_values;

    Values input_values = this->constant_values;
    input_values.set(this->initial(), initial_value);
    input_values.set(this->final(), final_value);
    input_values.set(this->lam(), 0.0);

    for (int i = 0; i < lambda_values.count(); ++i)
    {
        const auto lambda_value = lambda_values[i];

        lever_values["λ"][i] = lambda_value;

        const auto resolved = this->resolve_lambda(lambda_value);
        const int stage = std::get<0>(resolved);
        input_values.set(this->lam(), std::get<1>(resolved));

        const auto equation = this->default_equations[stage];
        lever_values["default"][i] = equation(input_values);

        for (const auto &lever : levers)
        {
            auto parts = lever.split("::");

            const auto equation = this->_getEquation(stage, parts[0], parts[1]);

            if (lever.startsWith("*::"))
                lever_values[lever.mid(3)][i] = equation(input_values);
            else
                lever_values[lever][i] = equation(input_values);
        }
    }

    return lever_values;
}

/** Return the lever name and parameter values for that lever
 *  for the specified number of lambda values generated
 *  evenly between 0 and 1, assuming that a
 *  parameter for that lever has an initial value of
 *  `initial_value` and a final value of `final_value`. This
 *  is mostly useful for testing and graphing how this
 *  schedule would change some hyperthetical forcefield
 *  parameters for the specified lambda values.
 */
QHash<QString, QVector<double>> LambdaSchedule::getLeverValues(
    int nvalues,
    double initial_value, double final_value) const
{
    return this->getLeverValues(generate_lambdas(nvalues),
                                initial_value, final_value);
}

/** Return the parameters for the specified lever called `lever_name`
 *  in the force 'force'
 *  that have been morphed from the passed list of initial values
 *  (in `initial`) to the passed list of final values (in `final`)
 *  for the specified global value of :lambda: (in `lambda_value`).
 *
 *  The morphed parameters will be returned in the matching
 *  order to `initial` and `final`.
 *
 *  This morphs a single floating point parameters.
 */
double LambdaSchedule::morph(const QString &force,
                             const QString &lever,
                             double initial, double final,
                             double lambda_value) const
{
    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    Values input_values = this->constant_values;
    input_values.set(this->lam(), std::get<1>(resolved));

    input_values.set(this->initial(), initial);
    input_values.set(this->final(), final);

    return this->_getEquation(stage, force, lever)(input_values);
}

/** Return the parameters for the specified lever called `lever_name`
 *  in the specified force,
 *  that have been morphed from the passed list of initial values
 *  (in `initial`) to the passed list of final values (in `final`)
 *  for the specified global value of :lambda: (in `lambda_value`).
 *
 *  The morphed parameters will be returned in the matching
 *  order to `initial` and `final`.
 *
 *  This morphs floating point parameters. There is an overload
 *  of this function that morphs integer parameters, in which
 *  case the result would be rounded to the nearest integer.
 */
QVector<double> LambdaSchedule::morph(const QString &force,
                                      const QString &lever,
                                      const QVector<double> &initial,
                                      const QVector<double> &final,
                                      double lambda_value) const
{
    const int nparams = initial.count();

    if (final.count() != nparams)
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of initial and final parameters for lever %1 is not the same. "
                                                "%2 versus %3. They need to be the same.")
                                                .arg(lever)
                                                .arg(initial.count())
                                                .arg(final.count()),
                                            CODELOC);

    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto equation = this->_getEquation(stage, force, lever);

    QVector<double> morphed(nparams);
    auto morphed_data = morphed.data();
    const auto initial_data = initial.constData();
    const auto final_data = final.constData();

    if (equation == default_morph_equation)
    {
        for (int i = 0; i < nparams; ++i)
        {
            morphed_data[i] = (1.0 - lambda_value) * initial_data[i] +
                              lambda_value * final_data[i];
        }
    }
    else
    {
        Values input_values = this->constant_values;
        input_values.set(this->lam(), std::get<1>(resolved));

        for (int i = 0; i < nparams; ++i)
        {
            input_values.set(this->initial(), initial_data[i]);
            input_values.set(this->final(), final_data[i]);

            morphed_data[i] = equation(input_values);
        }
    }

    return morphed;
}

/** Return the parameters for the specified lever called `lever_name`
 *  for the specified 'force'
 *  that have been morphed from the passed list of initial values
 *  (in `initial`) to the passed list of final values (in `final`)
 *  for the specified global value of :lambda: (in `lambda_value`).
 *
 *  The morphed parameters will be returned in the matching
 *  order to `initial` and `final`.
 *
 *  This function morphs integer parameters. In this case,
 *  the result will be the rounded to the nearest integer.
 */
QVector<int> LambdaSchedule::morph(const QString &force,
                                   const QString &lever,
                                   const QVector<int> &initial,
                                   const QVector<int> &final,
                                   double lambda_value) const
{
    const int nparams = initial.count();

    if (final.count() != nparams)
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of initial and final parameters for lever %1 is not the same. "
                                                "%2 versus %3. They need to be the same.")
                                                .arg(lever)
                                                .arg(initial.count())
                                                .arg(final.count()),
                                            CODELOC);

    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto equation = this->_getEquation(stage, force, lever);

    QVector<int> morphed(nparams);

    auto morphed_data = morphed.data();
    const auto initial_data = initial.constData();
    const auto final_data = final.constData();

    if (equation == default_morph_equation)
    {
        for (int i = 0; i < nparams; ++i)
        {
            morphed_data[i] = int((1.0 - lambda_value) * initial_data[i] +
                                  lambda_value * final_data[i]);
        }
    }
    else
    {
        Values input_values = this->constant_values;
        input_values.set(this->lam(), std::get<1>(resolved));

        for (int i = 0; i < nparams; ++i)
        {
            input_values.set(this->initial(), double(initial_data[i]));
            input_values.set(this->final(), double(final_data[i]));

            morphed_data[i] = int(equation(input_values));
        }
    }

    return morphed;
}
