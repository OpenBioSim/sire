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

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;

static RegisterMetaType<LambdaSchedule> r_schedule;

QDataStream &operator<<(QDataStream &ds, const LambdaSchedule &schedule)
{
    writeHeader(ds, r_schedule, 2);

    SharedDataStream sds(ds);

    sds << schedule.constant_values
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

    if (v == 1 or v == 2)
    {
        SharedDataStream sds(ds);

        sds >> schedule.constant_values >>
            schedule.lever_names >> schedule.stage_names >>
            schedule.default_equations >> schedule.stage_equations;

        if (v == 2)
            sds >> schedule.mol_schedules;

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
        throw version_error(v, "1, 2", r_schedule, CODELOC);

    return ds;
}

LambdaSchedule::LambdaSchedule() : ConcreteProperty<LambdaSchedule, Property>()
{
}

LambdaSchedule::LambdaSchedule(const LambdaSchedule &other)
    : ConcreteProperty<LambdaSchedule, Property>(other),
      mol_schedules(other.mol_schedules),
      constant_values(other.constant_values),
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

        for (const auto &lever : this->stage_equations[i].keys())
        {
            lines.append(QString("    %1: %2")
                             .arg(lever)
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
    if (this->lever_names.contains(lever))
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
        if (not this->lever_names.contains(lever))
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
        this->setEquation(this->stage_names[i], "charge",
                          scl * this->stage_equations[i].value("charge", this->default_equations[i]));
    }

    // now prepend the decharging stage, and append the recharging stage
    this->prependStage(decharge_name, this->initial());
    this->appendStage(recharge_name, this->final());

    this->setEquation(decharge_name, "charge", (1.0 - ((1.0 - scl) * this->lam())) * this->initial());
    this->setEquation(recharge_name, "charge", (1.0 - ((1.0 - scl) * (1.0 - this->lam()))) * this->final());
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
    if (this->stage_names.contains(name))
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

/** Append a stage called 'name' which uses the passed 'equation'
 *  to the end of this schedule. The equation will be the default
 *  equation that scales all parameters (levers) that don't have
 *  a custom lever for this stage.
 */
void LambdaSchedule::addStage(const QString &name,
                              const Expression &equation)
{
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
void LambdaSchedule::setDefaultEquation(const QString &stage,
                                        const Expression &equation)
{
    auto e = equation;

    if (e == default_morph_equation)
        e = default_morph_equation;

    this->default_equations[this->find_stage(stage)] = e;
}

/** Set the custom equation used to control the specified
 *  `lever` at the stage `stage` to `equation`. This equation
 *  will only be used to control the parameters for the
 *  specified lever at the specified stage.
 */
void LambdaSchedule::setEquation(const QString &stage,
                                 const QString &lever,
                                 const Expression &equation)
{
    auto e = equation;

    if (e == default_morph_equation)
        e = default_morph_equation;

    auto &lever_expressions = this->stage_equations[this->find_stage(stage)];

    if (not this->lever_names.contains(lever))
        this->addLever(lever);

    lever_expressions[lever] = e;
}

QString _create_lever_name(const QString &force, const QString &lever)
{
    return force + "::" + lever;
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
    this->setEquation(stage, _create_lever_name(force, lever), equation);
}

/** Remove the custom equation for the specified `lever` at the
 *  specified `stage`. The lever will now use the default
 *  equation at this stage.
 */
void LambdaSchedule::removeEquation(const QString &stage,
                                    const QString &lever)
{
    if (not(this->lever_names.contains(lever) and this->stage_names.contains(stage)))
        return;

    int idx = this->stage_names.indexOf(stage);

    this->stage_equations[idx].remove(lever);
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
    this->removeEquation(stage, _create_lever_name(force, lever));
}

/** Return the default equation used to control the parameters for
 *  the stage `stage`.
 */
Expression LambdaSchedule::getEquation(const QString &stage) const
{
    const int idx = this->find_stage(stage);

    return this->default_equations[idx];
}

/** Return the equation used to control the specified `lever`
 *  at the specified `stage`. This will be a custom equation
 *  if that has been set for this lever, or else the
 *  default equation for this stage.
 */
Expression LambdaSchedule::getEquation(const QString &stage,
                                       const QString &lever) const
{
    if (not this->lever_names.contains(lever))
        return this->getEquation(stage);

    const int idx = this->find_stage(stage);

    const auto &lever_expressions = this->stage_equations[idx];

    return lever_expressions.value(lever, this->default_equations[idx]);
}

/** Return whether the force 'force' has a force-specific equation
 *  for the specified 'lever' at the specified 'stage'
 */
bool LambdaSchedule::hasForceSpecificEquation(const QString &stage,
                                              const QString &force,
                                              const QString &lever) const
{
    const auto force_lever = _create_lever_name(force, lever);

    if (not this->lever_names.contains(force_lever))
        return false;

    const int idx = this->find_stage(stage);

    const auto &lever_expressions = this->stage_equations[idx];

    return lever_expressions.contains(force_lever);
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
    if (this->hasForceSpecificEquation(stage, force, lever))
    {
        return this->getEquation(stage, _create_lever_name(force, lever));
    }
    else
    {
        return this->getEquation(stage, lever);
    }
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

/** Return the list of lever stages that are used for the passed list
 *  of lambda values. The lever names will be returned in the matching
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

/** Return the lever stages used for the list of `nvalue` lambda values
 *  generated for the global lambda value between 0 and 1 inclusive.
 */
QStringList LambdaSchedule::getLeverStages(int nvalues) const
{
    return this->getLeverStages(generate_lambdas(nvalues));
}

/** Return the lever name and parameter values for that lever
 *  for the specified list of lambda values, assuming that a
 *  parameter for that lever has an initial value of
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
    lever_values.reserve(this->lever_names.count() + 1);

    lever_values.insert("λ", values);

    lever_values.insert("default", values);

    for (const auto &lever_name : this->lever_names)
    {
        lever_values.insert(lever_name, values);
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

        for (const auto &lever_name : lever_names)
        {
            const auto equation = this->stage_equations[stage].value(
                lever_name, this->default_equations[stage]);

            lever_values[lever_name][i] = equation(input_values);
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
 *  that have been morphed from the passed list of initial values
 *  (in `initial`) to the passed list of final values (in `final`)
 *  for the specified global value of :lambda: (in `lambda_value`).
 *
 *  The morphed parameters will be returned in the matching
 *  order to `initial` and `final`.
 *
 *  This morphs a single floating point parameters.
 */
double LambdaSchedule::morph(const QString &lever_name,
                             double initial, double final,
                             double lambda_value) const
{
    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto equation = this->stage_equations[stage].value(
        lever_name, this->default_equations[stage]);

    Values input_values = this->constant_values;
    input_values.set(this->lam(), std::get<1>(resolved));

    input_values.set(this->initial(), initial);
    input_values.set(this->final(), final);

    return equation(input_values);
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
                             const QString &lever_name,
                             double initial, double final,
                             double lambda_value) const
{
    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto force_lever = _create_lever_name(force, lever_name);

    const auto equation = this->stage_equations[stage].value(
        force_lever,
        this->stage_equations[stage].value(
            lever_name, this->default_equations[stage]));

    Values input_values = this->constant_values;
    input_values.set(this->lam(), std::get<1>(resolved));

    input_values.set(this->initial(), initial);
    input_values.set(this->final(), final);

    return equation(input_values);
}

/** Return the parameters for the specified lever called `lever_name`
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
QVector<double> LambdaSchedule::morph(const QString &lever_name,
                                      const QVector<double> &initial,
                                      const QVector<double> &final,
                                      double lambda_value) const
{
    const int nparams = initial.count();

    if (final.count() != nparams)
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of initial and final parameters for lever %1 is not the same. "
                                                "%2 versus %3. They need to be the same.")
                                                .arg(lever_name)
                                                .arg(initial.count())
                                                .arg(final.count()),
                                            CODELOC);

    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto equation = this->stage_equations[stage].value(
        lever_name, this->default_equations[stage]);

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
                                      const QString &lever_name,
                                      const QVector<double> &initial,
                                      const QVector<double> &final,
                                      double lambda_value) const
{
    const int nparams = initial.count();

    if (final.count() != nparams)
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of initial and final parameters for lever %1 is not the same. "
                                                "%2 versus %3. They need to be the same.")
                                                .arg(lever_name)
                                                .arg(initial.count())
                                                .arg(final.count()),
                                            CODELOC);

    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto force_lever = _create_lever_name(force, lever_name);

    const auto equation = this->stage_equations[stage].value(
        force_lever, this->stage_equations[stage].value(
                         lever_name, this->default_equations[stage]));

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
QVector<int> LambdaSchedule::morph(const QString &lever_name,
                                   const QVector<int> &initial,
                                   const QVector<int> &final,
                                   double lambda_value) const
{
    const int nparams = initial.count();

    if (final.count() != nparams)
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of initial and final parameters for lever %1 is not the same. "
                                                "%2 versus %3. They need to be the same.")
                                                .arg(lever_name)
                                                .arg(initial.count())
                                                .arg(final.count()),
                                            CODELOC);

    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto equation = this->stage_equations[stage].value(
        lever_name, this->default_equations[stage]);

    Values input_values = this->constant_values;
    input_values.set(this->lam(), std::get<1>(resolved));

    QVector<int> morphed(nparams);

    auto morphed_data = morphed.data();
    const auto initial_data = initial.constData();
    const auto final_data = final.constData();

    for (int i = 0; i < nparams; ++i)
    {
        input_values.set(this->initial(), double(initial_data[i]));
        input_values.set(this->final(), double(final_data[i]));

        // the result is the resulting float rounded to the nearest
        // integer
        morphed_data[i] = int(std::floor(equation(input_values) + 0.5));
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
                                   const QString &lever_name,
                                   const QVector<int> &initial,
                                   const QVector<int> &final,
                                   double lambda_value) const
{
    const int nparams = initial.count();

    if (final.count() != nparams)
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of initial and final parameters for lever %1 is not the same. "
                                                "%2 versus %3. They need to be the same.")
                                                .arg(lever_name)
                                                .arg(initial.count())
                                                .arg(final.count()),
                                            CODELOC);

    if (this->nStages() == 0)
        // just return the initial parameters as we don't know how to morph
        return initial;

    auto force_lever = _create_lever_name(force, lever_name);

    const auto resolved = this->resolve_lambda(lambda_value);
    const int stage = std::get<0>(resolved);

    const auto equation = this->stage_equations[stage].value(
        force_lever, this->stage_equations[stage].value(
                         lever_name, this->default_equations[stage]));

    Values input_values = this->constant_values;
    input_values.set(this->lam(), std::get<1>(resolved));

    QVector<int> morphed(nparams);

    auto morphed_data = morphed.data();
    const auto initial_data = initial.constData();
    const auto final_data = final.constData();

    for (int i = 0; i < nparams; ++i)
    {
        input_values.set(this->initial(), double(initial_data[i]));
        input_values.set(this->final(), double(final_data[i]));

        // the result is the resulting float rounded to the nearest
        // integer
        morphed_data[i] = int(std::floor(equation(input_values) + 0.5));
    }

    return morphed;
}
