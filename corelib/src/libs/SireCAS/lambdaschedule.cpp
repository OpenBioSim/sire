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
    writeHeader(ds, r_schedule, 1);

    SharedDataStream sds(ds);

    sds << schedule.lever_names << schedule.stage_names
        << schedule.default_equations
        << schedule.stage_equations
        << static_cast<const Property &>(schedule);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, LambdaSchedule &schedule)
{
    VersionID v = readHeader(ds, r_schedule);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> schedule.lever_names >> schedule.stage_names >>
            schedule.default_equations >> schedule.stage_equations >>
            static_cast<Property &>(schedule);
    }
    else
        throw version_error(v, "1", r_schedule, CODELOC);

    return ds;
}

LambdaSchedule::LambdaSchedule() : ConcreteProperty<LambdaSchedule, Property>()
{
}

LambdaSchedule::LambdaSchedule(const LambdaSchedule &other)
    : ConcreteProperty<LambdaSchedule, Property>(other),
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
    return lever_names == other.lever_names and
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

    return QObject::tr("LambdaSchedule( num_stages=%1 num_levers=%2 )")
        .arg(this->nStages())
        .arg(this->nLevers());
}

Symbol LambdaSchedule::lambda_symbol("λ");
Symbol LambdaSchedule::initial_symbol("initial");
Symbol LambdaSchedule::final_symbol("final");

Symbol LambdaSchedule::lam()
{
    return lambda_symbol;
}

Symbol LambdaSchedule::initial()
{
    return initial_symbol;
}

Symbol LambdaSchedule::final()
{
    return final_symbol;
}

void LambdaSchedule::addLever(const QString &lever)
{
    if (this->lever_names.contains(lever))
        return;

    this->lever_names.append(lever);
}

void LambdaSchedule::addLevers(const QStringList &levers)
{
    for (const auto &lever : levers)
    {
        if (not this->lever_names.contains(lever))
            this->lever_names.append(lever);
    }
}

void LambdaSchedule::removeLever(const QString &lever)
{
    if (not this->lever_names.contains(lever))
        return;

    int idx = this->lever_names.indexOf(lever);

    this->default_equations.removeAt(idx);
    this->stage_equations.removeAt(idx);
    this->stage_names.removeAt(idx);
}

void LambdaSchedule::removeLevers(const QStringList &levers)
{
    for (const auto &lever : levers)
    {
        this->removeLever(lever);
    }
}

int LambdaSchedule::nLevers() const
{
    return this->lever_names.count();
}

QStringList LambdaSchedule::getLevers() const
{
    return this->lever_names;
}

int LambdaSchedule::nStages() const
{
    return this->stage_names.count();
}

QStringList LambdaSchedule::getStages() const
{
    return this->stage_names;
}

double LambdaSchedule::clamp(double lambda_value) const
{
    return std::max(0.0, std::min(lambda_value, 1.0));
}

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

QString LambdaSchedule::getStage(double lambda_value) const
{
    auto resolved = this->resolve_lambda(lambda_value);

    return this->stage_names[std::get<0>(resolved)];
}

double LambdaSchedule::getLambdaInStage(double lambda_value) const
{
    if (this->nStages() == 0)
        return this->clamp(lambda_value);

    auto resolved = this->resolve_lambda(lambda_value);

    return std::get<1>(resolved);
}

void LambdaSchedule::clear()
{
    this->stage_names.clear();
    this->stage_equations.clear();
    this->default_equations.clear();
}

void LambdaSchedule::addStage(const QString &name,
                              const Expression &default_equation)
{
    if (this->stage_names.contains(name))
        return;

    this->stage_names.append(name);
    this->default_equations.append(default_equation);
    this->stage_equations.append(QHash<QString, Expression>());
}

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

void LambdaSchedule::setDefaultEquation(const QString &stage,
                                        const Expression &equation)
{
    this->default_equations[this->find_stage(stage)] = equation;
}

void LambdaSchedule::setEquation(const QString &stage,
                                 const QString &lever,
                                 const Expression &equation)
{
    if (not this->lever_names.contains(lever))
        throw SireError::invalid_key(QObject::tr(
                                         "There is no lever called '%1'. Valid levers are [ %2 ]")
                                         .arg(lever)
                                         .arg(this->lever_names.join(", ")),
                                     CODELOC);

    auto &lever_expressions = this->stage_equations[this->find_stage(stage)];
    lever_expressions[lever] = equation;
}

void LambdaSchedule::removeEquation(const QString &stage,
                                    const QString &lever)
{
    if (not(this->lever_names.contains(lever) and this->stage_names.contains(stage)))
        return;

    int idx = this->stage_names.indexOf(stage);

    this->stage_equations[idx].remove(lever);
}

Expression LambdaSchedule::getEquation(const QString &stage,
                                       const QString &lever) const
{
    if (not this->lever_names.contains(lever))
        throw SireError::invalid_key(QObject::tr(
                                         "There is no lever called '%1'. Valid levers are [ %2 ]")
                                         .arg(lever)
                                         .arg(this->lever_names.join(", ")),
                                     CODELOC);

    const int idx = this->find_stage(stage);

    const auto &lever_expressions = this->stage_equations[idx];

    return lever_expressions.value(lever, this->default_equations[idx]);
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

QStringList LambdaSchedule::getLeverStages(int nvalues) const
{
    return this->getLeverStages(generate_lambdas(nvalues));
}

QHash<QString, QVector<double>> LambdaSchedule::getLeverValues(
    const QVector<double> &lambda_values,
    double initial_value, double final_value) const
{
    QVector<double> values(lambda_values.count(), NAN);

    QHash<QString, QVector<double>> lever_values;
    lever_values.reserve(this->lever_names.count() + 1);

    lever_values.insert("λ", values);

    for (const auto &stage_name : this->stage_names)
    {
        lever_values.insert(QString("λ_{%1}").arg(stage_name), values);
    }

    for (const auto &lever_name : this->lever_names)
    {
        lever_values.insert(lever_name, values);
    }

    if (this->nStages() == 0)
        return lever_values;

    Values input_values;
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

        lever_values[QString("λ_{%1}").arg(this->stage_names[stage])][i] = std::get<1>(resolved);

        for (const auto &lever_name : lever_names)
        {
            const auto equation = this->stage_equations[stage].value(
                lever_name, this->default_equations[stage]);

            lever_values[lever_name][i] = equation(input_values);
        }
    }

    return lever_values;
}

QHash<QString, QVector<double>> LambdaSchedule::getLeverValues(
    int nvalues,
    double initial_value, double final_value) const
{
    return this->getLeverValues(generate_lambdas(nvalues),
                                initial_value, final_value);
}

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

    Values input_values;
    input_values.set(this->lam(), std::get<1>(resolved));

    QVector<double> morphed(nparams);

    auto morphed_data = morphed.data();
    const auto initial_data = initial.constData();
    const auto final_data = final.constData();

    for (int i = 0; i < nparams; ++i)
    {
        input_values.set(this->initial(), initial_data[i]);
        input_values.set(this->final(), final_data[i]);

        morphed_data[i] = equation(input_values);
    }

    return morphed;
}

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

    Values input_values;
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
