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
    : SireBase::ConcreteProperty<LambdaLever, SireBase::Property>(other)
{
}

LambdaLever::~LambdaLever()
{
}

LambdaLever &LambdaLever::operator=(const LambdaLever &other)
{
    Property::operator=(other);
    return *this;
}

bool LambdaLever::operator==(const LambdaLever &other) const
{
    return true;
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

void LambdaLever::set_lambda(OpenMM::Context &context,
                             double lambda_value) const
{
    qDebug() << "set lambda to" << lambda_value;
}

void LambdaLever::set_force_index(const QString &force, int index)
{
}

void LambdaLever::add_perturbable_molecule(const OpenMMMolecule &molecule,
                                           const QHash<QString, qint32> &start_indicies)
{
}

Symbol LambdaLever::lambda_symbol("λ");
Symbol LambdaLever::initial_symbol("initial");
Symbol LambdaLever::final_symbol("final");

Symbol LambdaLever::lambda()
{
    return lambda_symbol;
}

Symbol LambdaLever::initial()
{
    return initial_symbol;
}

Symbol LambdaLever::final()
{
    return final_symbol;
}

void LambdaLever::add_lever(const QString &lever)
{
    if (this->lever_names.contains(lever))
        return;

    this->lever_names.append(lever);
}

void LambdaLever::add_levers(const QStringList &levers)
{
    for (const auto &lever : levers)
    {
        if (not this->lever_names.contains(lever))
            this->lever_names.append(lever);
    }
}

int LambdaLever::num_levers() const
{
    return this->lever_names.count();
}

QStringList LambdaLever::get_levers() const
{
    return this->lever_names;
}

int LambdaLever::num_stages() const
{
    return this->stage_names.count();
}

QStringList LambdaLever::get_stages() const
{
    return this->stage_names;
}

double LambdaLever::clamp(double lambda_value) const
{
    return std::max(0.0, std::min(lambda_value, 1.0));
}

std::tuple<int, double> LambdaLever::resolve_lambda(double lambda_value) const
{
    if (this->num_stages() == 0)
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
        return std::tuple<int, double>(this->num_stages() - 1, 1.0);
    }

    double stage_width = 1.0 / this->num_stages();

    double resolved = lambda_value / stage_width;

    double stage = std::floor(resolved);

    return std::tuple<int, double>(int(stage), resolved - stage);
}

QString LambdaLever::get_stage(double lambda_value) const
{
    auto resolved = this->resolve_lambda(lambda_value);

    return this->stage_names[std::get<0>(resolved)];
}

double LambdaLever::get_lambda_in_stage(double lambda_value) const
{
    if (this->num_stages() == 0)
        return this->clamp(lambda_value);

    auto resolved = this->resolve_lambda(lambda_value);

    return std::get<1>(resolved);
}

void LambdaLever::clear()
{
    this->stage_names.clear();
    this->stage_equations.clear();
    this->default_equations.clear();
}

void LambdaLever::add_stage(const QString &name,
                            const Expression &default_equation)
{
    if (this->stage_names.contains(name))
        return;

    this->stage_names.append(name);
    this->default_equations.append(default_equation);
    this->stage_equations.append(QHash<QString, Expression>());
}

int LambdaLever::find_stage(const QString &stage) const
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

void LambdaLever::set_default_equation(const QString &stage,
                                       const Expression &equation)
{
    this->default_equations[this->find_stage(stage)] = equation;
}

void LambdaLever::set_equation(const QString &stage,
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

Expression LambdaLever::get_equation(const QString &stage,
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

QStringList LambdaLever::get_lever_stages(const QVector<double> &lambda_values) const
{
    QStringList stages;

    if (this->num_stages() == 0)
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

QHash<QString, QVector<double>> LambdaLever::get_lever_values(
    const QVector<double> &lambda_values,
    double initial_value, double final_value) const
{
    QVector<double> values(lambda_values.count(), 0.0);

    QHash<QString, QVector<double>> lever_values;
    lever_values.reserve(this->lever_names.count());

    for (const auto &lever_name : this->lever_names)
    {
        lever_values.insert(lever_name, values);
    }

    if (this->num_stages() == 0)
        return lever_values;

    Values input_values;
    input_values.add(this->initial() == initial_value);
    input_values.add(this->final() == final_value);
    input_values.add(this->lambda() == 0.0);

    for (int i = 0; i < lambda_values.count(); ++i)
    {
        const auto lambda_value = lambda_values[i];

        const auto resolved = this->resolve_lambda(lambda_value);
        const int stage = std::get<0>(resolved);
        input_values.set(this->lambda(), std::get<1>(resolved));

        for (const auto &lever_name : lever_names)
        {
            const auto equation = this->stage_equations[stage].value(
                lever_name, this->default_equations[stage]);

            lever_values[lever_name][i] = equation(input_values);
        }
    }

    return lever_values;
}
