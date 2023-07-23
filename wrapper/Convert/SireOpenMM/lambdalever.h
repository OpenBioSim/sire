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

#ifndef SIREOPENMM_LAMBDALEVER_H
#define SIREOPENMM_LAMBDALEVER_H

#include "openmmmolecule.h"

#include "SireCAS/expression.h"
#include "SireCAS/symbol.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{

    /** This is a lever that is used to change the parameters in an OpenMM
     *  context according to a lambda value. This is actually a collection
     *  of levers, each of which is controlled by the main lever.
     *
     *  You can use SireCAS expressions to control how each lever changes
     *  each parameter
     */
    class LambdaLever : public SireBase::ConcreteProperty<LambdaLever, SireBase::Property>
    {
    public:
        LambdaLever();
        LambdaLever(const LambdaLever &other);

        virtual ~LambdaLever();

        LambdaLever &operator=(const LambdaLever &other);

        bool operator==(const LambdaLever &other) const;
        bool operator!=(const LambdaLever &other) const;

        LambdaLever *clone() const;

        const char *what() const;
        static const char *typeName();

        void set_lambda(OpenMM::Context &context, double lam_val) const;

        void set_force_index(const QString &force, int index);

        void add_perturbable_molecule(const OpenMMMolecule &molecule,
                                      const QHash<QString, qint32> &start_indicies);

        static SireCAS::Symbol lambda();
        static SireCAS::Symbol initial();
        static SireCAS::Symbol final();

        void add_lever(const QString &lever);
        void add_levers(const QStringList &levers);

        int num_levers() const;

        QStringList get_levers() const;

        int num_stages() const;

        QStringList get_stages() const;

        QString get_stage(double lambda_value) const;

        double get_lambda_in_stage(double lambda_value) const;

        void add_stage(const QString &stage,
                       const SireCAS::Expression &equation);

        void clear();

        void set_equation(const QString &stage,
                          const QString &lever,
                          const SireCAS::Expression &equation);

        void set_default_equation(const QString &stage,
                                  const SireCAS::Expression &equation);

        SireCAS::Expression get_equation(const QString &stage,
                                         const QString &lever) const;

        QHash<QString, QVector<double>> get_lever_values(const QVector<double> &lambda_values,
                                                         double initial_value,
                                                         double final_value) const;

        QStringList get_lever_stages(const QVector<double> &lambda_values) const;

    protected:
        int find_stage(const QString &stage) const;
        double clamp(double lambda_value) const;

        std::tuple<int, double> resolve_lambda(double lambda) const;

        /** Map from a forcefield name to its index in the associated System */
        QHash<QString, int> name_to_ffidx;

        /** The names of all of the levers provided by the forcefields */
        QStringList lever_names;

        /** The names of each stage of the lever */
        QStringList stage_names;

        /** The default equations used for each stage */
        QVector<SireCAS::Expression> default_equations;

        /** Any over-ridden equations for a particular lever in a
            particular stage */
        QVector<QHash<QString, SireCAS::Expression>> stage_equations;

        /** The symbol used to represent the lambda value */
        static SireCAS::Symbol lambda_symbol;

        /** The symbol used to represent the initial value */
        static SireCAS::Symbol initial_symbol;

        /** The symbol used to represent the final value */
        static SireCAS::Symbol final_symbol;
    };
}

Q_DECLARE_METATYPE(SireOpenMM::LambdaLever)

SIRE_END_HEADER

#endif
