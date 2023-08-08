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

#include "SireCAS/lambdaschedule.h"

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

        double setLambda(OpenMM::Context &context, double lam_val) const;

        void setForceIndex(const QString &force, int index);

        int addPerturbableMolecule(const OpenMMMolecule &molecule,
                                   const QHash<QString, qint32> &start_indicies);

        void setExceptionIndicies(int idx, const QString &ff,
                                  const QVector<int> &exception_idxs);

        void setSchedule(const SireCAS::LambdaSchedule &schedule);

        SireCAS::LambdaSchedule getSchedule() const;

    protected:
        /** Map from a forcefield name to its index in the associated System */
        QHash<QString, int> name_to_ffidx;

        /** The schedule used to set lambda */
        SireCAS::LambdaSchedule lambda_schedule;

        /** The list of perturbable molecules */
        QVector<OpenMMMolecule> perturbable_mols;

        /** The start indicies of the parameters in each named
            forcefield for each perturbable moleucle */
        QVector<QHash<QString, qint32>> start_indicies;
    };
}

Q_DECLARE_METATYPE(SireOpenMM::LambdaLever)

SIRE_END_HEADER

#endif
