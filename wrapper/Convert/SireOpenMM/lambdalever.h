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

#ifndef SIREOPENMM_LAMBDALEVER_H
#define SIREOPENMM_LAMBDALEVER_H

#include "openmmmolecule.h"

#include "SireCAS/lambdaschedule.h"

#include <QReadWriteLock>

#include <memory>

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    class MolLambdaCache
    {
    public:
        MolLambdaCache();
        MolLambdaCache(double lam_val);
        MolLambdaCache(const MolLambdaCache &other);
        ~MolLambdaCache();

        MolLambdaCache &operator=(const MolLambdaCache &other);

        const QVector<double> &morph(const SireCAS::LambdaSchedule &schedule,
                                     const QString &force, const QString &key,
                                     const QVector<double> &initial,
                                     const QVector<double> &final) const;

    private:
        QHash<QString, QVector<double>> cache;
        QReadWriteLock lock;
        double lam_val;
    };

    class LeverCache
    {
    public:
        LeverCache();
        LeverCache(const LeverCache &other);
        ~LeverCache();

        LeverCache &operator=(const LeverCache &other);

        const MolLambdaCache &get(int molidx, double lam_val) const;

        void clear();

    private:
        QHash<int, QHash<double, MolLambdaCache>> cache;
    };

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

        double setLambda(OpenMM::Context &context, double lam_val,
                         bool update_constraints = true) const;

        void setForceIndex(const QString &force, int index);

        void addRestraintIndex(const QString &force, int index);

        int addPerturbableMolecule(const OpenMMMolecule &molecule,
                                   const QHash<QString, qint32> &start_indicies);

        void setExceptionIndicies(int idx, const QString &ff,
                                  const QVector<std::pair<int, int>> &exception_idxs);

        void setConstraintIndicies(int idx, const QVector<qint32> &constraint_idxs);

        void setSchedule(const SireCAS::LambdaSchedule &schedule);

        void addLever(const QString &lever_name);

        bool hasLever(const QString &lever_name);

        QHash<SireMol::MolNum, SireBase::PropertyMap> getPerturbableMoleculeMaps() const;

        SireCAS::LambdaSchedule getSchedule() const;

        template <class T>
        T *getForce(const QString &name, OpenMM::System &system) const;

        QList<OpenMM::Force *> getRestraints(const QString &name,
                                             OpenMM::System &system) const;

        int getForceIndex(const QString &name) const;
        QString getForceType(const QString &name,
                             const OpenMM::System &system) const;

    protected:
        void updateRestraintInContext(OpenMM::Force &ff, double rho,
                                      OpenMM::Context &context) const;

        /** Map from a forcefield name to its index in the associated System */
        QHash<QString, qint32> name_to_ffidx;

        /** Map from a restraint name to its index in the associated System.
         *  Note that multiple restraints can have the same name */
        QMultiHash<QString, qint32> name_to_restraintidx;

        /** The schedule used to set lambda */
        SireCAS::LambdaSchedule lambda_schedule;

        /** The list of perturbable molecules */
        QVector<PerturbableOpenMMMolecule> perturbable_mols;

        /** The start indicies of the parameters in each named
            forcefield for each perturbable moleucle */
        QVector<QHash<QString, qint32>> start_indicies;

        /** All of the property maps for the perturbable molecules */
        QHash<SireMol::MolNum, SireBase::PropertyMap> perturbable_maps;

        /** Cache of the parameters for different lambda values */
        LeverCache lambda_cache;
    };

#ifndef SIRE_SKIP_INLINE_FUNCTION

    template <class T>
    inline QString _get_typename()
    {
        return QString::fromStdString(T().getName());
    }

    template <>
    inline QString _get_typename<OpenMM::CustomNonbondedForce>()
    {
        return "OpenMM::CustomNonbondedForce";
    }

    template <>
    inline QString _get_typename<OpenMM::CustomBondForce>()
    {
        return "OpenMM::CustomBondForce";
    }

    // LESTER - UNCOMMENT BELOW FOR FEATURE_EMLE

    /*template <>
    inline QString _get_typename<SireOpenMM::QMMMForce>()
    {
        return "SireOpenMM::QMMMForce";
    }*/

    /** Return the OpenMM::Force (of type T) that is called 'name'
     *  from the passed OpenMM::System. This returns 0 if the force
     *  doesn't exist
     */
    template <class T>
    T *LambdaLever::getForce(const QString &name, OpenMM::System &system) const
    {
        auto it = name_to_ffidx.constFind(name);

        if (it == name_to_ffidx.constEnd())
        {
            return 0;
        }

        const int idx = it.value();

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
                                              "Cannot cast the force called '%1' to a %2. We think it is a %3.")
                                              .arg(name)
                                              .arg(_get_typename<T>())
                                              .arg(QString::fromStdString(force.getName())),
                                          CODELOC);
        }

        return t_force;
    }

#endif
}

Q_DECLARE_METATYPE(SireOpenMM::LambdaLever)

SIRE_EXPOSE_CLASS(SireOpenMM::LambdaLever)

SIRE_END_HEADER

#endif
