/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include "calculate_energy.h"

#include "SireMM/cljshiftfunction.h"
#include "SireMM/cljrffunction.h"
#include "SireMM/interff.h"
#include "SireMM/intergroupff.h"
#include "SireMM/internalff.h"
#include "SireMM/internalgroupff.h"
#include "SireMM/intraff.h"
#include "SireMM/intragroupff.h"
#include "SireMM/mmdetail.h"

#include "SireMM/intragroupff.h"

#include "SireBase/generalunitproperty.h"
#include "SireBase/parallel.h"

#include "SireSystem/forcefieldinfo.h"

#include "SireUnits/units.h"

#include "SireMol/core.h"

using namespace SireSystem;
using namespace SireMM;
using namespace SireFF;
using namespace SireMol;
using namespace SireVol;

using namespace SireBase;

using namespace SireUnits;
using namespace SireUnits::Dimension;

#include <QDebug>

namespace SireSystem
{
    ForceFieldInfo get_ffinfo(const MoleculeView &mol, const PropertyMap &map)
    {
        return ForceFieldInfo(mol, map);
    }

    template <class T>
    T get_cljfunc(const ForceFieldInfo &ffinfo)
    {
        const auto &mm = ffinfo.detail().asA<MMDetail>();

        if (not(mm.isAmberStyle() or mm.isOPLS()))
        {
            throw SireError::incompatible_error(QObject::tr("Calculating energies of forcefields that are not Amber- or "
                                                            "OPLS-style is currently not supported. The forcefield style "
                                                            "is %1.")
                                                    .arg(mm.toString()),
                                                CODELOC);
        }

        T cljfunc(ffinfo.cutoff());

        if (mm.usesArithmeticCombiningRules())
            cljfunc.setArithmeticCombiningRules(true);
        else if (mm.usesGeometricCombiningRules())
            cljfunc.setGeometricCombiningRules(true);

        return cljfunc;
    }

    CLJFunctionPtr get_intra_cljfunc(const ForceFieldInfo &ffinfo)
    {
        const auto cutoff_type = ffinfo.cutoffType();

        if (cutoff_type == "CUTOFF" or cutoff_type == "NO_CUTOFF" or cutoff_type == "SHIFT_ELECTROSTATICS")
        {
            return CLJFunctionPtr(get_cljfunc<CLJIntraShiftFunction>(ffinfo));
        }
        else if (cutoff_type == "REACTION_FIELD")
        {
            auto func = get_cljfunc<CLJIntraRFFunction>(ffinfo);
            func.setDielectric(ffinfo.getParameter("dielectric").value());
            return CLJFunctionPtr(func);
        }
        else
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "Cannot calculate energy as cutoff type %1 is not supported.")
                                                    .arg(cutoff_type),
                                                CODELOC);
        }
    }

    CLJFunctionPtr get_inter_cljfunc(const ForceFieldInfo &ffinfo)
    {
        const auto cutoff_type = ffinfo.cutoffType();

        if (cutoff_type == "CUTOFF" or cutoff_type == "NO_CUTOFF" or cutoff_type == "SHIFT_ELECTROSTATICS")
        {
            return CLJFunctionPtr(get_cljfunc<CLJShiftFunction>(ffinfo));
        }
        else if (cutoff_type == "REACTION_FIELD")
        {
            auto func = get_cljfunc<CLJRFFunction>(ffinfo);
            func.setDielectric(ffinfo.getParameter("dielectric").value());
            return CLJFunctionPtr(func);
        }
        else
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "Cannot calculate energy as cutoff type %1 is not supported.")
                                                    .arg(cutoff_type),
                                                CODELOC);
        }
    }

    SireUnits::Dimension::Length get_cutoff(const ForceFieldInfo &ffinfo)
    {
        return ffinfo.cutoff();
    }

    SIRESYSTEM_EXPORT ForceFields create_forcefield(const MoleculeView &mol, const PropertyMap &map)
    {
        auto ffinfo = get_ffinfo(mol, map);

        const auto &mm = ffinfo.detail().asA<MMDetail>();

        ForceFields ffields;

        InternalFF internalff("internal");
        internalff.setStrict(true);
        internalff.enable14Calculation();

        if (mm.usesArithmeticCombiningRules())
            internalff.setArithmeticCombiningRules(true);
        else if (mm.usesGeometricCombiningRules())
            internalff.setGeometricCombiningRules(true);

        internalff.add(mol, map);

        // internal calculation, so should not use a cutoff?
        ffinfo.setNoCutoff();

        IntraFF intraff("intraff");
        intraff.setCLJFunction(get_intra_cljfunc(ffinfo).read().asA<CLJIntraFunction>());
        intraff.add(mol, map);

        ffields.add(intraff);
        ffields.add(internalff);

        return ffields;
    }

    SIRESYSTEM_EXPORT ForceFields create_forcefield(const SireMol::Molecules &mols, const SireBase::PropertyMap &map)
    {
        if (mols.isEmpty())
            return ForceFields();
        else if (mols.nMolecules() == 1)
            return create_forcefield(mols.first(), map);

        const auto ffinfo = get_ffinfo(*(mols.constBegin()), map);

        auto intra_ffinfo = ffinfo;
        intra_ffinfo.setNoCutoff();

        const auto &mm = ffinfo.detail().asA<MMDetail>();

        ForceFields ffields;

        InternalFF internalff("internal");
        internalff.setStrict(true);
        internalff.enable14Calculation();

        if (mm.usesArithmeticCombiningRules())
            internalff.setArithmeticCombiningRules(true);
        else if (mm.usesGeometricCombiningRules())
            internalff.setGeometricCombiningRules(true);

        internalff.add(mols, map);

        InterFF interff("interff");
        interff.setCLJFunction(get_inter_cljfunc(ffinfo).read());
        interff.setProperty("space", ffinfo.space());
        interff.add(mols, map);

        IntraFF intraff("intraff");
        intraff.setCLJFunction(get_intra_cljfunc(intra_ffinfo).read().asA<CLJIntraFunction>());
        intraff.add(mols, map);

        ffields.add(interff);
        ffields.add(intraff);
        ffields.add(internalff);

        return ffields;
    }

    SIRESYSTEM_EXPORT ForceFields create_forcefield(const MoleculeView &mol0, const MoleculeView &mol1, const PropertyMap &map)
    {
        return create_forcefield(Molecules(mol0), Molecules(mol1), map);
    }

    SIRESYSTEM_EXPORT ForceFields create_forcefield(const MoleculeView &mol0, const Molecules &mols1, const PropertyMap &map)
    {
        return create_forcefield(Molecules(mol0), mols1, map);
    }

    SIRESYSTEM_EXPORT ForceFields create_forcefield(const Molecules &mols0, const Molecules &mols1, const PropertyMap &map)
    {
        if (mols0.isEmpty() or mols1.isEmpty())
            return ForceFields();

        auto ffinfo0 = get_ffinfo(*(mols0.constBegin()), map);
        auto ffinfo1 = get_ffinfo(*(mols1.constBegin()), map);

        const auto &mm0 = ffinfo0.detail().asA<MMDetail>();
        const auto &mm1 = ffinfo1.detail().asA<MMDetail>();

        if (not mm0.isCompatibleWith(mm1))
        {
            throw SireError::incompatible_error(QObject::tr("Cannot calculate the energy as the molecules have different "
                                                            "types of MM forcefield: %1 versus %2")
                                                    .arg(mm0.toString())
                                                    .arg(mm1.toString()),
                                                CODELOC);
        }

        ForceFields ffields;

        InterGroupFF interff("interff");

        // check to see if this is a pure intramolecular energy
        const auto molnums0 = mols0.molNums();
        const auto molnums1 = mols1.molNums();

        const bool is_single_mol = (molnums0.count() == 1 and molnums1.count() == 1 and molnums0 == molnums1);

        if (is_single_mol)
        {
            ffinfo0.setNoCutoff();
        }
        else
        {
            interff.setCLJFunction(get_inter_cljfunc(ffinfo0).read());
            interff.setProperty("space", ffinfo0.space());
            interff.add(mols0, MGIdx(0), map);
            interff.add(mols1, MGIdx(1), map);
        }

        auto intra_ffinfo = ffinfo0;
        intra_ffinfo.setNoCutoff();

        IntraGroupFF intraff("intraff");
        intraff.setCLJFunction(get_intra_cljfunc(intra_ffinfo).read().asA<CLJIntraFunction>());
        intraff.add(mols0, MGIdx(0), map);
        intraff.add(mols1, MGIdx(1), map);

        InternalGroupFF internalff("internalff");
        internalff.enable14Calculation();

        if (mm0.usesArithmeticCombiningRules())
            internalff.setArithmeticCombiningRules(true);
        else if (mm0.usesGeometricCombiningRules())
            internalff.setGeometricCombiningRules(true);

        internalff.add(mols0, MGIdx(0), map);
        internalff.add(mols1, MGIdx(1), map);

        if (not is_single_mol)
            ffields.add(interff);

        ffields.add(intraff);
        ffields.add(internalff);

        return ffields;
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(ForceFields &ffields)
    {
        auto nrgs = ffields.energies();

        double total = ffields.energy().value();

        GeneralUnit nrg;

        for (const auto &ffield : ffields.forceFields())
        {
            if (ffield->isA<InterFF>() or ffield->isA<InterGroupFF>())
            {
                const auto &clj = ffield.read().components().asA<MultiCLJComponent>();
                nrg.addComponent("coulomb", nrgs[clj.coulomb()] * kcal_per_mol);
                nrg.addComponent("LJ", nrgs[clj.lj()] * kcal_per_mol);
            }
            else if (ffield->isA<IntraFF>() or ffield->isA<IntraGroupFF>())
            {
                const auto &clj = ffield.read().components().asA<MultiCLJComponent>();
                nrg.addComponent("intra_coulomb", nrgs[clj.coulomb()] * kcal_per_mol);
                nrg.addComponent("intra_LJ", nrgs[clj.lj()] * kcal_per_mol);
            }
            else
            {
                const auto &comps = ffield.read().components();

                if (comps.isA<InternalComponent>())
                {
                    const auto &internal = comps.asA<InternalComponent>();

                    nrg.addComponent("bond", nrgs[internal.bond()] * kcal_per_mol);
                    nrg.addComponent("angle", nrgs[internal.angle()] * kcal_per_mol);
                    nrg.addComponent("dihedral", nrgs[internal.dihedral()] * kcal_per_mol);
                    nrg.addComponent("improper", nrgs[internal.improper()] * kcal_per_mol);
                    nrg.addComponent("urey-bradley", nrgs[internal.ureyBradley()] * kcal_per_mol);
                    nrg.addComponent("1-4_coulomb", nrgs[internal.intra14Coulomb()] * kcal_per_mol);
                    nrg.addComponent("1-4_LJ", nrgs[internal.intra14LJ()] * kcal_per_mol);
                }
                else if (comps.isA<MultiCLJComponent>())
                {
                    const auto &clj = comps.asA<MultiCLJComponent>();

                    nrg.addComponent("coulomb", nrgs[clj.coulomb()] * kcal_per_mol);
                    nrg.addComponent("LJ", nrgs[clj.lj()] * kcal_per_mol);
                }
            }
        }

        double delta = total - nrg.value();

        if (std::abs(delta) > 1e-8)
        {
            nrg.addComponent("other", (total - nrg.value()) * kcal_per_mol);
        }

        return nrg;
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol)
    {
        return calculate_energy(mol, PropertyMap());
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol, const SireBase::PropertyMap &map)
    {
        auto ff = create_forcefield(mol, map);
        return calculate_energy(ff);
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols)
    {
        return calculate_energy(mols, PropertyMap());
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols, const SireBase::PropertyMap &map)
    {
        auto ff = create_forcefield(mols, map);
        return calculate_energy(ff);
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0, const SireMol::MoleculeView &mol1)
    {
        return calculate_energy(mol0, mol1, PropertyMap());
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0, const SireMol::MoleculeView &mol1,
                                                   const SireBase::PropertyMap &map)
    {
        auto ff = create_forcefield(mol0, mol1, map);
        return calculate_energy(ff);
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0, const SireMol::Molecules &mols1)
    {
        return calculate_energy(mol0, mols1, PropertyMap());
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::MoleculeView &mol0, const SireMol::Molecules &mols1,
                                                   const SireBase::PropertyMap &map)
    {
        auto ff = create_forcefield(mol0, mols1, map);
        return calculate_energy(ff);
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols0, const SireMol::Molecules &mols1)
    {
        return calculate_energy(mols0, mols1, PropertyMap());
    }

    SIRESYSTEM_EXPORT GeneralUnit calculate_energy(const SireMol::Molecules &mols0, const SireMol::Molecules &mols1,
                                                   const SireBase::PropertyMap &map)
    {
        auto ff = create_forcefield(mols0, mols1, map);
        return calculate_energy(ff);
    }

    SIRESYSTEM_EXPORT QVector<GeneralUnit> calculate_trajectory_energy(const ForceFields &ff, const QList<qint64> &frames,
                                                                       const PropertyMap &map)
    {
        QVector<GeneralUnit> nrgs;

        if (frames.isEmpty())
            return nrgs;

        nrgs.resize(frames.count());

        auto nrgs_data = nrgs.data();

        QVector<qint64> local_frames = frames.toVector();
        auto frames_data = local_frames.constData();

        tbb::parallel_for(tbb::blocked_range<int>(0, frames.count()), [&](const tbb::blocked_range<int> &r)
                          {
        ForceFields local_ff = ff;

        for (int i = r.begin(); i < r.end(); ++i)
        {
            local_ff.loadFrame(frames_data[i], map);
            nrgs_data[i] = calculate_energy(local_ff);
        } });

        return nrgs;
    }

    SIRESYSTEM_EXPORT QVector<QVector<GeneralUnit>> calculate_trajectory_energies(const QVector<ForceFields> &ffs,
                                                                                  const QList<qint64> &frames,
                                                                                  const PropertyMap &map)
    {
        QVector<QVector<GeneralUnit>> nrgs;

        if (frames.isEmpty() or ffs.isEmpty())
            return nrgs;

        nrgs.resize(ffs.count());

        auto nrgs_data = nrgs.data();

        auto local_ffs = ffs.constData();

        QVector<qint64> local_frames = frames.toVector();
        auto frame_data = local_frames.constData();
        const int nframes = local_frames.count();

        tbb::parallel_for(tbb::blocked_range<int>(0, ffs.count()), [&](const tbb::blocked_range<int> &r)
                          {
        for (int i = r.begin(); i < r.end(); ++i)
        {
            QVector<GeneralUnit> ff_nrgs(nframes);
            auto ff_nrgs_data = ff_nrgs.data();

            tbb::parallel_for(tbb::blocked_range<int>(0, nframes), [&](const tbb::blocked_range<int> &r2) {
                ForceFields local_ff = local_ffs[i];

                for (int j = r2.begin(); j < r2.end(); ++j)
                {
                    local_ff.loadFrame(frame_data[j], map);
                    ff_nrgs_data[j] = calculate_energy(local_ff);
                }
            });

            nrgs_data[i] = ff_nrgs;
        } });

        return nrgs;
    }

} // end of namespace SireMM
