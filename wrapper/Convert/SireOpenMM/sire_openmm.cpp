
#include "sire_openmm.h"

#include <OpenMM.h>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireSystem/forcefieldinfo.h"

#include "SireMol/core.h"
#include "SireMol/moleditor.h"
#include "SireMol/atomelements.h"
#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomproperty.hpp"
#include "SireMol/connectivity.h"
#include "SireMol/bondid.h"
#include "SireMol/bondorder.h"
#include "SireMol/atomvelocities.h"

#include "SireMM/atomljs.h"
#include "SireMM/selectorbond.h"
#include "SireMM/amberparams.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireCAS/lambdaschedule.h"

#include "SireMaths/vector.h"

#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "tostring.h"

#include "openmmmolecule.h"

#include <QDebug>

using SireBase::PropertyMap;
using SireCAS::LambdaSchedule;
using SireMol::Molecule;
using SireMol::MolNum;
using SireMol::SelectorMol;
using SireSystem::ForceFieldInfo;

namespace SireOpenMM
{
    bool use_parallel(int n, const SireBase::PropertyMap &map)
    {
        if (n <= 16)
            return false;

        if (map["parallel"].hasValue())
        {
            return map["parallel"].value().asABoolean();
        }

        return true;
    }

    ////
    //// Implementation of OpenMMMetaData
    ////

    OpenMMMetaData::OpenMMMetaData()
    {
    }

    OpenMMMetaData::OpenMMMetaData(const SireMol::SelectorM<SireMol::Atom> &i,
                                   std::shared_ptr<std::vector<OpenMM::Vec3>> c,
                                   std::shared_ptr<std::vector<OpenMM::Vec3>> v,
                                   std::shared_ptr<std::vector<OpenMM::Vec3>> b,
                                   const LambdaLever &l)
        : atom_index(i), coords(c), vels(v), boxvecs(b), lambda_lever(l)
    {
    }

    OpenMMMetaData::~OpenMMMetaData()
    {
    }

    SireMol::SelectorM<SireMol::Atom> OpenMMMetaData::index() const
    {
        return atom_index;
    }

    LambdaLever OpenMMMetaData::lambdaLever() const
    {
        return lambda_lever;
    }

    bool OpenMMMetaData::hasCoordinates() const
    {
        return coords.get() != 0;
    }

    bool OpenMMMetaData::hasVelocities() const
    {
        return vels.get() != 0;
    }

    bool OpenMMMetaData::hasBoxVectors() const
    {
        return boxvecs.get() != 0;
    }

    const std::vector<OpenMM::Vec3> &OpenMMMetaData::coordinates() const
    {
        if (coords.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no coordinates available!"),
                                                CODELOC);

        return *coords;
    }

    const std::vector<OpenMM::Vec3> &OpenMMMetaData::velocities() const
    {
        if (vels.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no velocities available!"),
                                                CODELOC);

        return *vels;
    }

    const std::vector<OpenMM::Vec3> &OpenMMMetaData::boxVectors() const
    {
        if (boxvecs.get() == 0)
            throw SireError::incompatible_error(QObject::tr(
                                                    "There are no box vectors available!"),
                                                CODELOC);

        return *boxvecs;
    }

    ////
    //// Implementation of standalone functions
    ////

    SelectorMol openmm_system_to_sire(const OpenMM::System &mols,
                                      const PropertyMap &map)
    {
        throw SireError::incomplete_code(QObject::tr(
                                             "Still need to write openmm_to_sire"),
                                         CODELOC);

        return SelectorMol();
    }

    OpenMMMetaData sire_to_openmm_system(OpenMM::System &system,
                                         const SelectorMol &mols,
                                         const PropertyMap &map)
    {
        // we can assume that an empty system has been passed to us

        // we will get parameters from the map (e.g. cutoffs etc)
        ForceFieldInfo ffinfo(mols, map);

        // extract the data from all of the molecules
        const int nmols = mols.count();

        if (nmols == 0)
        {
            // nothing to do
            return OpenMMMetaData();
        }

        // whether or not to ignore perturbations
        bool ignore_perturbations = false;
        bool any_perturbable = false;

        if (map.specified("ignore_perturbations"))
        {
            ignore_perturbations = map["ignore_perturbations"].value().asABoolean();
        }

        // Extract all of the data needed by OpenMM from the Sire
        // molecules into some temporary OpenMMMolecule objects
        QVector<OpenMMMolecule> openmm_mols(nmols);
        auto openmm_mols_data = openmm_mols.data();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, mols.count()), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    openmm_mols_data[i] = OpenMMMolecule(mols[i], map);
                } });
        }
        else
        {
            for (int i = 0; i < nmols; ++i)
            {
                openmm_mols_data[i] = OpenMMMolecule(mols[i], map);
            }
        }

        // check to see if there are any perturbable molecules
        if (not ignore_perturbations)
        {
            for (int i = 0; i < nmols; ++i)
            {
                if (openmm_mols_data[i].isPerturbable())
                {
                    any_perturbable = true;
                    break;
                }
            }
        }

        // create all the OpenMM forcefields - ownership is taken by 'system'
        // (memory leak until we give this to 'system')

        // CLJ energy between all non-perturbable atoms
        OpenMM::NonbondedForce *cljff = new OpenMM::NonbondedForce();

        bool use_dispersion_correction = false;

        if (map.specified("use_dispersion_correction"))
        {
            use_dispersion_correction = map["use_dispersion_correction"].value().asABoolean();
        }

        // note that this will be very slow for perturbable systems, as
        // it needs recalculating for every change of lambda
        cljff->setUseDispersionCorrection(use_dispersion_correction);

        // create the periodic box vectors
        std::shared_ptr<std::vector<OpenMM::Vec3>> boxvecs;

        if (ffinfo.space().isPeriodic())
        {
            boxvecs.reset(new std::vector<OpenMM::Vec3>(3));
            auto boxvecs_data = boxvecs->data();

            if (ffinfo.space().isA<SireVol::PeriodicBox>())
            {
                const auto &space = ffinfo.space().asA<SireVol::PeriodicBox>();

                const double x = space.dimensions()[0] * OpenMM::NmPerAngstrom;
                const double y = space.dimensions()[1] * OpenMM::NmPerAngstrom;
                const double z = space.dimensions()[2] * OpenMM::NmPerAngstrom;

                boxvecs_data[0] = OpenMM::Vec3(x, 0, 0);
                boxvecs_data[1] = OpenMM::Vec3(0, y, 0);
                boxvecs_data[2] = OpenMM::Vec3(0, 0, z);
            }
            else if (ffinfo.space().isA<SireVol::TriclinicBox>())
            {
                const auto &space = ffinfo.space().asA<SireVol::TriclinicBox>();

                // Get the three triclinic box vectors.
                const auto v0 = space.vector0();
                const auto v1 = space.vector1();
                const auto v2 = space.vector2();

                // Get cell matrix components in nm.
                const double xx = v0.x() * OpenMM::NmPerAngstrom;
                const double xy = v0.y() * OpenMM::NmPerAngstrom;
                const double xz = v0.z() * OpenMM::NmPerAngstrom;
                const double yx = v1.x() * OpenMM::NmPerAngstrom;
                const double yy = v1.y() * OpenMM::NmPerAngstrom;
                const double yz = v1.z() * OpenMM::NmPerAngstrom;
                const double zx = v2.x() * OpenMM::NmPerAngstrom;
                const double zy = v2.y() * OpenMM::NmPerAngstrom;
                const double zz = v2.z() * OpenMM::NmPerAngstrom;

                boxvecs_data[0] = OpenMM::Vec3(xx, xy, xz);
                boxvecs_data[1] = OpenMM::Vec3(yx, yy, yz);
                boxvecs_data[2] = OpenMM::Vec3(zx, zy, zz);
            }

            system.setDefaultPeriodicBoxVectors(boxvecs_data[0],
                                                boxvecs_data[1],
                                                boxvecs_data[2]);
        }

        if (ffinfo.hasCutoff())
        {
            const auto typ = ffinfo.cutoffType();

            // need to set cutoff for perturbable forcefields

            if (typ == "PME" or typ == "EWALD")
            {
                if (not ffinfo.space().isPeriodic())
                {
                    throw SireError::incompatible_error(QObject::tr(
                                                            "You cannot use Ewald or PME with the non-periodic space %1.")
                                                            .arg(ffinfo.space().toString()),
                                                        CODELOC);
                }

                auto nbmethod = OpenMM::NonbondedForce::PME;

                if (typ != "PME")
                    nbmethod = OpenMM::NonbondedForce::Ewald;

                cljff->setNonbondedMethod(nbmethod);

                double tolerance = ffinfo.getParameter("tolerance").value();

                if (tolerance <= 0)
                    tolerance = 0.001;

                cljff->setEwaldErrorTolerance(tolerance);
            }
            else if (typ == "REACTION_FIELD")
            {
                auto nbmethod = OpenMM::NonbondedForce::CutoffPeriodic;

                if (not ffinfo.space().isPeriodic())
                {
                    nbmethod = OpenMM::NonbondedForce::CutoffNonPeriodic;
                }

                auto dielectric = ffinfo.getParameter("dielectric").value();

                if (dielectric <= 0)
                    dielectric = 78.3;

                cljff->setNonbondedMethod(nbmethod);
                cljff->setReactionFieldDielectric(dielectric);
            }
            else if (typ == "CUTOFF")
            {
                // use reaction field for non-periodic spaces, and PME for periodic
                if (ffinfo.space().isPeriodic())
                {
                    const auto nbmethod = OpenMM::NonbondedForce::PME;
                    cljff->setNonbondedMethod(nbmethod);

                    double tolerance = ffinfo.getParameter("tolerance").value();

                    if (tolerance <= 0)
                        tolerance = 0.001;

                    cljff->setEwaldErrorTolerance(tolerance);
                }
                else
                {
                    const auto nbmethod = OpenMM::NonbondedForce::CutoffNonPeriodic;
                    cljff->setNonbondedMethod(nbmethod);

                    double dielectric = ffinfo.getParameter("dielectric").value();

                    if (dielectric <= 0)
                        dielectric = 78.3;

                    cljff->setReactionFieldDielectric(dielectric);
                }
            }

            const auto cutoff = ffinfo.cutoff().to(SireUnits::nanometers);
            cljff->setCutoffDistance(cutoff);
        }

        // also populate a LambaLever for any perturbable molecules
        LambdaLever lambda_lever;

        if (any_perturbable)
        {
            if (map.specified("schedule"))
            {
                lambda_lever.setSchedule(
                    map["schedule"].value().asA<LambdaSchedule>());
            }
            else
            {
                lambda_lever.setSchedule(
                    LambdaSchedule::standard_morph());
            }
        }

        // make sure that we tell the lever the index of the named
        // forcefields when they are added
        lambda_lever.setForceIndex("clj", system.addForce(cljff));

        OpenMM::HarmonicBondForce *bondff = new OpenMM::HarmonicBondForce();
        lambda_lever.setForceIndex("bond", system.addForce(bondff));

        OpenMM::HarmonicAngleForce *angff = new OpenMM::HarmonicAngleForce();
        lambda_lever.setForceIndex("angle", system.addForce(angff));

        OpenMM::PeriodicTorsionForce *dihff = new OpenMM::PeriodicTorsionForce();
        lambda_lever.setForceIndex("torsion", system.addForce(dihff));

        // now let's define the forces used for managing the creation
        // or deletion of ghost atoms
        OpenMM::CustomBondForce *ghost_14ff = 0;
        OpenMM::CustomNonbondedForce *ghost_ghostff = 0;
        OpenMM::CustomNonbondedForce *ghost_nonghostff = 0;

        if (any_perturbable)
        {
            double shift_delta = 1.0;

            if (map.specified("shift_delta"))
            {
                shift_delta = map["shift_delta"].value().asADouble();
            }

            if (shift_delta < 0)
                shift_delta = 0;

            int coulomb_power = 0;

            if (map.specified("coulomb_power"))
            {
                coulomb_power = map["coulomb_power"].value().asAnInteger();
            }

            if (coulomb_power < 0)
                coulomb_power = 0;
            else if (coulomb_power > 4)
                coulomb_power = 4;

            auto coulomb_power_expression = [](const QString &alpha, int power)
            {
                if (power == 0)
                    return QString("1");
                else if (power == 1)
                    return QString("(1-%1)").arg(alpha);
                else if (power == 2)
                    return QString("(1-%1)*(1-%1)").arg(alpha);
                else
                    return QString("(1-%1)^%2").arg(alpha).arg(power);
            };

            // see below for the description of this energy expression
            const auto nb14_expression = QString(
                                             "coul_nrg+lj_nrg;"
                                             "coul_nrg=138.9354558466661*q*(((%1)/sqrt(alpha+r^2))-(1.0/r));"
                                             "lj_nrg=four_epsilon*((sig6^2)-sig6);"
                                             "sig6=(sigma^6)/(((sigma*delta) + r^2)^3);"
                                             "delta=%2*alpha;")
                                             .arg(coulomb_power_expression("alpha", coulomb_power))
                                             .arg(shift_delta)
                                             .toStdString();

            ghost_14ff = new OpenMM::CustomBondForce(nb14_expression);

            ghost_14ff->addPerBondParameter("q");
            ghost_14ff->addPerBondParameter("sigma");
            ghost_14ff->addPerBondParameter("four_epsilon");
            ghost_14ff->addPerBondParameter("alpha");

            // short-range intramolecular term that should not use
            // periodic boundaries or cutoffs
            ghost_14ff->setUsesPeriodicBoundaryConditions(false);

            // this uses the following potentials
            //            Zacharias and McCammon, J. Chem. Phys., 1994, and also,
            //            Michel et al., JCTC, 2007
            //
            //   V_{LJ}(r) = 4 epsilon [ ( sigma^12 / (delta*sigma + r^2)^6 ) -
            //                           ( sigma^6  / (delta*sigma + r^2)^3 ) ]
            //
            //   delta = shift_delta * alpha
            //
            //   V_{coul}(r) = (1-alpha)^n q_i q_j / 4 pi eps_0 (alpha+r^2)^(1/2)
            //
            // Note that we pre-calculate delta as a forcefield parameter,
            // and also supply half_sigma and two_sqrt_epsilon to save some
            // cycles
            //
            // Note also that we subtract the normal coulomb energy as this
            // is calculated during the standard NonbondedForce
            //
            // 138.9354558466661 is the constant needed to get energies in
            // kJ mol-1 given the units of charge (|e|) and distance (nm)
            //
            const auto clj_expression = QString("coul_nrg+lj_nrg;"
                                                "coul_nrg=138.9354558466661*q1*q2*(((%1)/sqrt(max_alpha+r^2))-(1.0/r));"
                                                "lj_nrg=two_sqrt_epsilon1*two_sqrt_epsilon2*((sig6^2)-sig6);"
                                                "sig6=(sigma^6)/(((sigma*delta) + r^2)^3);"
                                                "delta=%2*max_alpha;"
                                                "max_alpha=max(alpha1, alpha2);"
                                                "sigma=half_sigma1+half_sigma2;")
                                            .arg(coulomb_power_expression("max_alpha", coulomb_power))
                                            .arg(shift_delta)
                                            .toStdString();

            ghost_ghostff = new OpenMM::CustomNonbondedForce(clj_expression);
            ghost_nonghostff = new OpenMM::CustomNonbondedForce(clj_expression);

            ghost_ghostff->addPerParticleParameter("q");
            ghost_ghostff->addPerParticleParameter("half_sigma");
            ghost_ghostff->addPerParticleParameter("two_sqrt_epsilon");
            ghost_ghostff->addPerParticleParameter("alpha");

            ghost_nonghostff->addPerParticleParameter("q");
            ghost_nonghostff->addPerParticleParameter("half_sigma");
            ghost_nonghostff->addPerParticleParameter("two_sqrt_epsilon");
            ghost_nonghostff->addPerParticleParameter("alpha");

            // this will be slow if switched on, as it needs recalculating
            // for every change in parameters
            ghost_ghostff->setUseLongRangeCorrection(use_dispersion_correction);
            ghost_nonghostff->setUseLongRangeCorrection(use_dispersion_correction);

            if (ffinfo.hasCutoff())
            {
                if (ffinfo.space().isPeriodic())
                {
                    ghost_ghostff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
                    ghost_nonghostff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
                }
                else
                {
                    ghost_ghostff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
                    ghost_nonghostff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
                }

                ghost_ghostff->setCutoffDistance(ffinfo.cutoff().to(SireUnits::nanometers));
                ghost_nonghostff->setCutoffDistance(ffinfo.cutoff().to(SireUnits::nanometers));
            }
            else
            {
                ghost_ghostff->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
                ghost_nonghostff->setNonbondedMethod(OpenMM::CustomNonbondedForce::NoCutoff);
            }

            lambda_lever.setForceIndex("ghost/ghost", system.addForce(ghost_ghostff));
            lambda_lever.setForceIndex("ghost/non-ghost", system.addForce(ghost_nonghostff));
            lambda_lever.setForceIndex("ghost-14", system.addForce(ghost_14ff));
        }

        // now the forcefields relating to restraints - first positional
        // restraints - harmonic potentials keeping atoms close to
        // fixed anchor points
        OpenMM::HarmonicBondForce *positional_restraintff = 0;

        if (map.specified("positional_restraints"))
        {
            positional_restraintff = new OpenMM::HarmonicBondForce();
            positional_restraintff->setUsesPeriodicBoundaryConditions(true);

            lambda_lever.setForceIndex("positional_restraint",
                                       system.addForce(positional_restraintff));
        }

        // now classic distance, angle or torsion restraints between
        // any pairs, triples or quads of atoms (including intermolecular)
        OpenMM::HarmonicBondForce *distance_restraintff = 0;
        OpenMM::HarmonicAngleForce *angle_restraintff = 0;
        OpenMM::PeriodicTorsionForce *torsion_restraintff = 0;

        if (map.specified("distance_restraints"))
        {
            distance_restraintff = new OpenMM::HarmonicBondForce();
            distance_restraintff->setUsesPeriodicBoundaryConditions(true);

            lambda_lever.setForceIndex("distance_restraint",
                                       system.addForce(distance_restraintff));
        }

        if (map.specified("angle_restraints"))
        {
            angle_restraintff = new OpenMM::HarmonicAngleForce();
            angle_restraintff->setUsesPeriodicBoundaryConditions(true);

            lambda_lever.setForceIndex("angle_restraint",
                                       system.addForce(angle_restraintff));
        }

        if (map.specified("torsion_restraints"))
        {
            torsion_restraintff = new OpenMM::PeriodicTorsionForce();
            torsion_restraintff->setUsesPeriodicBoundaryConditions(true);

            lambda_lever.setForceIndex("torsion_restraint",
                                       system.addForce(torsion_restraintff));
        }

        // now the Boresch restraints that can be used to hold a ligand
        // in a defined position and orientation relative to a receptor
        OpenMM::HarmonicBondForce *boresch_distance_restraintff = 0;
        OpenMM::HarmonicAngleForce *boresch_angle_restraintff = 0;
        OpenMM::CustomTorsionForce *boresch_torsion_restraintff = 0;

        if (map.specified("boresch_restraints"))
        {
            boresch_distance_restraintff = new OpenMM::HarmonicBondForce();
            boresch_angle_restraintff = new OpenMM::HarmonicAngleForce();

            const auto torsion_expression = QString(
                                                "force_const*min(dtheta, 2*pi-dtheta)^2;"
                                                "dtheta = abs(theta-equil_val);"
                                                "pi = 3.1415926535;")
                                                .toStdString();

            boresch_torsion_restraintff = new OpenMM::CustomTorsionForce(torsion_expression);

            boresch_distance_restraintff->setUsesPeriodicBoundaryConditions(true);
            boresch_angle_restraintff->setUsesPeriodicBoundaryConditions(true);
            boresch_torsion_restraintff->setUsesPeriodicBoundaryConditions(true);

            lambda_lever.setForceIndex("boresch_distance",
                                       system.addForce(boresch_distance_restraintff));
            lambda_lever.setForceIndex("boresch_angle",
                                       system.addForce(boresch_angle_restraintff));
            lambda_lever.setForceIndex("boresch_torsion",
                                       system.addForce(boresch_torsion_restraintff));
        }

        // Now copy data from the temporary OpenMMMolecule objects
        // into these forcefields
        // (will deal with restraints, light atoms and virtual sites later)

        int start_index = 0;

        std::vector<std::pair<int, int>> bond_pairs;
        std::vector<std::tuple<int, int, double, double, double>> exception_params;

        // get the 1-4 scaling factors from the first molecule
        const double coul_14_scl = openmm_mols_data[0].ffinfo.electrostatic14ScaleFactor();
        const double lj_14_scl = openmm_mols_data[0].ffinfo.vdw14ScaleFactor();

        QVector<int> start_indexes(nmols);

        // save the atoms in the order they are added to the system
        SireMol::SelectorM<SireMol::Atom> atom_index;

        // the index to the perturbable molecule for the specified molecule
        QHash<int, int> idx_to_pert_idx;

        // just a holder for all of the custom parameters
        // (prevents us having to continually re-allocate it)
        std::vector<double> custom_params = {0.0, 0.0, 0.0, 0.0};

        // the sets of ghost atoms and non-ghost atoms
        std::set<int> ghost_atoms;
        std::set<int> non_ghost_atoms;

        // the set of all ghost atoms, with the value
        // indicating if this is a from-ghost (true) or
        // a to-ghost (false)
        QHash<int, bool> ghost_is_from;

        for (int i = 0; i < nmols; ++i)
        {
            start_indexes[i] = start_index;
            const auto &mol = openmm_mols_data[i];

            // add all of the atoms to the index of atoms
            // (we guarantee to add them in atomidx order)
            atom_index += mol.atoms;

            if (std::abs(mol.ffinfo.electrostatic14ScaleFactor() - coul_14_scl) > 0.001 or
                std::abs(mol.ffinfo.vdw14ScaleFactor() - lj_14_scl) > 0.001)
            {
                throw SireError::incompatible_error(QObject::tr(
                                                        "We cannot create the OpenMM system because the forcefields "
                                                        "of %1 and %2 (%3 and %4) are not compatible. They have "
                                                        "different 1-4 non-bonded scale factors.")
                                                        .arg(mols[0].toString())
                                                        .arg(mols[i].toString())
                                                        .arg(openmm_mols_data[0].ffinfo.toString())
                                                        .arg(mol.ffinfo.toString()),
                                                    CODELOC);
            }

            QHash<QString, qint32> start_indicies;

            // is this a perturbable molecule (and we haven't disabled perturbations)?
            if (any_perturbable and mol.isPerturbable())
            {
                // add a perturbable molecule, recording the start index
                // for each of the forcefields
                start_indicies.reserve(7);

                start_indicies.insert("clj", start_index);
                start_indicies.insert("ghost/ghost", start_index);
                start_indicies.insert("ghost/non-ghost", start_index);
                start_indicies.insert("bond", bondff->getNumBonds());
                start_indicies.insert("angle", angff->getNumAngles());
                start_indicies.insert("torsion", dihff->getNumTorsions());

                auto pert_idx = lambda_lever.addPerturbableMolecule(mol,
                                                                    start_indicies);

                idx_to_pert_idx.insert(i, pert_idx);
            }

            // first the atom parameters
            auto masses_data = mol.masses.constData();
            auto cljs_data = mol.cljs.constData();
            auto alphas_data = mol.alphas.constData();

            if (any_perturbable and mol.isPerturbable())
            {
                for (int j = 0; j < mol.molinfo.nAtoms(); ++j)
                {
                    const bool is_from_ghost = mol.from_ghost_idxs.contains(j);
                    const bool is_to_ghost = mol.to_ghost_idxs.contains(j);

                    system.addParticle(masses_data[j]);
                    const int atom_index = start_index + j;

                    const auto &clj = cljs_data[j];

                    // reduced_q
                    custom_params[0] = std::get<0>(clj);
                    // half_sigma
                    custom_params[1] = 0.5 * std::get<1>(clj);
                    // two_sqrt_epsilon
                    custom_params[2] = 2.0 * std::sqrt(std::get<2>(clj));
                    // alpha
                    custom_params[3] = alphas_data[j];

                    ghost_ghostff->addParticle(custom_params);
                    ghost_nonghostff->addParticle(custom_params);

                    if (is_from_ghost or is_to_ghost)
                    {
                        ghost_atoms.insert(atom_index);
                        ghost_is_from.insert(atom_index, is_from_ghost);

                        // don't include the LJ energy as this will be
                        // calculated using the ghost forcefields
                        // (the ghost forcefields include a coulomb term
                        //  that subtracts from whatever was calculated here)
                        cljff->addParticle(std::get<0>(clj), 0.0, 0.0);
                    }
                    else
                    {
                        cljff->addParticle(std::get<0>(clj), std::get<1>(clj),
                                           std::get<2>(clj));
                        non_ghost_atoms.insert(atom_index);
                    }
                }
            }
            else
            {
                for (int j = 0; j < mol.molinfo.nAtoms(); ++j)
                {
                    system.addParticle(masses_data[j]);
                    const int atom_index = start_index + j;

                    const auto &clj = cljs_data[j];

                    cljff->addParticle(std::get<0>(clj), std::get<1>(clj), std::get<2>(clj));

                    if (any_perturbable)
                    {
                        // reduced charge
                        custom_params[0] = std::get<0>(clj);
                        // half_sigma
                        custom_params[1] = 0.5 * std::get<1>(clj);
                        // two_sqrt_epsilon
                        custom_params[2] = 2.0 * std::sqrt(std::get<2>(clj));
                        // alpha - is zero for non-ghost atoms
                        custom_params[3] = 0.0;
                        ghost_ghostff->addParticle(custom_params);
                        ghost_nonghostff->addParticle(custom_params);
                        non_ghost_atoms.insert(atom_index);
                    }
                }
            }

            if (not(any_perturbable and mol.isPerturbable()))
            {
                // only save the bond pairs for non-perturbable molecules.
                // This is because we want to know which exceptions
                // are created for perturbable molecules, as we will
                // need to manage them ourselves
                for (const auto &bond : mol.bond_pairs)
                {
                    bond_pairs.push_back(std::make_pair(std::get<0>(bond) + start_index,
                                                        std::get<1>(bond) + start_index));
                }
            }

            // now bond parameters
            for (const auto &bond : mol.bond_params)
            {
                bondff->addBond(std::get<0>(bond) + start_index,
                                std::get<1>(bond) + start_index,
                                std::get<2>(bond), std::get<3>(bond));
            }

            // now the angles
            for (const auto &ang : mol.ang_params)
            {
                angff->addAngle(std::get<0>(ang) + start_index,
                                std::get<1>(ang) + start_index,
                                std::get<2>(ang) + start_index,
                                std::get<3>(ang), std::get<4>(ang));
            }

            // now the dihedrals and impropers
            for (const auto &dih : mol.dih_params)
            {
                dihff->addTorsion(std::get<0>(dih) + start_index,
                                  std::get<1>(dih) + start_index,
                                  std::get<2>(dih) + start_index,
                                  std::get<3>(dih) + start_index,
                                  std::get<4>(dih), std::get<5>(dih), std::get<6>(dih));
            }

            // now constraints
            if (any_perturbable and mol.isPerturbable())
            {
                // we may want to select out constraints that involve
                // perturbing atoms... For now we will include them all,
                // but will use the current coordinates to set the
                // constrained values, rather than using the
                // forcefield parameters. This will work as long
                // as the user has minimised the system at the desired
                // lambda value before running dynamics...
                for (const auto &constraint : mol.constraints)
                {
                    const auto coords0 = mol.coords[std::get<0>(constraint)];
                    const auto coords1 = mol.coords[std::get<1>(constraint)];

                    const auto delta = coords1 - coords0;

                    const auto length = std::sqrt((delta[0] * delta[0]) +
                                                  (delta[1] * delta[1]) +
                                                  (delta[2] * delta[2]));

                    system.addConstraint(
                        std::get<0>(constraint) + start_index,
                        std::get<1>(constraint) + start_index,
                        length);
                }
            }
            else
            {
                for (const auto &constraint : mol.constraints)
                {
                    system.addConstraint(std::get<0>(constraint) + start_index,
                                         std::get<1>(constraint) + start_index,
                                         std::get<2>(constraint));
                }
            }

            start_index += mol.masses.count();
        }

        // add the exceptions automatically for non-perturbable molecules
        cljff->createExceptionsFromBonds(bond_pairs, coul_14_scl, lj_14_scl);

        // add additional exceptions, including all exceptions for
        // perturbable molecules (perturbable molecules are handled
        // completely because the exceptions may change during
        // the perturbation)
        for (int i = 0; i < nmols; ++i)
        {
            int start_index = start_indexes[i];
            const auto &mol = openmm_mols_data[i];

            QVector<std::pair<int, int>> exception_idxs;

            const bool is_perturbable = any_perturbable and mol.isPerturbable();

            if (is_perturbable)
            {
                exception_idxs = QVector<std::pair<int, int>>(mol.exception_params.count(),
                                                              std::make_pair(-1, -1));
            }

            for (int j = 0; j < mol.exception_params.count(); ++j)
            {
                const auto &param = mol.exception_params[j];

                const auto atom0 = std::get<0>(param);
                const auto atom1 = std::get<1>(param);
                const auto coul_14_scale = std::get<2>(param);
                const auto lj_14_scale = std::get<3>(param);

                auto p = mol.getException(atom0, atom1,
                                          start_index,
                                          coul_14_scale,
                                          lj_14_scale);

                if (is_perturbable)
                {
                    const bool atom0_is_ghost = mol.isGhostAtom(atom0);
                    const bool atom1_is_ghost = mol.isGhostAtom(atom1);

                    int idx = -1;
                    int nbidx = -1;

                    if (atom0_is_ghost or atom1_is_ghost)
                    {
                        // don't include the LJ term, as this is calculated
                        // elsewhere - note that we need to use 1e-9 to
                        // make sure that OpenMM doesn't eagerly remove
                        // this, and cause "changed excluded atoms" warnings
                        idx = cljff->addException(std::get<0>(p), std::get<1>(p),
                                                  std::get<2>(p), 1e-9,
                                                  1e-9, true);

                        if (coul_14_scl != 0 or lj_14_scl != 0)
                        {
                            // this is a 1-4 interaction that should be added
                            // to the ghost-14 forcefield
                            if (ghost_14ff != 0)
                            {
                                // parameters are q, sigma, four_epsilon and alpha
                                std::vector<double> params14 =
                                    {std::get<2>(p), std::get<3>(p),
                                     4.0 * std::get<4>(p), 0.0};

                                if (params14[0] == 0)
                                    // cannot use zero params in case they are
                                    // eagerly removed
                                    params14[0] = 1e-9;

                                if (params14[1] == 0)
                                    params14[1] = 1e-9;

                                nbidx = ghost_14ff->addBond(std::get<0>(p),
                                                            std::get<1>(p),
                                                            params14);
                            }
                        }
                    }
                    else
                    {
                        idx = cljff->addException(std::get<0>(p), std::get<1>(p),
                                                  std::get<2>(p), std::get<3>(p),
                                                  std::get<4>(p), true);
                    }

                    exception_idxs[j] = std::make_pair(idx, nbidx);

                    // remove this interaction from the ghost forcefields
                    if (ghost_ghostff != 0)
                    {
                        ghost_ghostff->addExclusion(std::get<0>(p), std::get<1>(p));
                        ghost_nonghostff->addExclusion(std::get<0>(p), std::get<1>(p));
                    }
                }
                else
                {
                    cljff->addException(std::get<0>(p), std::get<1>(p),
                                        std::get<2>(p), std::get<3>(p),
                                        std::get<4>(p), true);
                }
            }

            if (is_perturbable)
            {
                auto pert_idx = idx_to_pert_idx.value(i, openmm_mols.count() + 1);
                lambda_lever.setExceptionIndicies(pert_idx,
                                                  "clj", exception_idxs);
            }
        }

        if (ghost_ghostff != 0)
        {
            // set up the interaction groups - ghost / non-ghost
            //                                 ghost / ghost
            ghost_ghostff->addInteractionGroup(ghost_atoms, ghost_atoms);
            ghost_nonghostff->addInteractionGroup(ghost_atoms, non_ghost_atoms);
        }

        // set up all of the restraints
        if (positional_restraintff != 0)
        {
            // we have some positional restraints
            const auto restraints = map["positional_restraints"].value().asA<PositionalRestraints>();

            // go through each restraint and create an anchor atom, then add
            // the restraint parameters
        }

        if (distance_restraintff != 0)
        {
            // we have some distance restraints
        }

        if (angle_restraintff != 0)
        {
            // we have some angle restraints
        }

        if (torsion_restraintff != 0)
        {
            // we have some torsion restraints
        }

        if (boresch_distance_restraintff != 0)
        {
            // we have some boresch restraints - these are a set
            // that will always use all three boresch restraint forces
        }

        // see if we want to remove COM motion
        const auto com_remove_prop = map["com_reset_frequency"];

        if (com_remove_prop.hasValue())
        {
            const int freq = com_remove_prop.value().asAnInteger();

            if (freq > 0)
            {
                OpenMM::CMMotionRemover *com_remover = new OpenMM::CMMotionRemover(freq);
                system.addForce(com_remover);
            }
        }

        // now get the coordinates and velocities
        std::shared_ptr<std::vector<OpenMM::Vec3>> coords, vels;

        const int natoms = start_index;

        coords.reset(new std::vector<OpenMM::Vec3>(natoms));
        vels.reset(new std::vector<OpenMM::Vec3>(natoms));

        auto coords_data = coords->data();
        auto vels_data = vels->data();

        const int *start_indexes_data = start_indexes.constData();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const int start_index = start_indexes_data[i];
                    const auto &mol = openmm_mols_data[i];
                    mol.copyInCoordsAndVelocities(coords_data + start_index,
                                                  vels_data + start_index);
                } });
        }
        else
        {
            for (int i = 0; i < nmols; ++i)
            {
                const int start_index = start_indexes_data[i];
                const auto &mol = openmm_mols_data[i];
                mol.copyInCoordsAndVelocities(coords_data + start_index,
                                              vels_data + start_index);
            }
        }

        // now copy in the positional restraint coordinates (if any)
        // ...

        return OpenMMMetaData(atom_index, coords, vels, boxvecs, lambda_lever);
    }

    void set_openmm_coordinates_and_velocities(OpenMM::Context &context,
                                               const OpenMMMetaData &metadata)
    {
        if (metadata.hasCoordinates())
        {
            context.setPositions(metadata.coordinates());
        }

        if (metadata.hasVelocities())
        {
            context.setVelocities(metadata.velocities());
        }

        if (metadata.hasBoxVectors())
        {
            const auto boxvecs = metadata.boxVectors();

            context.setPeriodicBoxVectors(boxvecs[0], boxvecs[1], boxvecs[2]);
        }
    }

    inline void _populate_coords(QVector<SireMaths::Vector> &coords,
                                 const OpenMM::Vec3 *omm_coords,
                                 int natoms)
    {
        coords.resize(natoms);
        auto coords_data = coords.data();

        const double nm_to_internal = (1 * SireUnits::nanometer).to(SireUnits::angstrom);

        for (int i = 0; i < natoms; ++i)
        {
            const auto &omm = omm_coords[i];
            coords_data[i] = SireMaths::Vector(omm[0] * nm_to_internal,
                                               omm[1] * nm_to_internal,
                                               omm[2] * nm_to_internal);
        }
    }

    inline void _populate_vels(QVector<SireMol::Velocity3D> &vels,
                               const OpenMM::Vec3 *omm_vels,
                               int natoms)
    {
        vels.resize(natoms);
        auto vels_data = vels.data();

        for (int i = 0; i < natoms; ++i)
        {
            const auto &omm = omm_vels[i];
            vels_data[i] = SireMol::Velocity3D(omm[0] * SireUnits::nanometers_per_ps,
                                               omm[1] * SireUnits::nanometers_per_ps,
                                               omm[2] * SireUnits::nanometers_per_ps);
        }
    }

    SireVol::SpacePtr extract_space(const OpenMM::State &state)
    {
        OpenMM::Vec3 a, b, c;

        try
        {
            state.getPeriodicBoxVectors(a, b, c);
        }
        catch (...)
        {
            return SireVol::SpacePtr(SireVol::Cartesian());
        }

        const double nm_to_internal = (1 * SireUnits::nanometer).to(SireUnits::angstrom);

        SireMaths::Vector x(a[0] * nm_to_internal,
                            a[1] * nm_to_internal,
                            a[2] * nm_to_internal);

        SireMaths::Vector y(b[0] * nm_to_internal,
                            b[1] * nm_to_internal,
                            b[2] * nm_to_internal);

        SireMaths::Vector z(c[0] * nm_to_internal,
                            c[1] * nm_to_internal,
                            c[2] * nm_to_internal);

        SireVol::TriclinicBox triclinic;

        try
        {
            triclinic = SireVol::TriclinicBox(x, y, z);
        }
        catch (...)
        {
            // this is not a valid space - could be an infinite space
            return SireVol::SpacePtr(SireVol::Cartesian());
        }

        if (triclinic.alpha() == 90 and triclinic.beta() == 90 and triclinic.gamma() == 90)
        {
            // this is a PeriodicBox?
            SireVol::PeriodicBox pbox(x.max(y).max(z) - x.min(y).min(z));

            if (std::abs(pbox.volume().value() - triclinic.volume().value()) < 0.001)
            {
                // yes - periodic box
                return SireVol::SpacePtr(pbox);
            }
        }

        return SireVol::SpacePtr(triclinic);
    }

    void set_context_platform_property(OpenMM::Context &context,
                                       const QString &key,
                                       const QString &value)
    {
        OpenMM::Platform &platform = context.getPlatform();

        platform.setPropertyValue(context,
                                  key.toStdString(),
                                  value.toStdString());

        QString new_value = QString::fromStdString(platform.getPropertyValue(context, key.toStdString()));

        if (new_value != value)
            throw SireError::incompatible_error(QObject::tr(
                                                    "Unable to change the value of property %1 to `%2` in the "
                                                    "platform %3. The property value is still '%4'.")
                                                    .arg(key)
                                                    .arg(value)
                                                    .arg(QString::fromStdString(platform.getName()))
                                                    .arg(new_value),
                                                CODELOC);
    }

    SelectorMol extract_coordinates(const OpenMM::State &state,
                                    const SireMol::SelectorMol &mols,
                                    const QHash<SireMol::MolNum, SireBase::PropertyMap> &perturbable_maps,
                                    const SireBase::PropertyMap &map)
    {
        const auto positions = state.getPositions();

        const int natoms = positions.size();
        const auto positions_data = positions.data();

        if (mols.nAtoms() != natoms)
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "Different number of atoms from OpenMM and sire. "
                                                    "%1 versus %2. Cannot extract the coordinates.")
                                                    .arg(natoms)
                                                    .arg(mols.nAtoms()),
                                                CODELOC);
        }

        const auto coords_prop = map["coordinates"];

        const int nmols = mols.count();

        QVector<Molecule> ret(nmols);
        auto ret_data = ret.data();

        QVector<int> offsets(nmols);

        int offset = 0;

        for (int i = 0; i < nmols; ++i)
        {
            offsets[i] = offset;
            offset += mols[i].nAtoms();
        }

        const auto offsets_data = offsets.constData();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                QVector<SireMaths::Vector> converted_coords;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    auto mol = mols[i].edit();
                    const int mol_natoms = mol.nAtoms();

                    _populate_coords(converted_coords, positions_data+offsets_data[i], mol_natoms);

                    auto my_coords_prop = coords_prop.source();

                    if (perturbable_maps.contains(mol.number()))
                    {
                        my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                    }

                    if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                        converted_coords, false))
                    {
                        SireMol::AtomCoords c(mol.data().info());
                        c.copyFrom(converted_coords);
                        mol.setProperty(my_coords_prop, c);
                    }

                    ret_data[i] = mol.commit();
                } });
        }
        else
        {
            QVector<SireMaths::Vector> converted_coords;

            for (int i = 0; i < nmols; ++i)
            {
                auto mol = mols[i].edit();
                const int mol_natoms = mol.nAtoms();

                _populate_coords(converted_coords, positions_data + offsets_data[i], mol_natoms);

                auto my_coords_prop = coords_prop.source();

                if (perturbable_maps.contains(mol.number()))
                {
                    my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                }

                if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                    converted_coords, false))
                {
                    SireMol::AtomCoords c(mol.data().info());
                    c.copyFrom(converted_coords);
                    mol.setProperty(my_coords_prop, c);
                }

                ret_data[i] = mol.commit();
            }
        }

        return SelectorMol(ret);
    }

    SelectorMol extract_coordinates_and_velocities(const OpenMM::State &state,
                                                   const SireMol::SelectorMol &mols,
                                                   const QHash<SireMol::MolNum, SireBase::PropertyMap> &perturbable_maps,
                                                   const SireBase::PropertyMap &map)
    {
        const auto positions = state.getPositions();
        const auto velocities = state.getVelocities();

        const int natoms = positions.size();
        const auto positions_data = positions.data();
        const auto velocities_data = velocities.data();

        if (mols.nAtoms() != natoms)
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "Different number of atoms from OpenMM and sire. "
                                                    "%1 versus %2. Cannot extract the coordinates.")
                                                    .arg(natoms)
                                                    .arg(mols.nAtoms()),
                                                CODELOC);
        }

        const auto coords_prop = map["coordinates"];
        const auto vels_prop = map["velocity"];

        const int nmols = mols.count();

        QVector<Molecule> ret(nmols);
        auto ret_data = ret.data();

        QVector<int> offsets(nmols);

        int offset = 0;

        for (int i = 0; i < nmols; ++i)
        {
            offsets[i] = offset;
            offset += mols[i].nAtoms();
        }

        const auto offsets_data = offsets.constData();

        if (use_parallel(nmols, map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](const tbb::blocked_range<int> &r)
                              {
                QVector<SireMaths::Vector> converted_coords;
                QVector<SireMol::Velocity3D> converted_vels;

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    auto mol = mols[i].edit();
                    const int mol_natoms = mol.nAtoms();

                    auto my_coords_prop = coords_prop.source();
                    auto my_vels_prop = vels_prop.source();

                    if (perturbable_maps.contains(mol.number()))
                    {
                        my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                        my_vels_prop = perturbable_maps[mol.number()]["velocity"].source();
                    }

                    _populate_coords(converted_coords, positions_data+offsets_data[i], mol_natoms);
                    _populate_vels(converted_vels, velocities_data+offsets_data[i], mol_natoms);

                    if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                        converted_coords, false))
                    {
                        SireMol::AtomCoords c(mol.data().info());
                        c.copyFrom(converted_coords);
                        mol.setProperty(my_coords_prop, c);
                    }

                    if (not mol.updatePropertyFrom<SireMol::AtomVelocities>(my_vels_prop,
                                                                            converted_vels, false))
                    {
                        SireMol::AtomVelocities v(mol.data().info());
                        v.copyFrom(converted_vels);
                        mol.setProperty(my_vels_prop, v);
                    }

                    ret_data[i] = mol.commit();
                } });
        }
        else
        {
            QVector<SireMaths::Vector> converted_coords;
            QVector<SireMol::Velocity3D> converted_vels;

            for (int i = 0; i < nmols; ++i)
            {
                auto mol = mols[i].edit();
                const int mol_natoms = mol.nAtoms();

                auto my_coords_prop = coords_prop.source();
                auto my_vels_prop = vels_prop.source();

                if (perturbable_maps.contains(mol.number()))
                {
                    my_coords_prop = perturbable_maps[mol.number()]["coordinates"].source();
                    my_vels_prop = perturbable_maps[mol.number()]["velocity"].source();
                }

                _populate_coords(converted_coords, positions_data + offsets_data[i], mol_natoms);
                _populate_vels(converted_vels, velocities_data + offsets_data[i], mol_natoms);

                if (not mol.updatePropertyFrom<SireMol::AtomCoords>(my_coords_prop,
                                                                    converted_coords, false))
                {
                    SireMol::AtomCoords c(mol.data().info());
                    c.copyFrom(converted_coords);
                    mol.setProperty(my_coords_prop, c);
                }

                if (not mol.updatePropertyFrom<SireMol::AtomVelocities>(my_vels_prop,
                                                                        converted_vels, false))
                {
                    SireMol::AtomVelocities v(mol.data().info());
                    v.copyFrom(converted_vels);
                    mol.setProperty(my_vels_prop, v);
                }

                ret_data[i] = mol.commit();
            }
        }

        return SelectorMol(ret);
    }

} // end of namespace SireOpenMM
