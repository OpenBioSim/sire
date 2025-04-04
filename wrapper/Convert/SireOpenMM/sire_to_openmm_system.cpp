
#include "sire_openmm.h"

#include <OpenMM.h>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireSystem/forcefieldinfo.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/atommasses.h"
#include "SireMol/atomproperty.hpp"
#include "SireMol/atomvelocities.h"
#include "SireMol/bondid.h"
#include "SireMol/bondorder.h"
#include "SireMol/connectivity.h"
#include "SireMol/core.h"
#include "SireMol/moleditor.h"

#include "SireMM/amberparams.h"
#include "SireMM/anglerestraints.h"
#include "SireMM/atomljs.h"
#include "SireMM/bondrestraints.h"
#include "SireMM/inversebondrestraints.h"
#include "SireMM/boreschrestraints.h"
#include "SireMM/dihedralrestraints.h"
#include "SireMM/positionalrestraints.h"
#include "SireMM/selectorbond.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireCAS/lambdaschedule.h"

#include "SireMaths/maths.h"
#include "SireMaths/vector.h"

#include "SireBase/generalunitproperty.h"
#include "SireBase/lengthproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "tostring.h"

#include "openmmmolecule.h"
#include "pyqm.h"

#include <QDebug>

using SireBase::PropertyMap;
using SireCAS::LambdaSchedule;
using SireMol::Molecule;
using SireMol::MolNum;
using SireMol::SelectorMol;
using SireSystem::ForceFieldInfo;

using namespace SireOpenMM;

/** Add all of the Boresch restraints from 'restraints' to the
 *  passed systen, which is being acted on by the passed LambdaLever.
 *  The number of real (non-anchor) atoms in the OpenMM::System is 'natoms'
 */
void _add_boresch_restraints(const SireMM::BoreschRestraints &restraints,
                             OpenMM::System &system, LambdaLever &lambda_lever,
                             int natoms)
{
    if (restraints.isEmpty())
        return;

    // energy expression of the Boresch restraint, which acts over
    // six atoms - this is a set of
    // one distance restraint, two angle restraints and three
    // torsion restraints between the three receptor (particles 0 to 2)
    // and three ligand atoms (particles 3 to 5)
    //
    // r is | Ligand1 - Receptor1 | = distance(P1, P4)
    // thetaA = angle(R2, R1, L1) = angle(P2, P1, P4)
    // thetaB = angle(R1, L1, L2) = angle(P1, P4, P5)
    // phiA = dihedral(R3, R2, R1, L1) = dihedral(P3, P2, P1, P4)
    // phiB = dihedral(R2, R1, L1, L2) = dihedral(P2, P1, P4, P5)
    // phiC = dihedral(R1, L1, L2, L3) = dihedral(P1, P4, P5, P6)
    //
    // Then the energies are
    //
    // e_restraint = rho * (e_bond + e_angle + e_torsion)
    // e_bond = kr (r - r0)^2
    // e_angle_i = ktheta_i (theta_i - theta0_i)^2
    // e_torsion_i = k_phi_i (min(dphi_i, 2pi-dphi_i))^2 where
    // dphi_i = abs(phi_i - phi0_i)
    //
    const auto energy_expression = QString(
                                       "rho*(e_bond + e_angle_A + e_angle_B + e_torsion_A + e_torsion_B + e_torsion_C);"
                                       "e_bond=kr*(r-r0)^2;"
                                       "e_angle_B=ktheta_B*(theta_B-theta0_B)^2;"
                                       "e_angle_A=ktheta_A*(theta_A-theta0_A)^2;"
                                       "e_torsion_C=kphi_C*(min(dphi_C, two_pi-dphi_C))^2;"
                                       "e_torsion_B=kphi_B*(min(dphi_B, two_pi-dphi_B))^2;"
                                       "e_torsion_A=kphi_A*(min(dphi_A, two_pi-dphi_A))^2;"
                                       "dphi_C=abs(phi_C-phi0_C);"
                                       "dphi_B=abs(phi_B-phi0_B);"
                                       "dphi_A=abs(phi_A-phi0_A);"
                                       "two_pi=6.283185307179586;"
                                       "phi_C=dihedral(p1, p4, p5, p6);"
                                       "phi_B=dihedral(p2, p1, p4, p5);"
                                       "phi_A=dihedral(p3, p2, p1, p4);"
                                       "theta_B=angle(p1, p4, p5);"
                                       "theta_A=angle(p2, p1, p4);"
                                       "r=distance(p1, p4);")
                                       .toStdString();

    auto *restraintff = new OpenMM::CustomCompoundBondForce(6, energy_expression);
    restraintff->setName("BoreschRestraintForce");

    restraintff->addPerBondParameter("rho");
    restraintff->addPerBondParameter("kr");
    restraintff->addPerBondParameter("r0");
    restraintff->addPerBondParameter("ktheta_A");
    restraintff->addPerBondParameter("theta0_A");
    restraintff->addPerBondParameter("ktheta_B");
    restraintff->addPerBondParameter("theta0_B");
    restraintff->addPerBondParameter("kphi_A");
    restraintff->addPerBondParameter("phi0_A");
    restraintff->addPerBondParameter("kphi_B");
    restraintff->addPerBondParameter("phi0_B");
    restraintff->addPerBondParameter("kphi_C");
    restraintff->addPerBondParameter("phi0_C");

    restraintff->setUsesPeriodicBoundaryConditions(true);

    lambda_lever.addRestraintIndex(restraints.name(),
                                   system.addForce(restraintff));

    const double internal_to_nm = (1 * SireUnits::angstrom).to(SireUnits::nanometer);
    const double internal_to_k = (1 * SireUnits::kcal_per_mol / (SireUnits::angstrom2)).to(SireUnits::kJ_per_mol / SireUnits::nanometer2);
    const double internal_to_ktheta = (1 * SireUnits::kcal_per_mol / (SireUnits::radian2)).to(SireUnits::kJ_per_mol / SireUnits::radian2);

    for (const auto &restraint : restraints.restraints())
    {
        if (restraint.isNull())
            continue;

        std::vector<int> particles;
        particles.resize(6);

        std::vector<double> parameters;
        parameters.resize(13);

        for (int i = 0; i < 3; ++i)
        {
            particles[i] = restraint.receptorAtoms()[i];
            particles[i + 3] = restraint.ligandAtoms()[i];

            if (particles[i] < 0 or particles[i] >= natoms or
                particles[i + 3] < 0 or particles[i + 3] >= natoms)
            {
                throw SireError::invalid_index(QObject::tr(
                                                   "Invalid particle indicies (ligand=%1, receptor=2) "
                                                   "for a Boresch restraint for %3 atoms.")
                                                   .arg(Sire::toString(restraint.ligandAtoms()))
                                                   .arg(Sire::toString(restraint.receptorAtoms()))
                                                   .arg(natoms),
                                               CODELOC);
            }
        }

        parameters[0] = 1.0;                                                // rho
        parameters[1] = restraint.kr().value() * internal_to_k;             // kr
        parameters[2] = restraint.r0().value() * internal_to_nm;            // r0
        parameters[3] = restraint.ktheta()[0].value() * internal_to_ktheta; // ktheta_A
        parameters[4] = restraint.theta0()[0].value();                      // theta0_A (already in radians)
        parameters[5] = restraint.ktheta()[1].value() * internal_to_ktheta; // ktheta_B
        parameters[6] = restraint.theta0()[1].value();                      // theta0_B
        parameters[7] = restraint.kphi()[0].value() * internal_to_ktheta;   // kphi_A
        parameters[8] = restraint.phi0()[0].value();                        // phi0_A
        parameters[9] = restraint.kphi()[1].value() * internal_to_ktheta;   // kphi_B
        parameters[10] = restraint.phi0()[1].value();                       // phi0_B
        parameters[11] = restraint.kphi()[2].value() * internal_to_ktheta;  // kphi_C
        parameters[12] = restraint.phi0()[2].value();                       // phi0_C

        restraintff->addBond(particles, parameters);
    }
}

/** Add all of the bond restraints from 'restraints' to the passed
 *  system, which is acted on by the passed LambdaLever. The number
 *  of real (non-anchor) atoms in the OpenMM::System is 'natoms'
 */
void _add_bond_restraints(const SireMM::BondRestraints &restraints,
                          OpenMM::System &system, LambdaLever &lambda_lever,
                          int natoms)
{
    if (restraints.isEmpty())
        return;

    if (restraints.hasCentroidRestraints())
    {
        throw SireError::unsupported(QObject::tr(
                                         "Centroid bond restraints aren't yet supported..."),
                                     CODELOC);
    }

    // energy expression of a harmonic bond potential, scaled by rho
    const auto energy_expression = QString(
                                       "rho*k*delta*delta;"
                                       "delta=(r-r0)")
                                       .toStdString();

    auto *restraintff = new OpenMM::CustomBondForce(energy_expression);
    restraintff->setName("BondRestraintForce");

    restraintff->addPerBondParameter("rho");
    restraintff->addPerBondParameter("k");
    restraintff->addPerBondParameter("r0");

    restraintff->setUsesPeriodicBoundaryConditions(true);

    lambda_lever.addRestraintIndex(restraints.name(),
                                   system.addForce(restraintff));

    const auto atom_restraints = restraints.atomRestraints();

    const double internal_to_nm = (1 * SireUnits::angstrom).to(SireUnits::nanometer);
    const double internal_to_k = (1 * SireUnits::kcal_per_mol / (SireUnits::angstrom2)).to(SireUnits::kJ_per_mol / (SireUnits::nanometer2));

    auto cljff = lambda_lever.getForce<OpenMM::NonbondedForce>("clj", system);

    std::vector<double> custom_params = {1.0, 0.0, 0.0};

    for (const auto &restraint : atom_restraints)
    {
        int atom0_index = restraint.atom0();
        int atom1_index = restraint.atom1();

        if (atom0_index < 0 or atom0_index >= natoms)
            throw SireError::invalid_index(QObject::tr(
                                               "Invalid particle index! %1 from %2")
                                               .arg(atom0_index)
                                               .arg(natoms),
                                           CODELOC);

        if (atom1_index < 0 or atom1_index >= natoms)
            throw SireError::invalid_index(QObject::tr(
                                               "Invalid particle index! %1 from %2")
                                               .arg(atom1_index)
                                               .arg(natoms),
                                           CODELOC);

        custom_params[0] = 1.0;                                     // rho - always equal to 1 (scaled by lever)
        custom_params[1] = restraint.k().value() * internal_to_k;   // k
        custom_params[2] = restraint.r0().value() * internal_to_nm; // rb

        restraintff->addBond(atom0_index, atom1_index, custom_params);
    }
}

void _add_inverse_bond_restraints(const SireMM::InverseBondRestraints &restraints,
                                OpenMM::System &system, LambdaLever &lambda_lever,
                                int natoms)
{
    if (restraints.isEmpty())
        return;

    if (restraints.hasCentroidRestraints())
    {
        throw SireError::unsupported(QObject::tr(
                                         "Centroid bond restraints aren't yet supported..."),
                                     CODELOC);
    }

    const auto energy_expression = QString(
        "rho*k*delta*delta*step;"
        "delta=(r-r0);"
        "step=max(0,min(1,(r0 - r)))")
        .toStdString();

    auto *restraintff = new OpenMM::CustomBondForce(energy_expression);
    restraintff->setName("InverseBondRestraintForce");

    restraintff->addPerBondParameter("rho");
    restraintff->addPerBondParameter("k");
    restraintff->addPerBondParameter("r0");

    restraintff->setUsesPeriodicBoundaryConditions(true);

    lambda_lever.addRestraintIndex(restraints.name(),
    system.addForce(restraintff));

    const auto atom_restraints = restraints.atomRestraints();

    const double internal_to_nm = (1 * SireUnits::angstrom).to(SireUnits::nanometer);
    const double internal_to_k = (1 * SireUnits::kcal_per_mol / (SireUnits::angstrom2)).to(SireUnits::kJ_per_mol / (SireUnits::nanometer2));

    auto cljff = lambda_lever.getForce<OpenMM::NonbondedForce>("clj", system);

    std::vector<double> custom_params = {1.0, 0.0, 0.0};

    for (const auto &restraint : atom_restraints)
    {
        int atom0_index = restraint.atom0();
        int atom1_index = restraint.atom1();

        if (atom0_index < 0 or atom0_index >= natoms)
        throw SireError::invalid_index(QObject::tr(
                        "Invalid particle index! %1 from %2")
                        .arg(atom0_index)
                        .arg(natoms),
                    CODELOC);

        if (atom1_index < 0 or atom1_index >= natoms)
        throw SireError::invalid_index(QObject::tr(
                        "Invalid particle index! %1 from %2")
                        .arg(atom1_index)
                        .arg(natoms),
                    CODELOC);

        custom_params[0] = 1.0;                                     // rho - always equal to 1 (scaled by lever)
        custom_params[1] = restraint.k().value() * internal_to_k;   // k
        custom_params[2] = restraint.r0().value() * internal_to_nm; // rb

        restraintff->addBond(atom0_index, atom1_index, custom_params);
    }
}
/** Add all of the positional restraints from 'restraints' to the passed
 *  system, which is acted on by the passed LambdaLever. All of the
 *  existing anchor atoms are in 'anchor_coords', which this function
 *  will add to if any more need adding. The number of real (non-anchor)
 *  atoms in the OpenMM::System is 'natoms'
 */
void _add_positional_restraints(const SireMM::PositionalRestraints &restraints,
                                OpenMM::System &system, LambdaLever &lambda_lever,
                                std::vector<OpenMM::Vec3> &anchor_coords,
                                int natoms)
{
    if (restraints.isEmpty())
        return;

    if (restraints.hasCentroidRestraints())
    {
        throw SireError::unsupported(QObject::tr(
                                         "Centroid positional restraints aren't yet supported..."),
                                     CODELOC);
    }

    // energy expression of a flat-bottom well potential, scaled by rho
    const auto energy_expression = QString(
                                       "rho*k*step(delta)*delta*delta;"
                                       "delta=(r-rb)")
                                       .toStdString();

    auto *restraintff = new OpenMM::CustomBondForce(energy_expression);
    restraintff->setName("PositionalRestraintForce");

    restraintff->addPerBondParameter("rho");
    restraintff->addPerBondParameter("k");
    restraintff->addPerBondParameter("rb");

    restraintff->setUsesPeriodicBoundaryConditions(true);

    lambda_lever.addRestraintIndex(restraints.name(),
                                   system.addForce(restraintff));

    const auto atom_restraints = restraints.atomRestraints();

    anchor_coords.reserve(anchor_coords.size() + atom_restraints.count());

    auto to_vec3 = [](const SireMaths::Vector &coords)
    {
        const double internal_to_nm = (1 * SireUnits::angstrom).to(SireUnits::nanometer);

        return OpenMM::Vec3(internal_to_nm * coords.x(),
                            internal_to_nm * coords.y(),
                            internal_to_nm * coords.z());
    };

    QVector<SireMaths::Vector> anchor_idxs;

    const double internal_to_nm = (1 * SireUnits::angstrom).to(SireUnits::nanometer);
    const double internal_to_k = (1 * SireUnits::kcal_per_mol / (SireUnits::angstrom2)).to(SireUnits::kJ_per_mol / (SireUnits::nanometer2));

    auto cljff = lambda_lever.getForce<OpenMM::NonbondedForce>("clj", system);
    auto ghost_ghostff = lambda_lever.getForce<OpenMM::CustomNonbondedForce>("ghost/ghost", system);
    auto ghost_nonghostff = lambda_lever.getForce<OpenMM::CustomNonbondedForce>("ghost/non-ghost", system);

    std::vector<double> custom_params = {1.0, 0.0, 0.0};
    // Define null parameters used to add these particles to the ghost forces (5 total)
    std::vector<double> custom_clj_params = {0.0, 0.0, 0.0, 0.0, 0.0};

    // we need to add all of the positions as anchor particles
    for (const auto &restraint : atom_restraints)
    {
        int atom_index = restraint.atom();

        if (atom_index < 0 or atom_index >= natoms)
            throw SireError::invalid_index(QObject::tr(
                                               "Invalid particle index! %1 from %2")
                                               .arg(atom_index)
                                               .arg(natoms),
                                           CODELOC);

        // find the anchor at this position
        int anchor_index = anchor_idxs.indexOf(restraint.position());

        if (anchor_index != -1)
        {
            anchor_index += natoms;
        }
        else if (anchor_index == -1)
        {
            // doesn't exist - create it as a massless particle
            // (won't be moved)
            anchor_index = system.addParticle(0.0);
            anchor_coords.push_back(to_vec3(restraint.position()));
            anchor_idxs.append(restraint.position());

            // add a null particle to the nonbonded forces
            if (cljff != 0)
            {
                cljff->addParticle(0, 0, 0);
            }

            if (ghost_ghostff != 0)
            {
                ghost_ghostff->addParticle(custom_clj_params);
            }

            if (ghost_nonghostff != 0)
            {
                ghost_nonghostff->addParticle(custom_clj_params);
            }
        }

        custom_params[0] = 1.0;                                     // rho - always equal to 1 (scaled by lever)
        custom_params[1] = restraint.k().value() * internal_to_k;   // k
        custom_params[2] = restraint.r0().value() * internal_to_nm; // rb

        if (cljff != 0)
        {
            // make sure to exclude interactions between
            // the atom being positionally restrained and
            // the anchor
            cljff->addException(anchor_index, atom_index, 0, 0, 0, true);
        }

        if (ghost_ghostff != 0)
        {
            // make sure to exclude interactions between
            // the atom being positionally restrained and
            // the anchor
            ghost_ghostff->addExclusion(anchor_index, atom_index);
        }

        if (ghost_nonghostff != 0)
        {
            // make sure to exclude interactions between
            // the atom being positionally restrained and
            // the anchor
            ghost_nonghostff->addExclusion(anchor_index, atom_index);
        }

        restraintff->addBond(anchor_index, atom_index, custom_params);
    }
}

/** Add all of the angle restraints from 'restraints' to the passed
 *  system, which is acted on by the passed LambdaLever. The number
 *  of real (non-anchor) atoms in the OpenMM::System is 'natoms'
 */

void _add_angle_restraints(const SireMM::AngleRestraints &restraints,
                           OpenMM::System &system, LambdaLever &lambda_lever,
                           int natoms)
{
    if (restraints.isEmpty())
        return;

    // energy expression of the angle restraint, which acts over three atoms
    //
    // theta = angle(P1, P2, P3)
    //
    // The energies are
    //
    // e_restraint = rho * (e_angle)
    // e_angle = ktheta(theta - theta0)^2

    const auto energy_expression = QString(
                                       "rho*k*(theta-theta0)^2;")
                                       .toStdString();

    auto *restraintff = new OpenMM::CustomAngleForce(energy_expression);

    restraintff->setName("AngleRestraintForce");
    restraintff->addPerAngleParameter("rho");
    restraintff->addPerAngleParameter("k");
    restraintff->addPerAngleParameter("theta0");

    restraintff->setUsesPeriodicBoundaryConditions(true);

    lambda_lever.addRestraintIndex(restraints.name(),
                                   system.addForce(restraintff));

    const double internal_to_ktheta = (1 * SireUnits::kcal_per_mol / (SireUnits::radian2)).to(SireUnits::kJ_per_mol / SireUnits::radian2);

    const auto atom_restraints = restraints.restraints();

    for (const auto &restraint : atom_restraints)
    {
        std::vector<int> particles;
        particles.resize(3);

        for (int i = 0; i < 3; ++i)
        {
            particles[i] = restraint.atoms()[i];
        }

        std::vector<double> parameters;
        parameters.resize(3);

        parameters[0] = 1.0;                                             // rho
        parameters[1] = restraint.ktheta().value() * internal_to_ktheta; // k
        parameters[2] = restraint.theta0().value();                      // theta0 (already in radians)

        // restraintff->addTorsion(particles, parameters);
        restraintff->addAngle(particles[0], particles[1], particles[2], parameters);
    }
}

void _add_dihedral_restraints(const SireMM::DihedralRestraints &restraints,
                              OpenMM::System &system, LambdaLever &lambda_lever,
                              int natoms)
{
    if (restraints.isEmpty())
        return;

    // energy expression of the dihedral restraint, which acts over four atoms
    //
    // phi = dihedral(P1, P2, P3, P4)
    //
    // The energies are
    //
    // e_restraint = rho * (e_torsion)
    // e_torsion = k_phi(dphi)^2 where
    // dphi = abs(phi - phi0)

    const auto energy_expression = QString(
                                       "rho*k*min(dtheta, two_pi-dtheta)^2;"
                                       "dtheta = abs(theta-theta0);"
                                       "two_pi=6.283185307179586;")
                                       .toStdString();

    auto *restraintff = new OpenMM::CustomTorsionForce(energy_expression);

    // it seems that OpenMM wants to call the torsion angle theta rather than phi
    // we need to rename our parameters accordingly
    restraintff->setName("TorsionRestraintForce");
    restraintff->addPerTorsionParameter("rho");
    restraintff->addPerTorsionParameter("k");
    restraintff->addPerTorsionParameter("theta0");

    restraintff->setUsesPeriodicBoundaryConditions(true);

    lambda_lever.addRestraintIndex(restraints.name(),
                                   system.addForce(restraintff));

    const double internal_to_ktheta = (1 * SireUnits::kcal_per_mol / (SireUnits::radian2)).to(SireUnits::kJ_per_mol / SireUnits::radian2);

    const auto atom_restraints = restraints.restraints();

    for (const auto &restraint : atom_restraints)
    {
        std::vector<int> particles;
        particles.resize(4);

        for (int i = 0; i < 4; ++i)
        {
            particles[i] = restraint.atoms()[i];
        }

        std::vector<double> parameters;
        parameters.resize(3);

        parameters[0] = 1.0;                                           // rho
        parameters[1] = restraint.kphi().value() * internal_to_ktheta; // k
        parameters[2] = restraint.phi0().value();                      // phi0 (already in radians) --> theta0

        // restraintff->addTorsion(particles, parameters);
        restraintff->addTorsion(particles[0], particles[1], particles[2], particles[3], parameters);
    }
}

/** Set the coulomb and LJ cutoff in the passed NonbondedForce,
 *  based on the information in the passed ForceFieldInfo.
 *  This sets the cutoff type (e.g. PME) and the actual
 *  cutoff length (if one is used)
 */
void _set_clj_cutoff(OpenMM::NonbondedForce &cljff,
                     const ForceFieldInfo &ffinfo)
{
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

            cljff.setNonbondedMethod(nbmethod);

            double tolerance = ffinfo.getParameter("tolerance").value();

            if (tolerance <= 0)
                tolerance = 0.001;

            cljff.setEwaldErrorTolerance(tolerance);
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

            cljff.setNonbondedMethod(nbmethod);
            cljff.setReactionFieldDielectric(dielectric);
        }
        else if (typ == "CUTOFF")
        {
            // use reaction field for non-periodic spaces, and PME for periodic
            if (ffinfo.space().isPeriodic())
            {
                const auto nbmethod = OpenMM::NonbondedForce::PME;
                cljff.setNonbondedMethod(nbmethod);

                double tolerance = ffinfo.getParameter("tolerance").value();

                if (tolerance <= 0)
                    tolerance = 0.001;

                cljff.setEwaldErrorTolerance(tolerance);
            }
            else
            {
                const auto nbmethod = OpenMM::NonbondedForce::CutoffNonPeriodic;
                cljff.setNonbondedMethod(nbmethod);

                double dielectric = ffinfo.getParameter("dielectric").value();

                if (dielectric <= 0)
                    dielectric = 78.3;

                cljff.setReactionFieldDielectric(dielectric);
            }
        }

        const auto cutoff = ffinfo.cutoff().to(SireUnits::nanometers);
        cljff.setCutoffDistance(cutoff);
    }
    else
    {
        cljff.setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);
    }
}

/** Set the periodic space box vectors in the system, returning
 *  them if they are set (we don't do anything for non-periodic
 *  spaces)
 */
std::shared_ptr<std::vector<OpenMM::Vec3>>
_set_box_vectors(OpenMM::System &system,
                 const ForceFieldInfo &ffinfo)
{
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

    return boxvecs;
}

class IndexPair
{
public:
    IndexPair(int atom0 = 0, int atom1 = 0) : _atom0(atom0), _atom1(atom1)
    {
        if (atom1 < atom0)
        {
            std::swap(_atom0, _atom1);
        }
    }

    bool operator==(const IndexPair &other) const
    {
        return _atom0 == other._atom0 and _atom1 == other._atom1;
    }

    int _atom0;
    int _atom1;
};

uint qHash(const IndexPair &pair)
{
    return qHash(pair._atom0) ^ qHash(pair._atom1);
}

/**

This is the (monster) function that converts a passed set of Sire
molecules (in the passed SelectorMols) into an OpenMM::System,
controlled via the properties in the passed PropertyMap.

The OpenMM::System is constructed in the passed (empty)
OpenMM::System that is the first argument. This is because this
function is called from Python, and this was the only way found to
have the resulting OpenMM::System make its way back up to the
Python layer.

This returns an extra set of metadata that doesn't fit into the
OpenMM::System. This metadata includes information about any
perturbations, the atom index, plus the coordinates and velocities
from the molecules (if these could be found)

This is a monster function, as it does need to do everything, and
the parts of not easily decomposable (they need information from
a prior part that could be passed as function arguments, but
would be messy).

This function is best read as a sequence of stages. These stages
are commented within the function. The stages are:

1. Initialisation - copying molecular data from sire into OpenMMMolecule

2. Create base forces (forcefields)

3. Define the LambdaLever and LambdaSchedule

4. Define the forces (forcefields) for the ghost atoms

5. Copy atomistic forcefield parameters to the OpenMM forces

6. Set up the nonbonded pair exceptions

7. Set up the restraints

8. Copy across all of the coordinates and velocities

*/
OpenMMMetaData SireOpenMM::sire_to_openmm_system(OpenMM::System &system,
                                                 const SelectorMol &mols,
                                                 const PropertyMap &map)
{
    ///
    /// Stage 1 - initialisation.
    ///
    /// Aim is to copy all of the molecular information from the
    /// sire SelectorMol and put it into a set of OpenMMMolecule
    /// objects ready for that information to then be fed into the
    /// OpenMM::System.
    ///
    /// We assume that a valid, empty OpenMM::System has been
    /// passed to us
    ///

    // We get the forcefield parameters from the map (e.g. cutoffs etc)
    ForceFieldInfo ffinfo(mols, map);

    // Some initialisation - do we have anything to convert?
    const int nmols = mols.count();

    if (nmols == 0)
    {
        // nothing to do
        return OpenMMMetaData();
    }

    // set the box vectors for periodic spaces
    auto boxvecs = _set_box_vectors(system, ffinfo);

    // Are any of the molecules perturbable, and if they are,
    // should we just ignore perturbations?
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

    // Create a vector containing the start index for the atoms in each molecule.
    QVector<int> start_atom_index(nmols);
    start_atom_index[0] = 0;
    for (int i = 1; i < nmols; ++i)
    {
        start_atom_index[i] = start_atom_index[i-1] + mols[i-1].nAtoms();
    }

    if (SireBase::should_run_in_parallel(nmols, map))
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, mols.count()), [&](const tbb::blocked_range<int> &r)
                          {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    openmm_mols_data[i] = OpenMMMolecule(mols[i], start_atom_index[i], map);
                } });
    }
    else
    {
        for (int i = 0; i < nmols; ++i)
        {
            openmm_mols_data[i] = OpenMMMolecule(mols[i], start_atom_index[i], map);
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

    QSet<int> fixed_atoms;

    if (map.specified("fixed"))
    {
        // this should be a list of indexes of atoms to fix
        const auto idxs = map["fixed"].value().asA<SireBase::IntegerArrayProperty>().toVector();

        if (not idxs.isEmpty())
        {
            fixed_atoms.reserve(idxs.count());

            for (const auto &idx : idxs)
            {
                fixed_atoms.insert(idx);
            }
        }
    }

    // End of stage 1 - we have now extracted all of the molecular information
    // and have worked out what parameters to use and whether any of the
    // molecules are perturbable

    ///
    /// Stage 2 - Base Forces
    ///
    /// Aim is to create the base forces, e.g. Nonbonded, Bond, Angle,
    /// Torsion. These forces are created and parameterised via the
    /// PropertyMap and ForceFieldInfo (e.g. cutoff type etc)
    ///

    // Force for the Coulomb and LJ (CLJ) energy between
    // all non-perturbable atoms
    OpenMM::NonbondedForce *cljff = new OpenMM::NonbondedForce();

    bool use_dispersion_correction = false;

    if (map.specified("use_dispersion_correction"))
    {
        use_dispersion_correction = map["use_dispersion_correction"].value().asABoolean();
    }

    // note that this will be very slow for perturbable systems, as
    // it needs recalculating for every change of lambda
    cljff->setUseDispersionCorrection(use_dispersion_correction);

    // set the non-bonded cutoff type and length based on
    // the infomation in ffinfo
    _set_clj_cutoff(*cljff, ffinfo);

    // now create the base bond, angle and torsion forcefields
    OpenMM::HarmonicBondForce *bondff = new OpenMM::HarmonicBondForce();
    OpenMM::HarmonicAngleForce *angff = new OpenMM::HarmonicAngleForce();
    OpenMM::PeriodicTorsionForce *dihff = new OpenMM::PeriodicTorsionForce();

    // now create the engine for computing QM forces on atoms
    QMForce *qmff = 0;

    if (map.specified("qm_engine"))
    {
        try
        {
            auto &engine = map["qm_engine"].value().asA<QMEngine>();
            qmff = engine.createForce();
            qmff->setName("QMForce");
        }
        catch (...)
        {
            throw SireError::incompatible_error(QObject::tr("Invalid QM engine!"), CODELOC);
        }
    }

    // end of stage 2 - we now have the base forces

    ///
    /// Stage 3 - define the LambdaLever and LambdaSchedule
    ///
    /// Aim is to create the lambda lever that perturbs forcefield
    /// potentials from the reference state to the perturbed state.
    /// This lever exists even if there are no perturbable molecules,
    /// as it provides useful metadata about the OpenMM::System.
    ///

    // First, create the LambdaLever
    LambdaLever lambda_lever;

    // Let the user supply their own schedule for any perturbations
    if (map.specified("schedule"))
    {
        lambda_lever.setSchedule(
            map["schedule"].value().asA<LambdaSchedule>());
    }
    else if (any_perturbable)
    {
        // use a standard morph if we have an alchemical or QM perturbation
        lambda_lever.setSchedule(
            LambdaSchedule::standard_morph());
    }

    // Add any QM force first so that we can guarantee that it is index zero.
    if (qmff != 0)
    {
        lambda_lever.setForceIndex("qmff", system.addForce(qmff));
        lambda_lever.addLever("qm_scale");
    }

    // We can now add the standard forces to the OpenMM::System.
    // We do this here, so that we can capture the index of the
    // force and associate it with a name in the lever.
    lambda_lever.setForceIndex("clj", system.addForce(cljff));

    // We also want to name the levers available for this force,
    // e.g. we can change the charge, sigma and epsilon parameters
    lambda_lever.addLever("charge");
    lambda_lever.addLever("sigma");
    lambda_lever.addLever("epsilon");

    // and the exceptions
    lambda_lever.addLever("charge_scale");
    lambda_lever.addLever("lj_scale");

    // Do the same for the bond, angle and torsion forces
    lambda_lever.setForceIndex("bond", system.addForce(bondff));
    lambda_lever.addLever("bond_length");
    lambda_lever.addLever("bond_k");

    lambda_lever.setForceIndex("angle", system.addForce(angff));
    lambda_lever.addLever("angle_size");
    lambda_lever.addLever("angle_k");

    lambda_lever.setForceIndex("torsion", system.addForce(dihff));
    lambda_lever.addLever("torsion_phase");
    lambda_lever.addLever("torsion_k");

    ///
    /// Stage 4 - define the forces for ghost atoms
    ///
    /// We now need to create the forces that enable the use of a
    /// softening potential that smooths the interactions of
    /// any ghost atoms that are created or deleted. We skip
    /// this bit if there are no ghost atoms
    ///

    OpenMM::CustomBondForce *ghost_14ff = 0;
    OpenMM::CustomNonbondedForce *ghost_ghostff = 0;
    OpenMM::CustomNonbondedForce *ghost_nonghostff = 0;

    if (any_perturbable)
    {
        lambda_lever.addLever("alpha");
        lambda_lever.addLever("kappa");

        // somd uses a default shift_delta of 2.0 A
        SireUnits::Dimension::Length shift_delta = 2.0 * SireUnits::angstrom;

        if (map.specified("shift_delta"))
        {
            const auto &value = map["shift_delta"].value();

            if (value.isA<SireBase::LengthProperty>())
                shift_delta = value.asA<SireBase::LengthProperty>().value();
            else
                shift_delta = value.asA<SireBase::GeneralUnitProperty>().toUnit<SireUnits::Dimension::Length>();
        }

        if (shift_delta.value() < 0)
            shift_delta = 0.0 * SireUnits::angstrom;

        // somd uses a default shift_coulomb of 1 A
        double shift_coulomb = SireMaths::pow_2((1 * SireUnits::angstrom).to(SireUnits::nanometer));

        if (map.specified("shift_coulomb"))
        {
            const auto &value = map["shift_coulomb"].value();

            if (value.isA<SireBase::LengthProperty>())
                shift_coulomb = SireMaths::pow_2(value.asA<SireBase::LengthProperty>().value().to(SireUnits::nanometer));
            else
                shift_coulomb = SireMaths::pow_2(
                    value.asA<SireBase::GeneralUnitProperty>().toUnit<SireUnits::Dimension::Length>().to(SireUnits::nanometer));
        }

        // use a Taylor LJ power of 1
        int taylor_power = 1;

        if (map.specified("taylor_power"))
        {
            taylor_power = map["taylor_power"].value().asAnInteger();
        }

        if (taylor_power < 0)
            taylor_power = 0;
        else if (taylor_power > 4)
            taylor_power = 4;

        // by default we use zacharias softening
        bool use_taylor_softening = false;

        if (map.specified("use_taylor_softening"))
        {
            use_taylor_softening = map["use_taylor_softening"].value().asABoolean();
        }

        if (map.specified("use_zacharias_softening"))
        {
            use_taylor_softening = not map["use_zacharias_softening"].value().asABoolean();
        }

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

        auto taylor_power_expression = [](const QString &alpha, int power)
        {
            if (power == 0)
                return QString("1");
            else if (power == 1)
                return QString("%1").arg(alpha);
            else if (power == 2)
                return QString("%1*%1").arg(alpha);
            else
                return QString("%1^%2").arg(alpha).arg(power);
        };

        // see below for the description of this energy expression
        std::string nb14_expression, clj_expression;

        if (use_taylor_softening)
        {
            nb14_expression = QString(
                                  "coul_nrg+lj_nrg;"
                                  "coul_nrg=138.9354558466661*q*(((%1)/sqrt((%2*alpha*alpha)+r_safe^2))-(kappa/r_safe));"
                                  "lj_nrg=four_epsilon*sig6*(sig6-1);"
                                  "sig6=(sigma^6)/(%3*sigma^6 + r_safe^6);"
                                  "r_safe=max(r, 0.001);")
                                  .arg(coulomb_power_expression("alpha", coulomb_power))
                                  .arg(shift_coulomb)
                                  .arg(taylor_power_expression("alpha", taylor_power))
                                  .toStdString();
        }
        else
        {
            nb14_expression = QString(
                                  "coul_nrg+lj_nrg;"
                                  "coul_nrg=138.9354558466661*q*(((%1)/sqrt((%2*alpha*alpha)+r_safe^2))-(kappa/r_safe));"
                                  "lj_nrg=four_epsilon*sig6*(sig6-1);"
                                  "sig6=(sigma^6)/(((sigma*delta) + r_safe^2)^3);"
                                  "r_safe=max(r, 0.001);"
                                  "delta=%3*alpha;")
                                  .arg(coulomb_power_expression("alpha", coulomb_power))
                                  .arg(shift_coulomb)
                                  .arg(shift_delta.to(SireUnits::nanometer))
                                  .toStdString();
        }

        ghost_14ff = new OpenMM::CustomBondForce(nb14_expression);
        ghost_14ff->setName("Ghost14BondForce");

        ghost_14ff->addPerBondParameter("q");
        ghost_14ff->addPerBondParameter("sigma");
        ghost_14ff->addPerBondParameter("four_epsilon");
        ghost_14ff->addPerBondParameter("alpha");
        ghost_14ff->addPerBondParameter("kappa");

        // short-range intramolecular term that should not use
        // periodic boundaries or cutoffs
        ghost_14ff->setUsesPeriodicBoundaryConditions(false);

        if (use_taylor_softening)
        {
            // this uses the following potentials
            //            Zacharias and McCammon, J. Chem. Phys., 1994, and also,
            //            Michel et al., JCTC, 2007
            //            LJ is Rich Taylor's softcore LJ
            //
            //   V_{LJ}(r) = 4 epsilon [ (sigma^12 / (alpha^m sigma^6 + r^6)^2) -
            //                           (sigma^6  / (alpha^m sigma^6 + r^6) ) ]
            //
            //   V_{coul}(r) = (1-alpha)^n q_i q_j / 4 pi eps_0 (delta+r^2)^(1/2)
            //
            //   delta = shift_coulomb^2 * alpha
            //
            // Note that we supply half_sigma and two_sqrt_epsilon to save some
            // cycles
            //
            // Note also that we subtract the normal coulomb energy as this
            // is calculated during the standard NonbondedForce
            //
            // 138.9354558466661 is the constant needed to get energies in
            // kJ mol-1 given the units of charge (|e|) and distance (nm)
            //
            clj_expression = QString("coul_nrg+lj_nrg;"
                                     "coul_nrg=138.9354558466661*q1*q2*(((%1)/sqrt((%2*max_alpha*max_alpha)+r_safe^2))-(max_kappa/r_safe));"
                                     "lj_nrg=two_sqrt_epsilon1*two_sqrt_epsilon2*sig6*(sig6-1);"
                                     "sig6=(sigma^6)/(%3*sigma^6 + r_safe^6);"
                                     "r_safe=max(r, 0.001);"
                                     "max_kappa=max(kappa1, kappa2);"
                                     "max_alpha=max(alpha1, alpha2);"
                                     "sigma=half_sigma1+half_sigma2;")
                                 .arg(coulomb_power_expression("max_alpha", coulomb_power))
                                 .arg(shift_coulomb)
                                 .arg(taylor_power_expression("max_alpha", taylor_power))
                                 .toStdString();
        }
        else
        {
            // this uses the following potentials
            //            Zacharias and McCammon, J. Chem. Phys., 1994, and also,
            //            Michel et al., JCTC, 2007
            //
            //   V_{LJ}(r) = 4 epsilon [ ( sigma^12 / (delta*sigma + r^2)^6 ) -
            //                           ( sigma^6  / (delta*sigma + r^2)^3 ) ]
            //
            //   delta = shift_delta * alpha
            //
            //   V_{coul}(r) = (1-alpha)^n q_i q_j / 4 pi eps_0 (delta+r^2)^(1/2)
            //
            //   delta = shift_coulomb^2 * alpha
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
            clj_expression = QString("coul_nrg+lj_nrg;"
                                     "coul_nrg=138.9354558466661*q1*q2*(((%1)/sqrt((%2*max_alpha*max_alpha)+r_safe^2))-(max_kappa/r_safe));"
                                     "lj_nrg=two_sqrt_epsilon1*two_sqrt_epsilon2*sig6*(sig6-1);"
                                     "sig6=(sigma^6)/(((sigma*delta) + r_safe^2)^3);"
                                     "delta=%3*max_alpha;"
                                     "r_safe=max(r, 0.001);"
                                     "max_kappa=max(kappa1, kappa2);"
                                     "max_alpha=max(alpha1, alpha2);"
                                     "sigma=half_sigma1+half_sigma2;")
                                 .arg(coulomb_power_expression("max_alpha", coulomb_power))
                                 .arg(shift_coulomb)
                                 .arg(shift_delta.to(SireUnits::nanometer))
                                 .toStdString();
        }

        ghost_ghostff = new OpenMM::CustomNonbondedForce(clj_expression);
        ghost_ghostff->setName("GhostGhostNonbondedForce");
        ghost_nonghostff = new OpenMM::CustomNonbondedForce(clj_expression);
        ghost_nonghostff->setName("GhostNonGhostNonbondedForce");

        ghost_ghostff->addPerParticleParameter("q");
        ghost_ghostff->addPerParticleParameter("half_sigma");
        ghost_ghostff->addPerParticleParameter("two_sqrt_epsilon");
        ghost_ghostff->addPerParticleParameter("alpha");
        ghost_ghostff->addPerParticleParameter("kappa");

        ghost_nonghostff->addPerParticleParameter("q");
        ghost_nonghostff->addPerParticleParameter("half_sigma");
        ghost_nonghostff->addPerParticleParameter("two_sqrt_epsilon");
        ghost_nonghostff->addPerParticleParameter("alpha");
        ghost_nonghostff->addPerParticleParameter("kappa");

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

    // Stage 4 is complete. We now have all(*) of the forces we need to run
    // a perturbable simulation. (*) well, we will define the restraint
    // forces in a much later stage after the particles have been added.

    ///
    /// Stage 5 - Copy in the atomic forcefield parameters
    ///
    /// We now go through and copy the atomic forcefield parameters from
    /// the data in the temporary OpenMMMolecule objects into the
    /// OpenMM forces that we've just created.
    ///

    // start_index keeps track of the index of the first atom in each molecule
    int start_index = 0;

    // this is the list of atoms added, in atom order
    QList<SireMol::Selector<SireMol::Atom>> order_of_added_atoms;

    // get the 1-4 scaling factors from the first molecule
    const double coul_14_scl = openmm_mols_data[0].ffinfo.electrostatic14ScaleFactor();
    const double lj_14_scl = openmm_mols_data[0].ffinfo.vdw14ScaleFactor();

    // this maps from molecule index to the start_index for the first
    // particle in that molecule
    QVector<int> start_indexes(nmols);

    // the index to the perturbable or QM molecule for the specified molecule
    // (i.e. the 5th perturbable molecule is the 10th molecule in the System)
    QHash<int, int> idx_to_pert_idx;
    QHash<int, int> idx_to_qm_idx;

    // just a holder for all of the custom parameters for the
    // ghost forces (prevents us having to continually re-allocate it)
    std::vector<double> custom_params = {0.0, 0.0, 0.0, 0.0, 0.0};

    // the sets of particle indexes for the ghost atoms and non-ghost atoms
    QVector<int> ghost_atoms;
    QVector<int> non_ghost_atoms;

    // count the number of atoms and ghost atoms
    int n_atoms = 0;
    int n_ghost_atoms = 0;

    for (int i = 0; i < nmols; ++i)
    {
        const auto &mol = openmm_mols_data[i];
        n_atoms += mol.nAtoms();
        n_ghost_atoms += mol.nGhostAtoms();
    }

    // there's probably lots of optimisations we can make if the
    // number of ghost atoms is zero...
    ghost_atoms.reserve(n_ghost_atoms);
    non_ghost_atoms.reserve(n_atoms - n_ghost_atoms);

    // the set of all ghost atoms, with the value
    // indicating if this is a from-ghost (true) or
    // a to-ghost (false)
    QVector<int> from_ghost_idxs;
    QVector<int> to_ghost_idxs;
    from_ghost_idxs.reserve(n_ghost_atoms);
    to_ghost_idxs.reserve(n_ghost_atoms);

    // loop over every molecule and add them one by one
    for (int i = 0; i < nmols; ++i)
    {
        // save the index in the OpenMM system of the first
        // particle for the first atom in this molecule
        start_indexes[i] = start_index;
        const auto &mol = openmm_mols_data[i];

        order_of_added_atoms.append(mol.atoms);

        // double-check that the molecule has a compatible forcefield with
        // the other molecules in this system
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

        // is this a perturbable molecule (and we haven't disabled perturbations)?
        if (any_perturbable and mol.isPerturbable())
        {
            // this hash holds the start indicies for the various
            // parameters for this molecule (e.g. bond, angle, CLJ parameters)
            QHash<QString, qint32> start_indicies;

            // add a perturbable molecule, recording the start index
            // for each of the forcefields
            start_indicies.reserve(7);

            start_indicies.insert("clj", start_index);
            start_indicies.insert("ghost/ghost", start_index);
            start_indicies.insert("ghost/non-ghost", start_index);

            // the start index for this molecules first bond, angle or
            // torsion parameters will be however many of these
            // parameters exist already (parameters are added
            // contiguously for each molecule)
            start_indicies.insert("bond", bondff->getNumBonds());
            start_indicies.insert("angle", angff->getNumAngles());
            start_indicies.insert("torsion", dihff->getNumTorsions());

            // we can now record this as a perturbable molecule
            // in the lambda lever. The returned index is the
            // index of this perturbable molecule in the list
            // of perturbable molecules (e.g. the first perturbable
            // molecule we find has index 0)
            auto pert_idx = lambda_lever.addPerturbableMolecule(mol,
                                                                start_indicies,
                                                                map);

            // and we can record the map from the molecule index
            // to the perturbable molecule index
            idx_to_pert_idx.insert(i, pert_idx);
        }

        // Copy in all of the atom (particle) parameters. These
        // are the masses, charge and LJ parameters.
        // There is a different code path depending on whether
        // or not there are any perturbable molecules, and
        // this particular molecule is perturbable
        auto masses_data = mol.masses.constData();
        auto cljs_data = mol.cljs.constData();
        auto alphas_data = mol.alphas.constData();
        auto kappas_data = mol.kappas.constData();

        if (any_perturbable and mol.isPerturbable())
        {
            // This is a perturbable molecule and we're modelling perturbations
            for (int j = 0; j < mol.molinfo.nAtoms(); ++j)
            {
                const bool is_from_ghost = mol.from_ghost_idxs.contains(j);
                const bool is_to_ghost = mol.to_ghost_idxs.contains(j);

                // add the particle - the OpenMMMolecule has already
                // ensured that the largest of the reference or perturbed
                // masses is used
                const int atom_index = start_index + j;

                // NEED TO UPDATE FIXED ATOMS WITH FIELD ATOMS INDEXES!!!
                if (fixed_atoms.contains(atom_index))
                {
                    // this is a fixed (zero mass) atom
                    system.addParticle(0.0);
                }
                else
                {
                    // this is a mobile atom (if its mass is > 0)
                    system.addParticle(masses_data[j]);
                }

                // now the reference CLJ parameters
                const auto &clj = cljs_data[j];

                // make sure that charges are added here - if all are zero,
                // the NonbondedForce will not include support for charge!
                double charge = boost::get<0>(clj);

                if (charge == 0.0)
                {
                    charge = 1.0e-6;
                }

                // reduced_q
                custom_params[0] = charge;
                // half_sigma
                custom_params[1] = 0.5 * boost::get<1>(clj);
                // two_sqrt_epsilon
                custom_params[2] = 2.0 * std::sqrt(boost::get<2>(clj));
                // alpha
                custom_params[3] = alphas_data[j];
                // kappa
                custom_params[4] = kappas_data[j];

                // Add the particle to the ghost and nonghost forcefields
                ghost_ghostff->addParticle(custom_params);
                ghost_nonghostff->addParticle(custom_params);

                if (is_from_ghost or is_to_ghost)
                {
                    // this is a ghost atom! We need to record this
                    // fact and make sure that we don't calculate
                    // the LJ energy using the standard cljff
                    ghost_atoms.append(atom_index);

                    if (is_from_ghost)
                    {
                        from_ghost_idxs.append(atom_index);
                    }
                    else
                    {
                        to_ghost_idxs.append(atom_index);
                    }

                    // don't include the LJ energy as this will be
                    // calculated using the ghost forcefields
                    // (the ghost forcefields include a coulomb term
                    //  that subtracts from whatever was calculated here)
                    cljff->addParticle(charge, 0.0, 0.0);
                }
                else
                {
                    // this isn't a ghost atom. Record this fact and
                    // just add it to the standard cljff as normal
                    cljff->addParticle(charge, boost::get<1>(clj),
                                       boost::get<2>(clj));
                    non_ghost_atoms.append(atom_index);
                }
            }
        }
        else
        {
            // Code path if this isn't a perturbable molecule or
            // we don't want to model perturbations
            for (int j = 0; j < mol.molinfo.nAtoms(); ++j)
            {
                // Add the particle to the system
                const int atom_index = start_index + j;

                // NEED TO UPDATE FIXED ATOMS WITH FIELD ATOMS INDEXES!!!
                if (fixed_atoms.contains(atom_index))
                {
                    // this is a fixed (zero mass) atom
                    system.addParticle(0.0);
                }
                else
                {
                    // this is a mobile atom (if its mass is > 0)
                    system.addParticle(masses_data[j]);
                }

                // Get the particle CLJ parameters
                const auto &clj = cljs_data[j];

                // Add the particle to the standard CLJ forcefield
                cljff->addParticle(boost::get<0>(clj), boost::get<1>(clj), boost::get<2>(clj));

                // We need to add this molecule to the ghost and ghost
                // forcefields if there are any perturbable molecules
                if (any_perturbable)
                {
                    // reduced charge
                    custom_params[0] = boost::get<0>(clj);
                    // half_sigma
                    custom_params[1] = 0.5 * boost::get<1>(clj);
                    // two_sqrt_epsilon
                    custom_params[2] = 2.0 * std::sqrt(boost::get<2>(clj));
                    // alpha - is zero for non-ghost atoms
                    custom_params[3] = 0.0;
                    // kappa - is 0 for non-ghost atoms
                    custom_params[4] = 0.0;
                    ghost_ghostff->addParticle(custom_params);
                    ghost_nonghostff->addParticle(custom_params);
                    non_ghost_atoms.append(atom_index);
                }
            }
        }

        // now add all of the bond parameters
        for (const auto &bond : mol.bond_params)
        {
            bondff->addBond(boost::get<0>(bond) + start_index,
                            boost::get<1>(bond) + start_index,
                            boost::get<2>(bond), boost::get<3>(bond));
        }

        // now add all of the angles
        for (const auto &ang : mol.ang_params)
        {
            angff->addAngle(boost::get<0>(ang) + start_index,
                            boost::get<1>(ang) + start_index,
                            boost::get<2>(ang) + start_index,
                            boost::get<3>(ang), boost::get<4>(ang));
        }

        // now add all of the dihedrals and impropers
        for (const auto &dih : mol.dih_params)
        {
            dihff->addTorsion(boost::get<0>(dih) + start_index,
                              boost::get<1>(dih) + start_index,
                              boost::get<2>(dih) + start_index,
                              boost::get<3>(dih) + start_index,
                              boost::get<4>(dih), boost::get<5>(dih), boost::get<6>(dih));
        }

        for (const auto &constraint : mol.constraints)
        {
            const auto atom0 = boost::get<0>(constraint);
            const auto atom1 = boost::get<1>(constraint);

            const auto mass0 = system.getParticleMass(atom0 + start_index);
            const auto mass1 = system.getParticleMass(atom1 + start_index);

            if (mass0 != 0 and mass1 != 0)
            {

                system.addConstraint(atom0 + start_index,
                                     atom1 + start_index,
                                     boost::get<2>(constraint));
            }
            // else we will need to think about how to constrain bonds
            // involving fixed atoms. Could we fix the other atom too?
        }

        start_index += mol.masses.count();
    }

    /// Finally tell the ghost forcefields about the ghost and non-ghost
    /// interaction groups, so that they can correctly calculate the
    /// ghost/ghost and ghost/non-ghost energies
    if (ghost_ghostff != 0 and ghost_nonghostff != 0)
    {
        // set up the interaction groups - ghost / non-ghost
        //                                 ghost / ghost
        std::set<int> ghost_atoms_set(ghost_atoms.begin(), ghost_atoms.end());
        std::set<int> non_ghost_atoms_set(non_ghost_atoms.begin(), non_ghost_atoms.end());
        ghost_ghostff->addInteractionGroup(ghost_atoms_set, ghost_atoms_set);
        ghost_nonghostff->addInteractionGroup(ghost_atoms_set, non_ghost_atoms_set);
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

    /// Stage 5 is complete. We have added all of the parameter data
    /// for the molecules to the OpenMM forces

    ///
    /// Stage 6 - Set up the exceptions and perturbable constraints
    ///
    /// We now have to add all of the exceptions to the non-bonded
    /// forces (including the ghost forces). Exceptions are overrides
    /// that replace or switch off pair-pair energies between
    /// pairs of particles. This part of the code is the slowest
    /// of the entire function, as it involves lots of atom-atom
    /// pair loops and can create a large exception list
    ///
    /// We will also add all of the perturbable constraints here
    ///
    /// (we need to remember which ghost-ghost interactions we have
    ///  excluded, so that we don't double-exclude them later)
    QSet<IndexPair> excluded_ghost_pairs;
    excluded_ghost_pairs.reserve((n_ghost_atoms * n_ghost_atoms) / 2);

    for (int i = 0; i < nmols; ++i)
    {
        int start_index = start_indexes[i];
        const auto &mol = openmm_mols_data[i];

        QVector<boost::tuple<int, int>> exception_idxs;
        QVector<int> constraint_idxs;

        const bool is_perturbable = any_perturbable and mol.isPerturbable();

        if (is_perturbable)
        {
            exception_idxs = QVector<boost::tuple<int, int>>(mol.exception_params.count(),
                                                             boost::make_tuple(-1, -1));
            constraint_idxs = QVector<int>(mol.perturbable_constraints.count(), -1);

            // do all of the perturbable constraints
            for (int j = 0; j < mol.perturbable_constraints.count(); ++j)
            {
                const auto &constraint = mol.perturbable_constraints[j];

                const auto atom0 = boost::get<0>(constraint);
                const auto atom1 = boost::get<1>(constraint);

                const auto mass0 = system.getParticleMass(atom0 + start_index);
                const auto mass1 = system.getParticleMass(atom1 + start_index);

                if (mass0 != 0 and mass1 != 0)
                {
                    // add the constraint using the reference state length
                    auto idx = system.addConstraint(atom0 + start_index,
                                                    atom1 + start_index,
                                                    boost::get<2>(constraint));

                    constraint_idxs[j] = idx;
                }
            }
        }

        for (int j = 0; j < mol.exception_params.count(); ++j)
        {
            const auto &param = mol.exception_params[j];

            const auto atom0 = boost::get<0>(param);
            const auto atom1 = boost::get<1>(param);
            const auto coul_14_scale = boost::get<2>(param);
            const auto lj_14_scale = boost::get<3>(param);

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
                    idx = cljff->addException(boost::get<0>(p), boost::get<1>(p),
                                              boost::get<2>(p), 1e-9,
                                              1e-9, true);

                    // Whether this is a to/from ghost interaction.
                    auto to_from_ghost = (from_ghost_idxs.contains(atom0) and to_ghost_idxs.contains(atom1)) or
                                         (from_ghost_idxs.contains(atom1) and to_ghost_idxs.contains(atom0));

                    if (not to_from_ghost and (coul_14_scl != 0 or lj_14_scl != 0))
                    {
                        // this is a 1-4 interaction that should be added
                        // to the ghost-14 forcefield
                        if (ghost_14ff != 0)
                        {
                            // parameters are q, sigma, four_epsilon, alpha(0) and kappa(1)
                            std::vector<double> params14 =
                                {boost::get<2>(p), boost::get<3>(p),
                                 4.0 * boost::get<4>(p), 0.0, 1.0};

                            if (params14[0] == 0)
                                // cannot use zero params in case they are
                                // eagerly removed
                                params14[0] = 1e-9;

                            if (params14[1] == 0)
                                params14[1] = 1e-9;

                            nbidx = ghost_14ff->addBond(boost::get<0>(p),
                                                        boost::get<1>(p),
                                                        params14);
                        }
                    }

                    if (atom0_is_ghost and atom1_is_ghost)
                    {
                        // remember that this ghost-ghost interaction is already excluded
                        excluded_ghost_pairs.insert(IndexPair(boost::get<0>(p), boost::get<1>(p)));
                    }
                }
                else
                {
                    idx = cljff->addException(boost::get<0>(p), boost::get<1>(p),
                                              boost::get<2>(p), boost::get<3>(p),
                                              boost::get<4>(p), true);
                }

                // these are the indexes of the exception in the
                // non-bonded forcefields and also the ghost-14 forcefield
                exception_idxs[j] = boost::make_tuple(idx, nbidx);
            }
            else
            {
                cljff->addException(boost::get<0>(p), boost::get<1>(p),
                                    boost::get<2>(p), boost::get<3>(p),
                                    boost::get<4>(p), true);
            }

            // we need to make sure that the list of exclusions in
            // the NonbondedForce match those in the CustomNonbondedForces
            if (ghost_ghostff != 0)
            {
                ghost_ghostff->addExclusion(boost::get<0>(p), boost::get<1>(p));
                ghost_nonghostff->addExclusion(boost::get<0>(p), boost::get<1>(p));
            }
        }

        if (is_perturbable)
        {
            auto pert_idx = idx_to_pert_idx.value(i, openmm_mols.count() + 1);
            lambda_lever.setExceptionIndicies(pert_idx,
                                             "clj", exception_idxs);
            lambda_lever.setConstraintIndicies(pert_idx,
                                               constraint_idxs);
        }
    }

    // go through all of the ghost atoms and exclude interactions
    // between from_ghosts and to_ghosts
    for (const auto &from_ghost_idx : from_ghost_idxs)
    {
        // work out the molecule index for the from ghost atom
        int mol_from = 0;
        while (start_indexes[mol_from] <= from_ghost_idx and mol_from < nmols)
            mol_from++;

        for (const auto &to_ghost_idx : to_ghost_idxs)
        {
            // work out the molecule index for the to ghost atom
            int mol_to = 0;
            while (start_indexes[mol_to] <= to_ghost_idx and mol_to < nmols)
                mol_to++;

            if (not excluded_ghost_pairs.contains(IndexPair(from_ghost_idx, to_ghost_idx)))
            {
                // only exclude if we haven't already excluded this pair
                // and if the two atoms are in the same molecule
                if (mol_from == mol_to)
                {
                    ghost_ghostff->addExclusion(from_ghost_idx, to_ghost_idx);
                    ghost_nonghostff->addExclusion(from_ghost_idx, to_ghost_idx);
                    cljff->addException(from_ghost_idx, to_ghost_idx,
                                        0.0, 1e-9, 1e-9, true);
                }
            }
        }
    }

    // Stage 6 is complete. We have set up all of the exceptions. The
    // total energy / force calculated for the system should now be
    // correct.

    ///
    /// Stage 7 - Set up the restraints
    ///
    /// In this stage we set up all of the user-specified restraints.
    /// These give extra forces that are not defined in the list of
    /// molecules, and for which we need to know the indexes of
    /// all of the atoms (hence why we do this at the end)
    ///

    /// Some of the restraints will depend on fixed "anchor" atoms.
    /// This is the list of the coordinates of all of the anchor
    /// atoms (we will use this later when we get all coordinates
    /// and velocities)
    std::vector<OpenMM::Vec3> anchor_coords;

    if (map.specified("restraints"))
    {
        // All restraints are provided via the "restraints" key in the map.
        // This should either be a single Restraints object, or a
        // PropertyList - the easiest thing is to turn this into a
        // PropertyList first
        auto all_restraints = map["restraints"].value().asAnArray();

        // loop over all of the restraints groups and add them
        for (const auto &prop : all_restraints.toList())
        {
            if (not prop.read().isA<SireMM::Restraints>())
                throw SireError::invalid_cast(QObject::tr(
                                                  "Cannot convert an object of type %1 to a SireMM::Restraints")
                                                  .arg(prop.read().what()),
                                              CODELOC);

            // we now need to choose what to do based on the type of restraint...
            if (prop.read().isA<SireMM::DihedralRestraints>())
            {
                _add_dihedral_restraints(prop.read().asA<SireMM::DihedralRestraints>(),
                                         system, lambda_lever, start_index);
            }
            else if (prop.read().isA<SireMM::AngleRestraints>())
            {
                _add_angle_restraints(prop.read().asA<SireMM::AngleRestraints>(),
                                      system, lambda_lever, start_index);
            }
            else if (prop.read().isA<SireMM::PositionalRestraints>())
            {
                _add_positional_restraints(prop.read().asA<SireMM::PositionalRestraints>(),
                                           system, lambda_lever, anchor_coords, start_index);
            }
            else if (prop.read().isA<SireMM::BondRestraints>())
            {
                _add_bond_restraints(prop.read().asA<SireMM::BondRestraints>(),
                                     system, lambda_lever, start_index);
            }
            else if (prop.read().isA<SireMM::InverseBondRestraints>())
            {
                _add_inverse_bond_restraints(prop.read().asA<SireMM::InverseBondRestraints>(),
                                     system, lambda_lever, start_index);
            }
            else if (prop.read().isA<SireMM::BoreschRestraints>())
            {
                _add_boresch_restraints(prop.read().asA<SireMM::BoreschRestraints>(),
                                        system, lambda_lever, start_index);
            }
        }
    }

    ///
    /// Stage 8 - Copy across the coordinates and velocities
    ///
    /// In this final(!) stage we copy out all of the atoms (and anchors)
    /// coordinates and velocities so that they can be returned
    /// as metadata to be added to the OpenMM integrator
    ///

    // shared pointer to these coorindates and velocities
    std::shared_ptr<std::vector<OpenMM::Vec3>> coords, vels;

    const int natoms = start_index;
    const int nanchors = anchor_coords.size();

    // allocate memory
    coords.reset(new std::vector<OpenMM::Vec3>(natoms + nanchors));
    vels.reset(new std::vector<OpenMM::Vec3>(natoms + nanchors));

    auto coords_data = coords->data();
    auto vels_data = vels->data();

    const int *start_indexes_data = start_indexes.constData();

    // now copy the atomic data into the arrays
    if (SireBase::should_run_in_parallel(nmols, map))
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

    // now copy in the positions of the anchors
    if (nanchors > 0)
    {
        for (int i = 0; i < nanchors; ++i)
        {
            coords_data[natoms + i] = anchor_coords[i];
            vels_data[natoms + i] = OpenMM::Vec3(0, 0, 0);
        }
    }

    // All done - we can return the metadata
    return OpenMMMetaData(SireMol::SelectorM<SireMol::Atom>(order_of_added_atoms),
                          coords, vels, boxvecs, lambda_lever);
}
