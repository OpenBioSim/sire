#ifndef SIREOPENMM_OPENMMMOLECULE_H
#define SIREOPENMM_OPENMMMOLECULE_H

#include <OpenMM.h>

#include "SireMol/moleculeinfo.h"
#include "SireMol/core.h"
#include "SireMol/atom.h"
#include "SireMol/selector.hpp"

#include "SireMM/mmdetail.h"
#include "SireMM/excludedpairs.h"
#include "SireMM/amberparams.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{

    /** Internal class used to hold all of the extracted information
     *  of an OpenMM Molecule. You should not use this outside
     *  of the sire_to_openmm_system function. It holds lots of scratch
     *  data that may change in the future.
     */
    class OpenMMMolecule
    {
    public:
        enum CONSTRAIN_TYPE
        {
            CONSTRAIN_NONE = 0x00000000,
            CONSTRAIN_BONDS = 0x00000001,
            CONSTRAIN_HBONDS = 0x00000010,
            CONSTRAIN_HANGLES = 0x00001000,
            CONSTRAIN_NOT_PERTURBED = 0x00010000
        };

        OpenMMMolecule();
        OpenMMMolecule(const SireMol::Molecule &mol,
                       const SireBase::PropertyMap &map);

        OpenMMMolecule(const SireMol::Molecule &mol,
                       const SireBase::PropertyMap &map0,
                       const SireBase::PropertyMap &map1);

        ~OpenMMMolecule();

        bool operator==(const OpenMMMolecule &other) const;
        bool operator!=(const OpenMMMolecule &other) const;

        void copyInCoordsAndVelocities(OpenMM::Vec3 *coords,
                                       OpenMM::Vec3 *velocities) const;

        bool isPerturbable() const;

        QVector<double> getCharges() const;
        QVector<double> getSigmas() const;
        QVector<double> getEpsilons() const;
        QVector<double> getAlphas() const;
        QVector<double> getKappas() const;

        QVector<double> getBondKs() const;
        QVector<double> getBondLengths() const;

        QVector<double> getAngleKs() const;
        QVector<double> getAngleSizes() const;

        QVector<int> getTorsionPeriodicities() const;
        QVector<double> getTorsionPhases() const;
        QVector<double> getTorsionKs() const;

        QVector<std::pair<int, int>> getExceptionAtoms() const;

        QVector<double> getChargeScales() const;
        QVector<double> getLJScales() const;

        bool isGhostAtom(int atom) const;

        std::tuple<int, int, double, double, double>
        getException(int atom0, int atom1,
                     int start_index,
                     double coul_14_scl,
                     double lj_14_scl) const;

        /** All the member data is public as this is an internal
         *  class. This class should not be used outside of
         *  this SireOpenMM converter library.
         *
         *  Any values are in the internal units of OpenMM
         */
        SireMol::MolNum number;

        /** The molecule info that contains metadata about the molecule */
        SireMol::MoleculeInfo molinfo;

        /** All of the atoms, in the order they should appear in the
         *  OpenMM context
         */
        SireMol::Selector<SireMol::Atom> atoms;

        /** The forcefield info that contains metadata about the parameters */
        SireMM::MMDetail ffinfo;

        /** Coordinates */
        QVector<OpenMM::Vec3> coords;

        /** Velocities */
        QVector<OpenMM::Vec3> vels;

        /** Masses */
        QVector<double> masses;

        /** Indexes of light atoms */
        QSet<int> light_atoms;

        /** Indexes of virtual sites */
        QSet<int> virtual_sites;

        /** Charge and LJ parameters (sigma / epsilon) */
        QVector<std::tuple<double, double, double>> cljs;

        /** Set of 1-4 or excluded pairs
            (with coulomb and LJ scaling factors) */
        QVector<std::tuple<int, int, double, double>> exception_params;

        /** All the bond parameters */
        QVector<std::tuple<int, int, double, double>> bond_params;

        /** All the angle parameters */
        QVector<std::tuple<int, int, int, double, double>> ang_params;

        /** All the dihedral and improper parameters */
        QVector<std::tuple<int, int, int, int, int, double, double>> dih_params;

        /** All the constraints */
        QVector<std::tuple<int, int, double>> constraints;

        /** All of the perturbable constraints - these include the r0 values */
        QVector<std::tuple<int, int, double, double>> perturbable_constraints;

        /** The molecule perturbed molecule, if this is perturbable */
        std::shared_ptr<OpenMMMolecule> perturbed;

        /** The indicies of the added exceptions - only populated
         *  if this is a peturbable molecule */
        QHash<QString, QVector<std::pair<int, int>>> exception_idxs;

        /** The property map used to get the perturbable properties -
         *  this is only non-default if the molecule is perturbable
         */
        SireBase::PropertyMap perturtable_map;

        /** The atoms that are missing any internal parameters */
        QSet<qint32> unbonded_atoms;

        /** Alpha values for all of the atoms. This is equal to zero for
         *  non-ghost atoms, and one for ghost atoms
         */
        QVector<double> alphas;

        /** Kappa values for all of the atoms. This is equal to zero for
         *  non-ghost atoms, and one for ghost atoms
         */
        QVector<double> kappas;

        /** The indexes of atoms that become ghosts in the
         *  perturbed state
         */
        QSet<int> to_ghost_idxs;

        /** The indexes of atoms that are ghosts in the reference
         *  state and are real in the perturbed state
         */
        QSet<int> from_ghost_idxs;

        /** What type of constraint to use */
        qint32 constraint_type;

        /** What type of constraint to use when the bond/angle involves
         *  perturbable atoms */
        qint32 perturbable_constraint_type;

    private:
        void constructFromAmber(const SireMol::Molecule &mol,
                                const SireMM::AmberParams &params,
                                const SireMM::AmberParams &params1,
                                const SireBase::PropertyMap &map,
                                bool is_perturbable);

        void buildExceptions(const SireMol::Molecule &mol,
                             QSet<qint64> &constrained_pairs,
                             const SireBase::PropertyMap &map);

        void alignInternals(const SireBase::PropertyMap &map);
    };

    /** This class holds all of the information of an OpenMM molecule
     *  that can be perturbed using a LambdaSchedule. The data is held
     *  in easy-to-access arrays, with guarantees that the arrays are
     *  compatible and the data is aligned.
     */
    class PerturbableOpenMMMolecule
    {
    public:
        PerturbableOpenMMMolecule();

        PerturbableOpenMMMolecule(const OpenMMMolecule &mol);
        PerturbableOpenMMMolecule(const PerturbableOpenMMMolecule &other);

        ~PerturbableOpenMMMolecule();

        bool operator==(const PerturbableOpenMMMolecule &other) const;
        bool operator!=(const PerturbableOpenMMMolecule &other) const;

        QVector<double> getAlphas0() const;
        QVector<double> getAlphas1() const;

        QVector<double> getKappas0() const;
        QVector<double> getKappas1() const;

        QVector<double> getCharges0() const;
        QVector<double> getCharges1() const;

        QVector<double> getSigmas0() const;
        QVector<double> getSigmas1() const;
        QVector<double> getEpsilons0() const;
        QVector<double> getEpsilons1() const;

        QVector<double> getBondKs0() const;
        QVector<double> getBondKs1() const;
        QVector<double> getBondLengths0() const;
        QVector<double> getBondLengths1() const;

        QVector<double> getAngleKs0() const;
        QVector<double> getAngleKs1() const;
        QVector<double> getAngleSizes0() const;
        QVector<double> getAngleSizes1() const;

        QVector<double> getTorsionKs0() const;
        QVector<double> getTorsionKs1() const;
        QVector<int> getTorsionPeriodicities0() const;
        QVector<int> getTorsionPeriodicities1() const;
        QVector<double> getTorsionPhases0() const;
        QVector<double> getTorsionPhases1() const;

        QVector<double> getChargeScales0() const;
        QVector<double> getChargeScales1() const;
        QVector<double> getLJScales0() const;
        QVector<double> getLJScales1() const;

        QSet<int> getToGhostIdxs() const;
        QSet<int> getFromGhostIdxs() const;

        bool isGhostAtom(int atom) const;

        QVector<std::pair<int, int>> getExceptionAtoms() const;

        QVector<std::pair<int, int>> getExceptionIndicies(const QString &name) const;

        void setExceptionIndicies(const QString &name,
                                  const QVector<std::pair<int, int>> &exception_idxs);

        void setConstraintIndicies(const QVector<int> &constraint_idxs);

        QVector<int> getConstraintIndicies() const;

        std::tuple<QVector<int>, QVector<double>, QVector<double>> getPerturbableConstraints() const;

    private:
        /** The array of parameters for the two end states, aligned
         *  so that they can be morphed via the LambdaLever
         */
        QVector<double> alpha0, alpha1;
        QVector<double> kappa0, kappa1;
        QVector<double> chg0, chg1;
        QVector<double> sig0, sig1;
        QVector<double> eps0, eps1;
        QVector<double> bond_k0, bond_k1;
        QVector<double> bond_r0, bond_r1;
        QVector<double> ang_k0, ang_k1;
        QVector<double> ang_t0, ang_t1;
        QVector<double> tors_k0, tors_k1;
        QVector<int> tors_periodicity0, tors_periodicity1;
        QVector<double> tors_phase0, tors_phase1;
        QVector<double> charge_scl0, charge_scl1;
        QVector<double> lj_scl0, lj_scl1;

        /** The indexes of atoms that become ghosts in the
         *  perturbed state
         */
        QSet<int> to_ghost_idxs;

        /** The indexes of atoms that are ghosts in the reference
         *  state and are real in the perturbed state
         */
        QSet<int> from_ghost_idxs;

        /** The indicies of the atoms in the exceptions, in exception order */
        QVector<std::pair<int, int>> exception_atoms;

        /** The indicies of the added exceptions - only populated
         *  if this is a peturbable molecule */
        QHash<QString, QVector<std::pair<int, int>>> exception_idxs;

        /** All of the perturbable constraints - these include the r0 values
         *  for both end states
         */
        QVector<std::tuple<int, int, double, double>> perturbable_constraints;

        /** The indicies of the added constraints - this should be equal
         *  to the number of perturbable constraints in the molecule
         */
        QVector<int> constraint_idxs;
    };

}

SIRE_END_HEADER

#endif
