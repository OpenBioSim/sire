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
     *  of an OpenMM Molecule
     */
    class OpenMMMolecule
    {
    public:
        enum CONSTRAIN_TYPE
        {
            CONSTRAIN_NONE = 0x0000,
            CONSTRAIN_BONDS = 0x0001,
            CONSTRAIN_HBONDS = 0x0010,
            CONSTRAIN_ANGLES = 0x0100
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

        QVector<double> getBondKs() const;
        QVector<double> getBondLengths() const;

        QVector<double> getAngleKs() const;
        QVector<double> getAngleSizes() const;

        QVector<int> getTorsionPeriodicities() const;
        QVector<double> getTorsionPhases() const;
        QVector<double> getTorsionKs() const;

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
        QList<int> light_atoms;

        /** Indexes of virtual sites */
        QList<int> virtual_sites;

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

        /** The molecule perturbed molecule, if this is perturbable */
        std::shared_ptr<OpenMMMolecule> perturbed;

        /** The indicies of the added exceptions - only populated
         *  if this is a peturbable molecule */
        QHash<QString, QVector<std::pair<int, int>>> exception_idxs;

        /** The property map used to get the perturbable properties -
         *  this is only non-default if the molecule is perturbable
         */
        SireBase::PropertyMap perturtable_map;

        /** Alpha values for all of the atoms. This is equal to zero for
         *  non-ghost atoms, and one for ghost atoms
         */
        QVector<double> alphas;

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

}

SIRE_END_HEADER

#endif
