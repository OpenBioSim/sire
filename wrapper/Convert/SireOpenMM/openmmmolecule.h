#ifndef SIREOPENMM_OPENMMMOLECULE_H
#define SIREOPENMM_OPENMMMOLECULE_H

#include <OpenMM.h>

#include "SireMol/moleculeinfo.h"
#include "SireMol/core.h"

#include "SireMM/mmdetail.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{

    /** Internal class used to hold all of the extracted information
     *  of an OpenMM Molecule
     */
    class OpenMMMolecule
    {
    public:
        OpenMMMolecule();
        OpenMMMolecule(const SireMol::Molecule &mol,
                       const SireBase::PropertyMap &map);
        ~OpenMMMolecule();

        void copyInCoordsAndVelocities(OpenMM::Vec3 *coords,
                                       OpenMM::Vec3 *velocities) const;

        /** All the member data is public as this is an internal
         *  class. This class should not be used outside of
         *  this SireOpenMM converter library.
         *
         *  Any values are in the internal units of OpenMM
         */

        /** The molecule info that contains metadata about the molecule */
        SireMol::MoleculeInfo molinfo;

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

        /** Indexes of all bond pairs */
        QVector<std::pair<int, int>> bond_pairs;

        /** Indexes of pairs with custom 1-4 interactions */
        QVector<std::tuple<int, int, double, double, double>> custom_pairs;

        /** All the bond parameters */
        QVector<std::tuple<int, int, double, double>> bond_params;

        /** All the angle parameters */
        QVector<std::tuple<int, int, int, double, double>> ang_params;

        /** All the dihedral and improper parameters */
        QVector<std::tuple<int, int, int, int, double, double, double>> dih_params;

    private:
        void constructFromAmber(const SireMol::Molecule &mol,
                                const SireBase::PropertyMap &map);
    };

}

SIRE_END_HEADER

#endif
