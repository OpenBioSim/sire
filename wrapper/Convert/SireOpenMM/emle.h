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

#ifndef SIREOPENMM_EMLE_H
#define SIREOPENMM_EMLE_H

#include "OpenMM.h"
#include "openmm/Force.h"
#ifdef SIRE_USE_CUSTOMCPPFORCE
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCPPForceImpl.h"
#endif

#include "boost/python.hpp"
#include <boost/tuple/tuple.hpp>

#include <QMap>
#include <QVector>

#include "sireglobal.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "qmmm.h"

namespace bp = boost::python;

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    // A callback wrapper class to allow use of electrostatic embedding of
    // machine learning potentials via emle-engine.
    class EMLECallback
    {
    public:
        //! Default constructor.
        EMLECallback();

        //! Constructor
        /*! \param py_object
                A Python object that contains the callback function.

            \param callback
                The name of a callback method that take the following arguments:
                    - numbers_qm: A list of atomic numbers for the atoms in the ML region.
                    - charges_mm: A list of the MM charges in mod electron charge.
                    - xyz_qm: A vector of positions for the atoms in the ML region in Angstrom.
                    - xyz_mm: A vector of positions for the atoms in the MM region in Angstrom.
         */
        EMLECallback(bp::object, QString callback="_sire_callback");

        //! Constructor
        /*! \param numbers_qm
                A vector of atomic numbers for the atoms in the ML region.

            \param charges_mm
                A vector of the charges on the MM atoms in mod electron charge.

            \param xyz_qm
                A vector of positions for the atoms in the ML region in Angstrom.

            \param xyz_mm
                A vector of positions for the atoms in the MM region in Angstrom.

            \returns
                A tuple containing:
                    - The energy in kJ/mol.
                    - A vector of forces for the QM atoms in kJ/mol/nm.
                    - A vector of forces for the MM atoms in kJ/mol/nm.
         */
        boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>> call(
                QVector<int> numbers_qm,
                QVector<double> charges_mm,
                QVector<QVector<double>> xyz_qm,
                QVector<QVector<double>> xyz_mm
        ) const;

    private:
        bp::object py_object;
        QString callback;
    };

    class EMLEEngine : public SireBase::ConcreteProperty<EMLEEngine, QMMMForce>
    {
    public:
        //! Default constructor.
        EMLEEngine();

        //! Constructor
        /*! \param py_object
                An EMLECalculator Python object.

            \param cutoff
                The ML cutoff distance.

            \param neighbour_list_frequency
                The frequency at which the neighbour list is updated. (Number of steps.)
                If zero, then no neighbour list is used.

            \param lambda
                The lambda weighting factor. This can be used to interpolate between
                potentials for end-state correction calculations.
         */
        EMLEEngine(
            bp::object,
            SireUnits::Dimension::Length cutoff=7.5*SireUnits::angstrom,
            int neighbour_list_frequency=20,
            double lambda=1.0
        );

        //! Copy constructor.
        EMLEEngine(const EMLEEngine &other);

        //! Assignment operator.
        EMLEEngine &operator=(const EMLEEngine &other);

        //! Set the callback object.
        /*! \param callback
                A Python object that contains the callback function.
         */
        void setCallback(EMLECallback callback);

        //! Get the callback object.
        /*! \returns
                A Python object that contains the callback function.
         */
        EMLECallback getCallback() const;

        //! Get the lambda weighting factor.
        /*! \returns
                The lambda weighting factor.
         */
        double getLambda() const;

        //! Set the lambda weighting factor.
        /*! \param lambda
                The lambda weighting factor.
         */
        void setLambda(double lambda);

        //! Get the QM cutoff distance.
        /*! \returns
                The QM cutoff distance.
         */
        SireUnits::Dimension::Length getCutoff() const;

        //! Set the QM cutoff distance.
        /*! \param cutoff
                The QM cutoff distance.
         */
        void setCutoff(SireUnits::Dimension::Length cutoff);

        //! Get the neighbour list frequency.
        /*! \returns
                The neighbour list frequency.
         */
        int getNeighbourListFrequency() const;

        //! Set the neighbour list frequency.
        /*! \param neighbour_list_frequency
                The neighbour list frequency.
         */
        void setNeighbourListFrequency(int neighbour_list_frequency);

        //! Get the indices of the atoms in the QM region.
        /*! \returns
                A vector of atom indices for the QM region.
         */
        QVector<int> getAtoms() const;

        //! Set the list of atom indices for the QM region.
        /*! \param atoms
                A vector of atom indices for the QM region.
         */
        void setAtoms(QVector<int> atoms);

        //! Get the link atoms associated with each QM atom.
        /*! \returns
                A tuple containing:

            mm1_to_qm
                A dictionary mapping link atom (MM1) indices to the QM atoms to
                which they are bonded.

            mm1_to_mm2
                A dictionary of link atoms indices (MM1) to a list of the MM
                atoms to which they are bonded (MM2).

            bond_scale_factors
                A dictionary of link atom indices (MM1) to a list of the bond
                length scale factors between the QM and MM1 atoms. The scale
                factors are the ratio of the equilibrium bond lengths for the
                QM-L (QM-link) atom and QM-MM1 atom, i.e. R0(QM-L) / R0(QM-MM1),
                taken from the MM force field parameters for the molecule.

        */
        boost::tuple<QMap<int, int>, QMap<int, QVector<int>>, QMap<int, double>> getLinkAtoms() const;

        //! Set the link atoms associated with each QM atom.
        /*! \param mm1_to_qm
                A dictionary mapping link atom (MM1) indices to the QM atoms to
                which they are bonded.

            \param mm1_to_mm2
                A dictionary of link atoms indices (MM1) to a list of the MM
                atoms to which they are bonded (MM2).

            \param bond_scale_factors
                A dictionary of link atom indices (MM1) to a list of the bond
                length scale factors between the QM and MM1 atoms. The scale
                factors are the ratio of the equilibrium bond lengths for the
                QM-L (QM-link) atom and QM-MM1 atom, i.e. R0(QM-L) / R0(QM-MM1),
                taken from the MM force field parameters for the molecule.

        */
        void setLinkAtoms(QMap<int, int> mm1_to_qm, QMap<int, QVector<int>> mm1_to_mm2, QMap<int, double> bond_scale_factors);

        //! Get the vector of MM2 atoms.
        /*! \returns
                A vector of MM2 atom indices.
         */
        QVector<int> getMM2Atoms() const;

        //! Get the atomic numbers for the atoms in the QM region.
        /*! \returns
                A vector of atomic numbers for the atoms in the QM region.
         */
        QVector<int> getNumbers() const;

        //! Set the atomic numbers for the atoms in the QM region.
        /*! \param numbers
                A vector of atomic numbers for the atoms in the QM region.
         */
        void setNumbers(QVector<int> numbers);

        //! Get the atomic charges of all atoms in the system.
        /*! \returns
                A vector of atomic charges for all atoms in the system.
         */
        QVector<double> getCharges() const;

        //! Set the atomic charges of all atoms in the system.
        /*! \param charges
                A vector of atomic charges for all atoms in the system.
         */
        void setCharges(QVector<double> charges);

        //! Return the C++ name for this class.
        static const char *typeName();

        //! Return the C++ name for this class.
        const char *what() const;

        //! Constructor
        /*! \param numbers_qm
                A vector of atomic numbers for the atoms in the ML region.

            \param charges_mm
                A vector of the charges on the MM atoms in mod electron charge.

            \param xyz_qm
                A vector of positions for the atoms in the ML region in Angstrom.

            \param xyz_mm
                A vector of positions for the atoms in the MM region in Angstrom.

            \returns
                A tuple containing:
                    - The energy in kJ/mol.
                    - A vector of forces for the QM atoms in kJ/mol/nm.
                    - A vector of forces for the MM atoms in kJ/mol/nm.
         */
        boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>> call(
                QVector<int> numbers_qm,
                QVector<double> charges_mm,
                QVector<QVector<double>> xyz_qm,
                QVector<QVector<double>> xyz_mm
        ) const;

    protected:
        OpenMM::ForceImpl *createImpl() const;

    private:
        EMLECallback callback;
        SireUnits::Dimension::Length cutoff;
        int neighbour_list_frequency;
        double lambda;
        QVector<int> atoms;
        QMap<int, int> mm1_to_qm;
        QMap<int, QVector<int>> mm1_to_mm2;
        QMap<int, double> bond_scale_factors;
        QVector<int> mm2_atoms;
        QVector<int> numbers;
        QVector<double> charges;
    };

#ifdef SIRE_USE_CUSTOMCPPFORCE
    class EMLEEngineImpl : public OpenMM::CustomCPPForceImpl
    {
    public:
        EMLEEngineImpl(const EMLEEngine &owner);

        ~EMLEEngineImpl();

        double computeForce(OpenMM::ContextImpl &context,
                            const std::vector<OpenMM::Vec3> &positions,
                            std::vector<OpenMM::Vec3> &forces);

        const EMLEEngine &getOwner() const;

    private:
        const EMLEEngine &owner;
        unsigned long long step_count=0;
        double cutoff;
        bool is_neighbour_list;
        int neighbour_list_frequency;
        double neighbour_list_cutoff;
        QSet<int> neighbour_list;
    };
#endif
}

Q_DECLARE_METATYPE(SireOpenMM::EMLEEngine)

SIRE_END_HEADER

#endif
