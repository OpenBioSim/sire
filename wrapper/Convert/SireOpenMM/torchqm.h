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

#ifndef SIREOPENMM_TORCHQM_H
#define SIREOPENMM_TORCHQM_H

#include "OpenMM.h"
#include "openmm/Force.h"
#ifdef SIRE_USE_CUSTOMCPPFORCE
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCPPForceImpl.h"
#endif

#include "boost/python.hpp"
#include <boost/tuple/tuple.hpp>

#include <torch/script.h>

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
    class TorchQMForce;
}

QDataStream &operator<<(QDataStream &, const SireOpenMM::TorchQMForce &);
QDataStream &operator>>(QDataStream &, SireOpenMM::TorchQMForce &);

namespace SireOpenMM
{
    class TorchQMForce : public QMForce
    {
        friend QDataStream & ::operator<<(QDataStream &, const TorchQMForce &);
        friend QDataStream & ::operator>>(QDataStream &, TorchQMForce &);

    public:
        //! Default constructor.
        TorchQMForce();

        //! Constructor.
        /*  \param module_path
                The path to the serialised TorchScript module.

            \param torch_module
                The TorchScript module.

            \param cutoff
                The ML cutoff distance.

            \param neighbour_list_frequency
                The frequency at which the neighbour list is updated. (Number of steps.)
                If zero, then no neighbour list is used.

            \param is_mechanical
                A flag to indicate if mechanical embedding is being used.

            \param lambda
                The lambda weighting factor. This can be used to interpolate between
                potentials for end-state correction calculations.

            \param atoms
                A vector of atom indices for the QM region.

            \param mm1_to_qm
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

            \param mm2_atoms
                A vector of MM2 atom indices.

            \param numbers
                A vector of atomic charges for all atoms in the system.

            \param charges
                A vector of atomic charges for all atoms in the system.
         */
        TorchQMForce(
            QString module_path,
            SireUnits::Dimension::Length cutoff,
            int neighbour_list_frequency,
            bool is_mechanical,
            double lambda,
            QVector<int> atoms,
            QMap<int, int> mm1_to_qm,
            QMap<int, QVector<int>> mm1_to_mm2,
            QMap<int, double> bond_scale_factors,
            QVector<int> mm2_atoms,
            QVector<int> numbers,
            QVector<double> charges
        );

        //! Copy constructor.
        TorchQMForce(const TorchQMForce &other);

        //! Assignment operator.
        TorchQMForce &operator=(const TorchQMForce &other);

        //! Get the path to the serialised TorchScript module.
        /*! \returns
                The path to the serialised TorchScript module.
         */
        QString getModulePath() const;

        //! Set the path to the serialised TorchScript module.
        /*! \param module_path
                The path to the serialised TorchScript module.
         */
        void setModulePath(QString module_path);

        //! Get the TorchScript module.
        /*! \returns
                The TorchScript module.
         */
        torch::jit::script::Module getTorchModule() const;

        //! Get the lambda weighting factor.
        /*! \returns
                The lambda weighting factor.
         */
        double getLambda() const;

        //! Set the lambda weighting factor
        /* \param lambda
                The lambda weighting factor.
         */
        void setLambda(double lambda);

        //! Get the QM cutoff distance.
        /*! \returns
                The QM cutoff distance.
         */
        SireUnits::Dimension::Length getCutoff() const;

        //! Get the neighbour list frequency.
        /*! \returns
                The neighbour list frequency.
         */
        int getNeighbourListFrequency() const;

        //! Get the mechanical embedding flag.
        /*! \returns
                A flag to indicate if mechanical embedding is being used.
         */
        bool getIsMechanical() const;

        //! Get the indices of the atoms in the QM region.
        /*! \returns
                A vector of atom indices for the QM region.
         */
        QVector<int> getAtoms() const;

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

        //! Get the atomic charges of all atoms in the system.
        /*! \returns
                A vector of atomic charges for all atoms in the system.
         */
        QVector<double> getCharges() const;

        //! Return the C++ name for this class.
        static const char *typeName();

        //! Return the C++ name for this class.
        const char *what() const;

    protected:
        OpenMM::ForceImpl *createImpl() const;

    private:
        QString module_path;
        torch::jit::script::Module torch_module;
        SireUnits::Dimension::Length cutoff;
        int neighbour_list_frequency;
        bool is_mechanical;
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
    class TorchQMForceImpl : public OpenMM::CustomCPPForceImpl
    {
    public:
        TorchQMForceImpl(const TorchQMForce &owner);

        ~TorchQMForceImpl();

        double computeForce(OpenMM::ContextImpl &context,
                            const std::vector<OpenMM::Vec3> &positions,
                            std::vector<OpenMM::Vec3> &forces);

        const TorchQMForce &getOwner() const;

    private:
        const TorchQMForce &owner;
        torch::jit::script::Module torch_module;
        unsigned long long step_count=0;
        double cutoff;
        bool is_neighbour_list;
        int neighbour_list_frequency;
        double neighbour_list_cutoff;
        QSet<int> neighbour_list;
        int max_num_mm=0;
    };
#endif

    class TorchQMEngine : public SireBase::ConcreteProperty<TorchQMEngine, QMEngine>
    {
    public:
        //! Default constructor.
        TorchQMEngine();

        //! Constructor
        /*! \param module_path
                The path to the serialised TorchScript module.

            \param cutoff
                The ML cutoff distance.

            \param neighbour_list_frequency
                The frequency at which the neighbour list is updated. (Number of steps.)
                If zero, then no neighbour list is used.

            \param is_mechanical
                A flag to indicate if mechanical embedding is being used.

            \param lambda
                The lambda weighting factor. This can be used to interpolate between
                potentials for end-state correction calculations.
         */
        TorchQMEngine(
            QString module_path,
            SireUnits::Dimension::Length cutoff=7.5*SireUnits::angstrom,
            int neighbour_list_frequency=20,
            bool is_mechanical=false,
            double lambda=1.0
        );

        //! Copy constructor.
        TorchQMEngine(const TorchQMEngine &other);

        //! Assignment operator.
        TorchQMEngine &operator=(const TorchQMEngine &other);

        //! Get the path to the serialised TorchScript module.
        /*! \returns
                The path to the serialised TorchScript module.
         */
        QString getModulePath() const;

        //! Set the path to the serialised TorchScript module.
        /*! \param module_path
                The path to the serialised TorchScript module.
         */
        void setModulePath(QString module_path);

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

        //! Get the mechanical embedding flag.
        /*! \returns
                A flag to indicate if mechanical embedding is being used.
         */
        bool getIsMechanical() const;

        //! Set the mechanical embedding flag.
        /*! \param is_mechanical
                A flag to indicate if mechanical embedding is being used.
         */
        void setIsMechanical(bool is_mechanical);

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

        //! Create an EMLE force object.
        QMForce* createForce() const;

    private:
        QString module_path;
        SireUnits::Dimension::Length cutoff;
        int neighbour_list_frequency;
        bool is_mechanical;
        double lambda;
        QVector<int> atoms;
        QMap<int, int> mm1_to_qm;
        QMap<int, QVector<int>> mm1_to_mm2;
        QMap<int, double> bond_scale_factors;
        QVector<int> mm2_atoms;
        QVector<int> numbers;
        QVector<double> charges;
    };
}

Q_DECLARE_METATYPE(SireOpenMM::TorchQMForce)
Q_DECLARE_METATYPE(SireOpenMM::TorchQMEngine)

SIRE_EXPOSE_CLASS(SireOpenMM::TorchQMForce)
SIRE_EXPOSE_CLASS(SireOpenMM::TorchQMEngine)

SIRE_END_HEADER

#endif
