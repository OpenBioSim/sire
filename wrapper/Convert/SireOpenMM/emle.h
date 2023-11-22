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
                    - charges_mm: A list of the MM charges.
                    - xyz_qm: A vector of positions for the atoms in the ML region.
                    - xyz_mm: A vector of positions for the atoms in the MM region.
         */
        EMLECallback(bp::object, QString callback="_sire_callback");

        //! Constructor
        /*! \param numbers_qm
                A vector of atomic numbers for the atoms in the ML region.

            \param charges_mm
                A vector of the charges on the MM atoms.

            \param xyz_qm
                A vector of positions for the atoms in the ML region.

            \param xyz_mm
                A vector of positions for the atoms in the MM region.

            \returns
                A vector of forces for the QM and MM atoms.
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

    class EMLEEngine : public SireBase::ConcreteProperty<EMLEEngine, QMMMEngine>
    {
    public:
        //! Default constructor.
        EMLEEngine();

        //! Constructor
        /*! \param py_object
                An EMLECalculator Python object.

            \param cutoff
                The ML cutoff distance.

            \param lambda
                The lambda weighting factor. This can be used to interpolate between
                potentials for end-state correction calculations.
         */
        EMLEEngine(
            bp::object,
            SireUnits::Dimension::Length cutoff=8.0*SireUnits::angstrom,
            double lambda=1.0
        );

        //! Copy constructor.
        EMLEEngine(const EMLEEngine &other);

        //! Assignment operator.
        EMLEEngine &operator=(const EMLEEngine &other);

        //! Set the callback object.
        void setCallback(EMLECallback callback);

        //! Get the callback object.
        EMLECallback getCallback() const;

        //! Get the lambda weighting factor.
        double getLambda() const;

        //! Set the lambda weighting factor.
        void setLambda(double lambda);

        //! Get the QM cutoff distance.
        SireUnits::Dimension::Length getCutoff() const;

        //! Set the QM cutoff distance.
        void setCutoff(SireUnits::Dimension::Length cutoff);

        //! Get the indices of the atoms in the QM region.
        QVector<int> getAtoms() const;

        //! Set the list of atom indices for the QM region.
        void setAtoms(QVector<int> atoms);

        //! Return the C++ name for this class.
        static const char *typeName();

        //! Return the C++ name for this class.
        const char *what() const;

        //! Constructor
        /*! \param numbers_qm
                A vector of atomic numbers for the atoms in the ML region.

            \param charges_mm
                A vector of the charges on the MM atoms.

            \param xyz_qm
                A vector of positions for the atoms in the ML region.

            \param xyz_mm
                A vector of positions for the atoms in the MM region.

            \returns
                A vector of forces for the QM and MM atoms.
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
        double lambda;
        QVector<int> atoms;
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
    };
#endif
}

Q_DECLARE_METATYPE(SireOpenMM::EMLEEngine)

SIRE_END_HEADER

#endif
