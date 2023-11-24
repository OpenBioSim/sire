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

#include "SireMaths/vector.h"
#include "SireVol/triclinicbox.h"

#include "emle.h"

using namespace SireMaths;
using namespace SireOpenMM;
using namespace SireVol;

class GILLock
{
public:
    GILLock()  { state_ = PyGILState_Ensure(); }
    ~GILLock() { PyGILState_Release(state_);   }
private:
    PyGILState_STATE state_;
};

EMLECallback::EMLECallback()
{
}

EMLECallback::EMLECallback(bp::object py_object, QString callback) :
    py_object(py_object), callback(callback)
{
}

boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>
EMLECallback::call(
    QVector<int> numbers_qm,
    QVector<double> charges_mm,
    QVector<QVector<double>> xyz_qm,
    QVector<QVector<double>> xyz_mm) const
{
    // Acquire GIL before calling Python code.
    GILLock lock;

    return bp::call_method<boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>>(
        this->py_object.ptr(),
        this->callback.toStdString().c_str(),
        numbers_qm,
        charges_mm,
        xyz_qm,
        xyz_mm
    );
}

EMLEEngine::EMLEEngine()
{
}

EMLEEngine::EMLEEngine(bp::object py_object, SireUnits::Dimension::Length cutoff, double lambda) :
    callback(py_object, "_sire_callback"),
    cutoff(cutoff),
    lambda(lambda)
{
}

EMLEEngine::EMLEEngine(const EMLEEngine &other) :
    callback(other.callback),
    cutoff(other.cutoff),
    lambda(other.lambda),
    atoms(other.atoms),
    numbers(other.numbers),
    charges(other.charges)
{
}

EMLEEngine &EMLEEngine::operator=(const EMLEEngine &other)
{
    this->callback = other.callback;
    this->cutoff = other.cutoff;
    this->lambda = other.lambda;
    this->atoms = other.atoms;
    this->numbers = other.numbers;
    this->charges = other.charges;
    return *this;
}

void EMLEEngine::setCallback(EMLECallback callback)
{
    this->callback = callback;
}

EMLECallback EMLEEngine::getCallback() const
{
    return this->callback;
}

void EMLEEngine::setLambda(double lambda)
{
    this->lambda = lambda;
}

double EMLEEngine::getLambda() const
{
    return this->lambda;
}

void EMLEEngine::setCutoff(SireUnits::Dimension::Length cutoff)
{
    this->cutoff = cutoff;
}

SireUnits::Dimension::Length EMLEEngine::getCutoff() const
{
    return this->cutoff;
}

QVector<int> EMLEEngine::getAtoms() const
{
    return this->atoms;
}

void EMLEEngine::setAtoms(QVector<int> atoms)
{
    this->atoms = atoms;
}

QVector<int> EMLEEngine::getNumbers() const
{
    return this->numbers;
}

void EMLEEngine::setNumbers(QVector<int> numbers)
{
    this->numbers = numbers;
}

QVector<double> EMLEEngine::getCharges() const
{
    return this->charges;
}

void EMLEEngine::setCharges(QVector<double> charges)
{
    this->charges = charges;
}

const char *EMLEEngine::typeName()
{
    return QMetaType::typeName(qMetaTypeId<EMLEEngine>());
}

const char *EMLEEngine::what() const
{
    return EMLEEngine::typeName();
}

OpenMM::ForceImpl *EMLEEngine::createImpl() const
{
#ifdef SIRE_USE_CUSTOMCPPFORCE
    return new EMLEEngineImpl(*this);
#else
    throw SireError::unsupported(QObject::tr(
                                     "Unable to create a EMLEEngine because OpenMM::CustomCPPForceImpl "
                                     "is not available. You need to use OpenMM 8.1 or later."),
                                 CODELOC);
    return 0;
#endif
}

EMLEEngineImpl::EMLEEngineImpl(const EMLEEngine &owner) :
    OpenMM::CustomCPPForceImpl(owner),
    owner(owner)
{
}

EMLEEngineImpl::~EMLEEngineImpl()
{
}

const EMLEEngine &EMLEEngineImpl::getOwner() const
{
    return this->owner;
}

double EMLEEngineImpl::computeForce(
    OpenMM::ContextImpl &context,
    const std::vector<OpenMM::Vec3> &positions,
    std::vector<OpenMM::Vec3> &forces)
{
    // Get the current box vectors in nanometers.
    OpenMM::Vec3 box_x, box_y, box_z;
    context.getPeriodicBoxVectors(box_x, box_y, box_z);

    // Create a triclinic space, converting to Angstrom.
    TriclinicBox space(
        Vector(10*box_x[0], 10*box_x[1], 10*box_x[2]),
        Vector(10*box_y[0], 10*box_y[1], 10*box_y[2]),
        Vector(10*box_z[0], 10*box_z[1], 10*box_z[2])
    );

    // Store the QM atom indices.
    const auto qm_atoms = this->owner.getAtoms();

    // Initialise a vector to hold the current positions for the QM atoms.
    QVector<QVector<double>> xyz_qm(qm_atoms.size());
    QVector<Vector> xyz_qm_vec(qm_atoms.size());

    // Loop over all atoms in the QM region, get the OpenMM position, then
    // store it in Angstroms. We also store in Sire Vector format so that
    // we can use the Sire space to compute the distances. Positions are
    // stored using the minimum image convention with respect to the center
    // of geometry of the QM atoms.

    // Get the reference position. (First QM atom.)
    const auto ref_omm_pos = positions[this->owner.getAtoms()[0]];

    // Convert to Sire Vector format.
    Vector ref_vec(10*ref_omm_pos[0], 10*ref_omm_pos[1], 10*ref_omm_pos[2]);

    // Work out the center of geometry of the QM atoms.
    int i = 0;
    Vector center;
    for (const auto &idx : qm_atoms)
    {
        const auto &pos = positions[idx];
        Vector qm_vec(10*pos[0], 10*pos[1], 10*pos[2]);
        center += space.getMinimumImage(qm_vec, ref_vec);
        i++;
    }
    center /= i;

    // Loop over all of the QM atoms and store the positions with respect to
    // the reference position.
    i = 0;
    for (const auto &idx : qm_atoms)
    {
        const auto &pos = positions[idx];
        Vector qm_vec(10*pos[0], 10*pos[1], 10*pos[2]);
        const auto qm_vec_min = space.getMinimumImage(qm_vec, center);
        xyz_qm[i] = QVector<double>({qm_vec_min[0], qm_vec_min[1], qm_vec_min[2]});
        xyz_qm_vec[i] = qm_vec_min;
        i++;
    }

    // Next we need to work out the position of the MM atoms within the cutoff,
    // along with their charges. Here the cutoff is applied by including MM atoms
    // within the cutoff distance from any QM atom. Again we store the positions
    // using the minimum image with respect to the reference position.

    // Initialise a vector to hold the current positions for the MM atoms.
    QVector<QVector<double>> xyz_mm;

    // Initialise a vector to hold the charges for the MM atoms.
    QVector<double> charges_mm;

    // Initialise a list to hold the indices of the MM atoms.
    QVector<int> idx_mm;

    // Store the cutoff as a double in Angstom.
    const auto cutoff = this->owner.getCutoff().value();

    // Loop over all of the OpenMM positions.
    i = 0;
    for (const auto &pos : positions)
    {
        // Exclude QM atoms.
        if (not qm_atoms.contains(i))
        {
            // Whether to add the atom, i.e. it is an MM atom and is
            // within the cutoff.
            bool to_add = false;

            // Current MM atom position in Sire Vector format.
            const Vector mm_vec(10*pos[0], 10*pos[1], 10*pos[2]);

            // Loop over all of the QM atoms.
            for (const auto &qm_vec : xyz_qm_vec)
            {
                if (space.calcDist(mm_vec, qm_vec) < cutoff)
                {
                    // The current MM atom is within the cutoff, add it.
                    to_add = true;
                    break;
                }
            }

            // Store the MM atom information.
            if (to_add)
            {
                // Work out the minimum image position with respect to the
                // reference positions and add to the vector.
                const auto mm_vec_min = space.getMinimumImage(mm_vec, center);
                QVector<double> xyz = {mm_vec_min[0], mm_vec_min[1], mm_vec_min[2]};
                xyz_mm.append(xyz);

                // Add the charge and index.
                charges_mm.append(this->owner.getCharges()[i]);
                idx_mm.append(i);
            }
        }

        // Update the atom index.
        i++;
    }

    // Call the callback.
    auto result = this->owner.call(
        this->owner.getNumbers(),
        charges_mm,
        xyz_qm,
        xyz_mm
    );

    // Extract the results. This will automatically be returned in
    // OpenMM units.
    auto energy = result.get<0>();
    auto forces_qm = result.get<1>();
    auto forces_mm = result.get<2>();

    // Now update the force vector.

    // First the QM atoms.
    i = 0;
    for (const auto &force : forces_qm)
    {
        // Get the index of the atom.
        const auto idx = this->owner.getAtoms()[i];

        // Convert to OpenMM format.
        OpenMM::Vec3 omm_force(force[0], force[1], force[2]);

        // Update the force vector.
        forces[idx] = omm_force;

        // Update the atom index.
        i++;
    }

    // Now the MM atoms.
    i = 0;
    for (const auto &force : forces_mm)
    {
        // Get the index of the atom.
        const auto idx = idx_mm[i];

        // Convert to OpenMM format.
        OpenMM::Vec3 omm_force(force[0], force[1], force[2]);

        // Update the force vector.
        forces[idx] = omm_force;

        // Update the atom index.
        i++;
    }

    // Finally, return the energy.
    return energy;
}

boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>
EMLEEngine::call(
    QVector<int> numbers_qm,
    QVector<double> charges_mm,
    QVector<QVector<double>> xyz_qm,
    QVector<QVector<double>> xyz_mm) const
{
    return this->callback.call(numbers_qm, charges_mm, xyz_qm, xyz_mm);
}
