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

#include "SireError/errors.h"
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

EMLEEngine::EMLEEngine(
    bp::object py_object,
    SireUnits::Dimension::Length cutoff,
    int neighbour_list_frequency,
    double lambda) :
    callback(py_object, "_sire_callback"),
    cutoff(cutoff),
    neighbour_list_frequency(neighbour_list_frequency),
    lambda(lambda)
{
    if (this->neighbour_list_frequency < 0)
    {
        neighbour_list_frequency = 0;
    }
    if (this->lambda < 0.0)
    {
        this->lambda = 0.0;
    }
    else if (this->lambda > 1.0)
    {
        this->lambda = 1.0;
    }
}

EMLEEngine::EMLEEngine(const EMLEEngine &other) :
    callback(other.callback),
    cutoff(other.cutoff),
    neighbour_list_frequency(other.neighbour_list_frequency),
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
    this->neighbour_list_frequency = other.neighbour_list_frequency;
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
    // Clamp the lambda value.
    if (lambda < 0.0)
    {
        lambda = 0.0;
    }
    else if (lambda > 1.0)
    {
        lambda = 1.0;
    }
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

int EMLEEngine::getNeighbourListFrequency() const
{
    return this->neighbour_list_frequency;
}

void EMLEEngine::setNeighbourListFrequency(int neighbour_list_frequency)
{
    // Assume anything less than zero means no neighbour list.
    if (neighbour_list_frequency < 0)
    {
        neighbour_list_frequency = 0;
    }
    this->neighbour_list_frequency = neighbour_list_frequency;
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
#ifdef SIRE_USE_CUSTOMCPPFORCE
    // If this is the first step, then setup information for the neighbour list.
    if (this->step_count == 0)
    {
        // Store the cutoff as a double in Angstom.
        this->cutoff = this->owner.getCutoff().value();

        // The neighbour list cutoff is 20% larger than the cutoff.
        this->neighbour_list_cutoff = 1.2*this->cutoff;

        // Store the neighbour list update frequency.
        this->neighbour_list_frequency = this->owner.getNeighbourListFrequency();

        // Flag whether a neighbour list is used.
        this->is_neighbour_list = this->neighbour_list_frequency > 0;
    }

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

    // First loop over all QM atoms and store the positions.
    int i = 0;
    for (const auto &idx : qm_atoms)
    {
        const auto &pos = positions[idx];
        Vector qm_vec(10*pos[0], 10*pos[1], 10*pos[2]);
        xyz_qm_vec[i] = qm_vec;
        i++;
    }

    // Next sure that the QM atoms are whole (unwrapped).
    xyz_qm_vec = space.makeWhole(xyz_qm_vec);

    // Get the center of the QM atoms. We will use this as a reference when
    // re-imaging the MM atoms. Also store the QM atoms in the xyz_qm vector.
    Vector center;
    i = 0;
    for (const auto &qm_vec : xyz_qm_vec)
    {
        xyz_qm[i] = QVector<double>({qm_vec[0], qm_vec[1], qm_vec[2]});
        center += qm_vec;
        i++;
    }
    center /= i;

    // Initialise a vector to hold the current positions for the MM atoms.
    QVector<QVector<double>> xyz_mm;

    // Initialise a vector to hold the charges for the MM atoms.
    QVector<double> charges_mm;

    // Initialise a list to hold the indices of the MM atoms.
    QVector<int> idx_mm;

    // Manually work out the MM point charges and build the neigbour list.
    if (not this->is_neighbour_list or this->step_count % this->neighbour_list_frequency == 0)
    {
        // Clear the neighbour list.
        if (this->is_neighbour_list)
        {
            this->neighbour_list.clear();
        }

        i = 0;
        // Loop over all of the OpenMM positions.
        for (const auto &pos : positions)
        {
            // Exclude QM atoms.
            if (not qm_atoms.contains(i))
            {
                // Store the MM atom position in Sire Vector format.
                Vector mm_vec(10*pos[0], 10*pos[1], 10*pos[2]);

                // Loop over all of the QM atoms.
                for (const auto &qm_vec : xyz_qm_vec)
                {
                    // Work out the distance between the current MM atom and QM atoms.
                    const auto dist = space.calcDist(mm_vec, qm_vec);

                    // The current MM atom is within the neighbour list cutoff.
                    if (this->is_neighbour_list and dist < this->neighbour_list_cutoff)
                    {
                        // Insert the MM atom index into the neighbour list.
                        this->neighbour_list.insert(i);
                    }

                    // The current MM atom is within the cutoff, add it.
                    if (dist < cutoff)
                    {
                        // Work out the minimum image position with respect to the
                        // reference position and add to the vector.
                        mm_vec = space.getMinimumImage(mm_vec, center);
                        xyz_mm.append(QVector<double>({mm_vec[0], mm_vec[1], mm_vec[2]}));

                        // Add the charge and index.
                        charges_mm.append(this->owner.getCharges()[i]);
                        idx_mm.append(i);

                        // Exit the inner loop.
                        break;
                    }
                }
            }

            // Update the atom index.
            i++;
        }
    }
    // Use the neighbour list.
    else
    {
        // Loop over the MM atoms in the neighbour list.
        for (const auto &idx : this->neighbour_list)
        {
            // Store the MM atom position in Sire Vector format.
            Vector mm_vec(10*positions[idx][0], 10*positions[idx][1], 10*positions[idx][2]);

            // Loop over all of the QM atoms.
            for (const auto &qm_vec : xyz_qm_vec)
            {
                // The current MM atom is within the cutoff, add it.
                if (space.calcDist(mm_vec, qm_vec) < cutoff)
                {
                    // Work out the minimum image position with respect to the
                    // reference position and add to the vector.
                    mm_vec = space.getMinimumImage(mm_vec, center);
                    xyz_mm.append(QVector<double>({mm_vec[0], mm_vec[1], mm_vec[2]}));

                    // Add the charge and index.
                    charges_mm.append(this->owner.getCharges()[idx]);
                    idx_mm.append(idx);

                    // Exit the inner loop.
                    break;
                }
            }

        }
    }

    // Call the callback.
    auto result = this->owner.call(
        this->owner.getNumbers(),
        charges_mm,
        xyz_qm,
        xyz_mm
    );

    // Extract the results. These will automatically be returned in OpenMM units.
    auto energy = result.get<0>();
    auto forces_qm = result.get<1>();
    auto forces_mm = result.get<2>();

    // Store the current lambda weighting factor.
    const auto lambda = this->owner.getLambda();

    // Now update the force vector.

    // First the QM atoms.
    i = 0;
    for (const auto &force : forces_qm)
    {
        // Get the index of the atom.
        const auto idx = qm_atoms[i];

        // Convert to OpenMM format.
        OpenMM::Vec3 omm_force(force[0], force[1], force[2]);

        // Update the force vector.
        forces[idx] = lambda * omm_force;

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
        forces[idx] = lambda * omm_force;

        // Update the atom index.
        i++;
    }

    // Update the step count.
    this->step_count++;

    // Finally, return the energy.
    return lambda * energy;
#endif
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
