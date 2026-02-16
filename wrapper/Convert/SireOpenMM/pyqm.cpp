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

#include <mutex>

#include <QHash>
#include <QUuid>

#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"

#include "SireError/errors.h"
#include "SireMaths/vector.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireVol/triclinicbox.h"

#include "pyqm.h"

using namespace SireMaths;
using namespace SireOpenMM;
using namespace SireStream;
using namespace SireVol;

// The delta used to place virtual point charges either side of the MM2
// atoms, in nanometers.
static const double VIRTUAL_PC_DELTA = 0.01;

class GILLock
{
public:
    GILLock()  { state_ = PyGILState_Ensure(); }
    ~GILLock() { PyGILState_Release(state_);   }
private:
    PyGILState_STATE state_;
};

/////////
///////// Implementation of PyQMCallback
/////////

// A registry to store Python callback objects.
QHash<QUuid, bp::object> py_object_registry;

// A mutex to protect the registry.
std::mutex py_object_mutex;

// Set a callback Python object in the registry using a mutex.
void setPyObject(bp::object cb, QString uuid)
{
    std::lock_guard<std::mutex> lock(py_object_mutex);
    py_object_registry[uuid] = cb;
};

// Get a callback object from the registry using a mutex.
bp::object getPythonObject(QString uuid)
{
    std::lock_guard<std::mutex> lock(py_object_mutex);

    if (not py_object_registry.contains(uuid))
    {
        throw SireError::invalid_key(QObject::tr(
            "Unable to find UUID %1 in the PyQMForce callback registry.").arg(uuid),
            CODELOC);
    }

    return py_object_registry[uuid];
}

static const RegisterMetaType<PyQMCallback> r_pyqmcallback(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const PyQMCallback &pyqmcallback)
{
    writeHeader(ds, r_pyqmcallback, 1);

    SharedDataStream sds(ds);

    // Generate a unique identifier for the callback.
    auto uuid = QUuid::createUuid().toString();

    sds << uuid << pyqmcallback.name << pyqmcallback.is_method;

    // Set the Python object in the registry.
    setPyObject(pyqmcallback.py_object, uuid);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PyQMCallback &pyqmcallback)
{
    VersionID v = readHeader(ds, r_pyqmcallback);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        QString uuid;

        // Get the UUID of the Python object and the callback name.
        sds >> uuid >> pyqmcallback.name >> pyqmcallback.is_method;

        // Set the Python object.
        pyqmcallback.py_object = getPythonObject(uuid);
    }
    else
        throw version_error(v, "1", r_pyqmcallback, CODELOC);

    return ds;
}

PyQMCallback::PyQMCallback()
{
}

PyQMCallback::PyQMCallback(bp::object py_object, QString name) :
    py_object(py_object), name(name)
{
    // Is this a method or free function.
    if (name.isEmpty())
    {
        this->is_method = false;
    }
}

boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>
PyQMCallback::call(
    QVector<int> numbers_qm,
    QVector<double> charges_mm,
    QVector<QVector<double>> xyz_qm,
    QVector<QVector<double>> xyz_mm,
    QVector<QVector<double>> cell,
    QVector<int> idx_mm) const
{

    // Acquire GIL before calling Python code.
    GILLock lock;

    if (this->is_method)
    {
        try
        {
            return bp::call_method<boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>>(
                this->py_object.ptr(),
                this->name.toStdString().c_str(),
                numbers_qm,
                charges_mm,
                xyz_qm,
                xyz_mm,
                cell,
                idx_mm
            );
        }
        catch (const bp::error_already_set &)
        {
            PyErr_Print();
            throw SireError::process_error(QObject::tr(
                "An error occurred when calling the QM Python callback method"),
                CODELOC);
        }
    }
    else
    {
        try
        {
            return bp::call<boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>>(
                this->py_object.ptr(),
                numbers_qm,
                charges_mm,
                xyz_qm,
                xyz_mm,
                cell,
                idx_mm
            );
        }
        catch (const bp::error_already_set &)
        {
            PyErr_Print();
            throw SireError::process_error(QObject::tr(
                "An error occurred when calling the QM Python callback method"),
                CODELOC);
        }
    }
}

const char *PyQMCallback::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PyQMCallback>());
}

const char *PyQMCallback::what() const
{
    return PyQMCallback::typeName();
}

/////////
///////// Implementation of PyQMForce
/////////

static const RegisterMetaType<PyQMForce> r_pyqmforce(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const PyQMForce &pyqmforce)
{
    writeHeader(ds, r_pyqmforce, 1);

    SharedDataStream sds(ds);

    sds << pyqmforce.callback << pyqmforce.cutoff << pyqmforce.neighbour_list_frequency
        << pyqmforce.is_mechanical <<  pyqmforce.lambda << pyqmforce.atoms
        << pyqmforce.mm1_to_qm << pyqmforce.mm1_to_mm2 << pyqmforce.bond_scale_factors
        << pyqmforce.mm2_atoms << pyqmforce.numbers << pyqmforce.charges;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PyQMForce &pyqmforce)
{
    VersionID v = readHeader(ds, r_pyqmforce);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> pyqmforce.callback >> pyqmforce.cutoff >> pyqmforce.neighbour_list_frequency
            >> pyqmforce.is_mechanical >> pyqmforce.lambda >> pyqmforce.atoms
            >> pyqmforce.mm1_to_qm >> pyqmforce.mm1_to_mm2 >> pyqmforce.bond_scale_factors
            >> pyqmforce.mm2_atoms >> pyqmforce.numbers >> pyqmforce.charges;
    }
    else
        throw version_error(v, "1", r_pyqmforce, CODELOC);

    return ds;
}

PyQMForce::PyQMForce()
{
}

PyQMForce::PyQMForce(
    PyQMCallback callback,
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
    QVector<double> charges) :
    callback(callback),
    cutoff(cutoff),
    neighbour_list_frequency(neighbour_list_frequency),
    is_mechanical(is_mechanical),
    lambda(lambda),
    atoms(atoms),
    mm1_to_qm(mm1_to_qm),
    mm1_to_mm2(mm1_to_mm2),
    bond_scale_factors(bond_scale_factors),
    mm2_atoms(mm2_atoms),
    numbers(numbers),
    charges(charges)
{
}

PyQMForce::PyQMForce(const PyQMForce &other) :
    callback(other.callback),
    cutoff(other.cutoff),
    neighbour_list_frequency(other.neighbour_list_frequency),
    is_mechanical(other.is_mechanical),
    lambda(other.lambda),
    atoms(other.atoms),
    mm1_to_qm(other.mm1_to_qm),
    mm1_to_mm2(other.mm1_to_mm2),
    mm2_atoms(other.mm2_atoms),
    bond_scale_factors(other.bond_scale_factors),
    numbers(other.numbers),
    charges(other.charges)
{
}

PyQMForce &PyQMForce::operator=(const PyQMForce &other)
{
    this->callback = other.callback;
    this->cutoff = other.cutoff;
    this->neighbour_list_frequency = other.neighbour_list_frequency;
    this->is_mechanical = other.is_mechanical;
    this->lambda = other.lambda;
    this->atoms = other.atoms;
    this->mm1_to_qm = other.mm1_to_qm;
    this->mm1_to_mm2 = other.mm1_to_mm2;
    this->mm2_atoms = other.mm2_atoms;
    this->bond_scale_factors = other.bond_scale_factors;
    this->numbers = other.numbers;
    this->charges = other.charges;
    return *this;
}

void PyQMForce::setCallback(PyQMCallback callback)
{
    this->callback = callback;
}

PyQMCallback PyQMForce::getCallback() const
{
    return this->callback;
}

void PyQMForce::setLambda(double lambda)
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

double PyQMForce::getLambda() const
{
    return this->lambda;
}

SireUnits::Dimension::Length PyQMForce::getCutoff() const
{
    return this->cutoff;
}

int PyQMForce::getNeighbourListFrequency() const
{
    return this->neighbour_list_frequency;
}

bool PyQMForce::getIsMechanical() const
{
    return this->is_mechanical;;
}

QVector<int> PyQMForce::getAtoms() const
{
    return this->atoms;
}

boost::tuple<QMap<int, int>, QMap<int, QVector<int>>, QMap<int, double>> PyQMForce::getLinkAtoms() const
{
    return boost::make_tuple(this->mm1_to_qm, this->mm1_to_mm2, this->bond_scale_factors);
}

QVector<int> PyQMForce::getMM2Atoms() const
{
    return this->mm2_atoms;
}

QVector<int> PyQMForce::getNumbers() const
{
    return this->numbers;
}

QVector<double> PyQMForce::getCharges() const
{
    return this->charges;
}

const char *PyQMForce::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PyQMForce>());
}

const char *PyQMForce::what() const
{
    return PyQMForce::typeName();
}

boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>
PyQMForce::call(
    QVector<int> numbers_qm,
    QVector<double> charges_mm,
    QVector<QVector<double>> xyz_qm,
    QVector<QVector<double>> xyz_mm,
    QVector<QVector<double>> cell,
    QVector<int> idx_mm) const
{
    return this->callback.call(numbers_qm, charges_mm, xyz_qm, xyz_mm, cell, idx_mm);
}

/////////
///////// OpenMM Serialization
/////////

namespace OpenMM
{
    class PyQMForceProxy : public SerializationProxy {
        public:
            PyQMForceProxy() : SerializationProxy("PyQMForce")
            {
            };

            void serialize(const void* object, SerializationNode& node) const
            {
                // Serialize the object.
                QByteArray data;
                QDataStream ds(&data, QIODevice::WriteOnly);
                PyQMForce pyqmforce = *static_cast<const PyQMForce*>(object);
                ds << pyqmforce;

                // Set the version.
                node.setIntProperty("version", 0);

                // Set the note attribute.
                node.setStringProperty("note",
                 "This force only supports partial serialization, so can only be used "
                 "within the same session and memory space.");

                // Set the data by converting the QByteArray to a hexidecimal string.
                node.setStringProperty("data", data.toHex().data());
            };

            void* deserialize(const SerializationNode& node) const
            {
                // Check the version.
                int version = node.getIntProperty("version");
                if (version != 0)
                {
                    throw OpenMM::OpenMMException("Unsupported version number");
                }

                // Get the data as a std::string.
                auto string = node.getStringProperty("data");

                // Convert to hexidecimal.
                auto hex = QByteArray::fromRawData(string.data(), string.size());

                // Convert to a QByteArray.
                auto data = QByteArray::fromHex(hex);

                // Deserialize the object.
                QDataStream ds(data);
                PyQMForce pyqmforce;

                try
                {
                    ds >> pyqmforce;
                }
                catch (...)
                {
                    throw OpenMM::OpenMMException("Unable to find UUID in the PyQMForce callback registry.");
                }

                return new PyQMForce(pyqmforce);
            };
    };

    // Register the PyQMForce serialization proxy.
    extern "C" void registerPyQMSerializationProxies() {
        SerializationProxy::registerProxy(typeid(PyQMForce), new PyQMForceProxy());
    }
};

/////////
///////// Implementation of PyQMForceImpl
/////////

OpenMM::ForceImpl *PyQMForce::createImpl() const
{
#ifdef SIRE_USE_CUSTOMCPPFORCE
    return new PyQMForceImpl(*this);
#else
    throw SireError::unsupported(QObject::tr(
                                     "Unable to create an PyQMForceImpl because OpenMM::CustomCPPForceImpl "
                                     "is not available. You need to use OpenMM 8.1 or later."),
                                 CODELOC);
    return 0;
#endif
}

#ifdef SIRE_USE_CUSTOMCPPFORCE
PyQMForceImpl::PyQMForceImpl(const PyQMForce &owner) :
    OpenMM::CustomCPPForceImpl(owner),
    owner(owner)
{
}

PyQMForceImpl::~PyQMForceImpl()
{
}

const PyQMForce &PyQMForceImpl::getOwner() const
{
    return this->owner;
}

double PyQMForceImpl::computeForce(
    OpenMM::ContextImpl &context,
    const std::vector<OpenMM::Vec3> &positions,
    std::vector<OpenMM::Vec3> &forces)
{
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

    // Store the cell vectors in Angstrom.
    QVector<QVector<double>> cell = {
        {10*box_x[0], 10*box_x[1], 10*box_x[2]},
        {10*box_y[0], 10*box_y[1], 10*box_y[2]},
        {10*box_z[0], 10*box_z[1], 10*box_z[2]}
    };

    // Store the QM atomic indices and numbers.
    auto qm_atoms = this->owner.getAtoms();
    auto numbers = this->owner.getNumbers();

    // Store the link atom info. Link atoms are handled using the Charge Shift
    // method. See: https://www.ks.uiuc.edu/Research/qmmm
    const auto link_atoms = this->owner.getLinkAtoms();
    const auto mm1_to_qm = link_atoms.get<0>();
    const auto mm1_to_mm2 = link_atoms.get<1>();
    const auto bond_scale_factors = link_atoms.get<2>();
    const auto mm2_atoms = this->owner.getMM2Atoms();

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

    // Initialise a vector to hold the current positions and charges for the MM atoms.
    QVector<QVector<double>> xyz_mm;
    QVector<double> charges_mm;

    // Initialise a list to hold the indices of the MM atoms.
    QVector<int> idx_mm;

    // Store the current number of MM atoms.
    unsigned int num_mm = 0;

    // If we are using electrostatic embedding, the work out the MM point charges and
    // build the neighbour list.
    if (not this->owner.getIsMechanical())
    {
        // Initialise a vector to hold the current positions and charges for the virtual
        // point charges.
        QVector<QVector<double>> xyz_virtual;
        QVector<double> charges_virtual;

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
                // Exclude QM atoms or link atoms, which are handled later.
                if (not qm_atoms.contains(i) and
                    not mm1_to_mm2.contains(i) and
                    not mm2_atoms.contains(i))
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

        // Handle link atoms via the Charge Shift method.
        // See: https://www.ks.uiuc.edu/Research/qmmm
        for (const auto &idx: mm1_to_mm2.keys())
        {
            // Get the QM atom to which the current MM atom is bonded.
            const auto qm_idx = mm1_to_qm[idx];

            // Store the MM1 position in Sire Vector format, along with the
            // position of the QM atom to which it is bonded.
            Vector mm1_vec(10*positions[idx][0], 10*positions[idx][1], 10*positions[idx][2]);
            Vector qm_vec(10*positions[qm_idx][0], 10*positions[qm_idx][1], 10*positions[qm_idx][2]);

            // Work out the minimum image positions with respect to the reference position.
            mm1_vec = space.getMinimumImage(mm1_vec, center);
            qm_vec = space.getMinimumImage(qm_vec, center);

            // Work out the position of the link atom. Here we use a bond length
            // scale factor taken from the MM bond potential, i.e. R0(QM-L) / R0(QM-MM1),
            // where R0(QM-L) is the equilibrium bond length for the QM and link (L)
            // elements, and R0(QM-MM1) is the equilibrium bond length for the QM
            // and MM1 elements.
            const auto link_vec = qm_vec + bond_scale_factors[idx]*(mm1_vec - qm_vec);

            // Add to the QM positions.
            xyz_qm.append(QVector<double>({link_vec[0], link_vec[1], link_vec[2]}));

            // Add the MM1 index to the QM atoms vector.
            qm_atoms.append(qm_idx);

            // Append a hydrogen element to the numbers vector.
            numbers.append(1);

            // Store the number of MM2 atoms.
            const auto num_mm2 = mm1_to_mm2[idx].size();

            // Store the fractional charge contribution to the MM2 atoms and
            // virtual point charges.
            const auto frac_charge = this->owner.getCharges()[idx] / num_mm2;

            // Loop over the MM2 atoms and perform charge shifting. Here the MM1
            // charge is redistributed over the MM2 atoms and two virtual point
            // charges are added either side of the MM2 atoms in order to preserve
            // the MM1-MM2 dipole.
            for (const auto& mm2_idx : mm1_to_mm2[idx])
            {
                // Store the MM2 position in Sire Vector format.
                Vector mm2_vec(10*positions[mm2_idx][0], 10*positions[mm2_idx][1], 10*positions[mm2_idx][2]);

                // Work out the minimum image position with respect to the reference position.
                mm2_vec = space.getMinimumImage(mm2_vec, center);

                // Add to the MM positions.
                xyz_mm.append(QVector<double>({mm2_vec[0], mm2_vec[1], mm2_vec[2]}));

                // Add the charge and index.
                charges_mm.append(this->owner.getCharges()[mm2_idx] + frac_charge);
                idx_mm.append(mm2_idx);

                // Now add the virtual point charges.

                // Compute the normal vector from the MM1 to MM2 atom.
                const auto normal = (mm2_vec - mm1_vec).normalise();

                // Positive direction. (Away from MM1 atom.)
                auto xyz = mm2_vec + VIRTUAL_PC_DELTA*normal;
                xyz_virtual.append(QVector<double>({xyz[0], xyz[1], xyz[2]}));
                charges_virtual.append(-frac_charge);

                // Negative direction (Towards MM1 atom.)
                xyz = mm2_vec - VIRTUAL_PC_DELTA*normal;
                xyz_virtual.append(QVector<double>({xyz[0], xyz[1], xyz[2]}));
                charges_virtual.append(frac_charge);
            }
        }

        // Store the current number of MM atoms.
        num_mm = xyz_mm.size();

        // If there are any virtual point charges, then add to the MM positions
        // and charges.
        if (xyz_virtual.size() > 0)
        {
            xyz_mm.append(xyz_virtual);
            charges_mm.append(charges_virtual);
        }
    }

    // Call the callback.
    auto result = this->owner.call(
        numbers,
        charges_mm,
        xyz_qm,
        xyz_mm,
        cell,
        idx_mm
    );

    // Extract the results. These will automatically be returned in OpenMM units.
    auto energy = result.get<0>();
    auto forces_qm = result.get<1>();
    auto forces_mm = result.get<2>();

    // The current interpolation (weighting) parameter.
    double lambda;

    // Try to get the "lambda_emle" global parameter from the context.
    try
    {
        lambda = context.getParameter("lambda_emle");
    }
    catch (...)
    {
        // Try to get the "lambda_interpolate" global parameter from the context.
        try
        {
            lambda = context.getParameter("lambda_interpolate");
        }
        // Fall back on the lambda value stored in the PyQMForce object.
        catch (...)
        {
            lambda = this->owner.getLambda();
        }
    }

    // Clamp the lambda value.
    if (lambda < 0.0)
    {
        lambda = 0.0;
    }
    else if (lambda > 1.0)
    {
        lambda = 1.0;
    }

    // Now update the force vector.

    // First the QM atoms.
    for (int i=0; i<qm_atoms.size(); i++)
    {
        // Get the index of the atom.
        const auto idx = qm_atoms[i];

        // Convert to OpenMM format.
        OpenMM::Vec3 omm_force(forces_qm[i][0], forces_qm[i][1], forces_qm[i][2]);

        // Update the force vector.
        forces[idx] = lambda * omm_force;
    }

    // Now the MM atoms.
    for (int i=0; i<num_mm; i++)
    {
        // Get the index of the atom.
        const auto idx = idx_mm[i];

        // Convert to OpenMM format.
        OpenMM::Vec3 omm_force(forces_mm[i][0], forces_mm[i][1], forces_mm[i][2]);

        // Update the force vector.
        forces[idx] = lambda * omm_force;
    }

    // Update the step count.
    this->step_count++;

    // Finally, return the energy.
    return lambda * energy;
}
#endif

/////////
///////// Implementation of PyQMEngine
/////////

PyQMEngine::PyQMEngine() : ConcreteProperty<PyQMEngine, QMEngine>()
{
    // Register the serialization proxies.
    OpenMM::registerPyQMSerializationProxies();
}

PyQMEngine::PyQMEngine(
    bp::object py_object,
    QString name,
    SireUnits::Dimension::Length cutoff,
    int neighbour_list_frequency,
    bool is_mechanical,
    double lambda) :
    ConcreteProperty<PyQMEngine, QMEngine>(),
    callback(py_object, name),
    cutoff(cutoff),
    neighbour_list_frequency(neighbour_list_frequency),
    is_mechanical(is_mechanical),
    lambda(lambda)
{
    // Register the serialization proxies.
    OpenMM::registerPyQMSerializationProxies();

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

PyQMEngine::PyQMEngine(const PyQMEngine &other) :
    callback(other.callback),
    cutoff(other.cutoff),
    neighbour_list_frequency(other.neighbour_list_frequency),
    is_mechanical(other.is_mechanical),
    lambda(other.lambda),
    atoms(other.atoms),
    mm1_to_qm(other.mm1_to_qm),
    mm1_to_mm2(other.mm1_to_mm2),
    mm2_atoms(other.mm2_atoms),
    bond_scale_factors(other.bond_scale_factors),
    numbers(other.numbers),
    charges(other.charges)
{
}

PyQMEngine &PyQMEngine::operator=(const PyQMEngine &other)
{
    this->callback = other.callback;
    this->cutoff = other.cutoff;
    this->neighbour_list_frequency = other.neighbour_list_frequency;
    this->is_mechanical = other.is_mechanical;
    this->lambda = other.lambda;
    this->atoms = other.atoms;
    this->mm1_to_qm = other.mm1_to_qm;
    this->mm1_to_mm2 = other.mm1_to_mm2;
    this->mm2_atoms = other.mm2_atoms;
    this->bond_scale_factors = other.bond_scale_factors;
    this->numbers = other.numbers;
    this->charges = other.charges;
    return *this;
}

void PyQMEngine::setCallback(PyQMCallback callback)
{
    this->callback = callback;
}

PyQMCallback PyQMEngine::getCallback() const
{
    return this->callback;
}

void PyQMEngine::setLambda(double lambda)
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

double PyQMEngine::getLambda() const
{
    return this->lambda;
}

void PyQMEngine::setCutoff(SireUnits::Dimension::Length cutoff)
{
    this->cutoff = cutoff;
}

SireUnits::Dimension::Length PyQMEngine::getCutoff() const
{
    return this->cutoff;
}

int PyQMEngine::getNeighbourListFrequency() const
{
    return this->neighbour_list_frequency;
}

void PyQMEngine::setNeighbourListFrequency(int neighbour_list_frequency)
{
    // Assume anything less than zero means no neighbour list.
    if (neighbour_list_frequency < 0)
    {
        neighbour_list_frequency = 0;
    }
    this->neighbour_list_frequency = neighbour_list_frequency;
}

bool PyQMEngine::getIsMechanical() const
{
    return this->is_mechanical;
}

void PyQMEngine::setIsMechanical(bool is_mechanical)
{
    this->is_mechanical = is_mechanical;
}

QVector<int> PyQMEngine::getAtoms() const
{
    return this->atoms;
}

void PyQMEngine::setAtoms(QVector<int> atoms)
{
    this->atoms = atoms;
}

boost::tuple<QMap<int, int>, QMap<int, QVector<int>>, QMap<int, double>> PyQMEngine::getLinkAtoms() const
{
    return boost::make_tuple(this->mm1_to_qm, this->mm1_to_mm2, this->bond_scale_factors);
}

void PyQMEngine::setLinkAtoms(
    QMap<int, int> mm1_to_qm,
    QMap<int, QVector<int>> mm1_to_mm2,
    QMap<int, double> bond_scale_factors)
{
    this->mm1_to_qm = mm1_to_qm;
    this->mm1_to_mm2 = mm1_to_mm2;
    this->bond_scale_factors = bond_scale_factors;

    // Build a vector of all of the MM2 atoms.
    this->mm2_atoms.clear();
    for (const auto &mm2 : this->mm1_to_mm2.values())
    {
        this->mm2_atoms.append(mm2);
    }
}

QVector<int> PyQMEngine::getMM2Atoms() const
{
    return this->mm2_atoms;
}

QVector<int> PyQMEngine::getNumbers() const
{
    return this->numbers;
}

void PyQMEngine::setNumbers(QVector<int> numbers)
{
    this->numbers = numbers;
}

QVector<double> PyQMEngine::getCharges() const
{
    return this->charges;
}

void PyQMEngine::setCharges(QVector<double> charges)
{
    this->charges = charges;
}

const char *PyQMEngine::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PyQMEngine>());
}

const char *PyQMEngine::what() const
{
    return PyQMEngine::typeName();
}

boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>
PyQMEngine::call(
    QVector<int> numbers_qm,
    QVector<double> charges_mm,
    QVector<QVector<double>> xyz_qm,
    QVector<QVector<double>> xyz_mm,
    QVector<QVector<double>> cell,
    QVector<int> idx_mm) const
{
    return this->callback.call(numbers_qm, charges_mm, xyz_qm, xyz_mm, cell, idx_mm);
}

QMForce* PyQMEngine::createForce() const
{
    return new PyQMForce(
        this->callback,
        this->cutoff,
        this->neighbour_list_frequency,
        this->is_mechanical,
        this->lambda,
        this->atoms,
        this->mm1_to_qm,
        this->mm1_to_mm2,
        this->bond_scale_factors,
        this->mm2_atoms,
        this->numbers,
        this->charges
    );
}
