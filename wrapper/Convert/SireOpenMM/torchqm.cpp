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

#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"

#include <torch/csrc/autograd/autograd.h>

#include "SireError/errors.h"
#include "SireMaths/vector.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireVol/triclinicbox.h"

#include "torchqm.h"

using namespace SireMaths;
using namespace SireOpenMM;
using namespace SireStream;
using namespace SireVol;

// The delta used to place virtual point charges either side of the MM2
// atoms, in nanometers.
static const double VIRTUAL_PC_DELTA = 0.01;

// Conversion factor from Hartree to kJ/mol.
static const double HARTREE_TO_KJ_MOL = 2625.499638755248;

/////////
///////// Implementation of TorchQMForce
/////////

static const RegisterMetaType<TorchQMForce> r_torchqmforce(NO_ROOT);

QDataStream &operator<<(QDataStream &ds, const TorchQMForce &torchqmforce)
{
    writeHeader(ds, r_torchqmforce, 1);

    SharedDataStream sds(ds);

    sds << torchqmforce.module_path << torchqmforce.cutoff << torchqmforce.neighbour_list_frequency
        << torchqmforce.lambda << torchqmforce.atoms << torchqmforce.mm1_to_qm
        << torchqmforce.mm1_to_mm2 << torchqmforce.bond_scale_factors << torchqmforce.mm2_atoms
        << torchqmforce.numbers << torchqmforce.charges;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, TorchQMForce &torchqmforce)
{
    VersionID v = readHeader(ds, r_torchqmforce);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> torchqmforce.module_path >> torchqmforce.cutoff >> torchqmforce.neighbour_list_frequency
            >> torchqmforce.lambda >> torchqmforce.atoms >> torchqmforce.mm1_to_qm
            >> torchqmforce.mm1_to_mm2 >> torchqmforce.bond_scale_factors >> torchqmforce.mm2_atoms
            >> torchqmforce.numbers >> torchqmforce.charges;

        // Re-load the Torch module.
        torchqmforce.setModulePath(torchqmforce.getModulePath());
    }
    else
        throw version_error(v, "1", r_torchqmforce, CODELOC);

    return ds;
}

TorchQMForce::TorchQMForce()
{
}

TorchQMForce::TorchQMForce(
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
    QVector<double> charges) :
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
    // Try to load the Torch module.
    this->setModulePath(module_path);
}

TorchQMForce::TorchQMForce(const TorchQMForce &other) :
    module_path(other.module_path),
    torch_module(other.torch_module),
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

TorchQMForce &TorchQMForce::operator=(const TorchQMForce &other)
{
    this->module_path = other.module_path;
    this->torch_module = other.torch_module;
    this->cutoff = other.cutoff;
    this->neighbour_list_frequency = other.neighbour_list_frequency;
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

void TorchQMForce::setModulePath(QString module_path)
{
    // Try to load the Torch module.
    try
    {
        torch::jit::getProfilingMode() = false;
        torch::jit::setGraphExecutorOptimize(false);
        this->torch_module = torch::jit::load(module_path.toStdString());
        this->torch_module.eval();
    }
    catch (const c10::Error& e)
    {
        throw SireError::io_error(
                QObject::tr(
                    "Unable to load the TorchScript module '%1'. The error was '%2'.")
                    .arg(module_path).arg(e.what()),
                    CODELOC);
    }

    this->module_path = module_path;
}

QString TorchQMForce::getModulePath() const
{
    return this->module_path;
}

torch::jit::script::Module TorchQMForce::getTorchModule() const
{
    return this->torch_module;
}

void TorchQMForce::setLambda(double lambda)
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

double TorchQMForce::getLambda() const
{
    return this->lambda;
}

SireUnits::Dimension::Length TorchQMForce::getCutoff() const
{
    return this->cutoff;
}

int TorchQMForce::getNeighbourListFrequency() const
{
    return this->neighbour_list_frequency;
}

bool TorchQMForce::getIsMechanical() const
{
    return this->is_mechanical;;
}

QVector<int> TorchQMForce::getAtoms() const
{
    return this->atoms;
}

boost::tuple<QMap<int, int>, QMap<int, QVector<int>>, QMap<int, double>> TorchQMForce::getLinkAtoms() const
{
    return boost::make_tuple(this->mm1_to_qm, this->mm1_to_mm2, this->bond_scale_factors);
}

QVector<int> TorchQMForce::getMM2Atoms() const
{
    return this->mm2_atoms;
}

QVector<int> TorchQMForce::getNumbers() const
{
    return this->numbers;
}

QVector<double> TorchQMForce::getCharges() const
{
    return this->charges;
}

const char *TorchQMForce::typeName()
{
    return QMetaType::typeName(qMetaTypeId<TorchQMForce>());
}

const char *TorchQMForce::what() const
{
    return TorchQMForce::typeName();
}

/////////
///////// OpenMM Serialization
/////////

namespace OpenMM
{
    class TorchQMForceProxy : public SerializationProxy {
        public:
            TorchQMForceProxy() : SerializationProxy("TorchQMForce")
            {
            };

            void serialize(const void* object, SerializationNode& node) const
            {
                // Serialize the object.
                QByteArray data;
                QDataStream ds(&data, QIODevice::WriteOnly);
                TorchQMForce torchqmforce = *static_cast<const TorchQMForce*>(object);
                ds << torchqmforce;

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
                TorchQMForce torchqmforce;

                try
                {
                    ds >> torchqmforce;
                }
                catch (...)
                {
                    throw OpenMM::OpenMMException("Unable deserialize TorchQMForce");
                }

                return new TorchQMForce(torchqmforce);
            };
    };

    // Register the TorchQMForce serialization proxy.
    extern "C" void registerTorchQMSerializationProxies() {
        SerializationProxy::registerProxy(typeid(TorchQMForce), new TorchQMForceProxy());
    }
};

/////////
///////// Implementation of TorchQMForceImpl
/////////

OpenMM::ForceImpl *TorchQMForce::createImpl() const
{
#ifdef SIRE_USE_CUSTOMCPPFORCE
    return new TorchQMForceImpl(*this);
#else
    throw SireError::unsupported(QObject::tr(
                                     "Unable to create an TorchQMForceImpl because OpenMM::CustomCPPForceImpl "
                                     "is not available. You need to use OpenMM 8.1 or later."),
                                 CODELOC);
    return 0;
#endif
}

TorchQMForceImpl::TorchQMForceImpl(const TorchQMForce &owner) :
    OpenMM::CustomCPPForceImpl(owner),
    owner(owner)
{
    this->torch_module = owner.getTorchModule();
}

TorchQMForceImpl::~TorchQMForceImpl()
{
}

const TorchQMForce &TorchQMForceImpl::getOwner() const
{
    return this->owner;
}

double TorchQMForceImpl::computeForce(
    OpenMM::ContextImpl &context,
    const std::vector<OpenMM::Vec3> &positions,
    std::vector<OpenMM::Vec3> &forces)
{
#ifdef SIRE_USE_CUSTOMCPPFORCE
    // Get the platform name from the context.
    const auto platform = context.getPlatform().getName();

    // Set the Torch device.
    auto device = torch::kCUDA;
    if (platform != "CUDA")
    {
        device = torch::kCPU;
    }

    // If this is the first step, then setup information for the neighbour list.
    if (this->step_count == 0)
    {
        // Move the Torch module to the correct device.
        this->torch_module.to(device);

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
    QVector<Vector> xyz_qm_vec(qm_atoms.size());
    std::vector<float> xyz_qm(3*qm_atoms.size());

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
        xyz_qm[3*i] = qm_vec[0];
        xyz_qm[3*i+1] = qm_vec[1];
        xyz_qm[3*i+2] = qm_vec[2];
        center += qm_vec;
        i++;
    }
    center /= i;

    // Initialise a vector to hold the current positions and charges for the MM atoms.
    std::vector<float> xyz_mm;
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
        std::vector<float> xyz_virtual;
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
                            xyz_mm.push_back(mm_vec[0]);
                            xyz_mm.push_back(mm_vec[1]);
                            xyz_mm.push_back(mm_vec[2]);

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
                        xyz_mm.push_back(mm_vec[0]);
                        xyz_mm.push_back(mm_vec[1]);
                        xyz_mm.push_back(mm_vec[2]);

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
            xyz_qm.push_back(link_vec[0]);
            xyz_qm.push_back(link_vec[1]);
            xyz_qm.push_back(link_vec[2]);

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
                xyz_mm.push_back(mm2_vec[0]);
                xyz_mm.push_back(mm2_vec[1]);
                xyz_mm.push_back(mm2_vec[2]);

                // Add the charge and index.
                charges_mm.append(this->owner.getCharges()[mm2_idx] + frac_charge);
                idx_mm.append(mm2_idx);

                // Now add the virtual point charges.

                // Compute the normal vector from the MM1 to MM2 atom.
                const auto normal = (mm2_vec - mm1_vec).normalise();

                // Positive direction. (Away from MM1 atom.)
                auto xyz = mm2_vec + VIRTUAL_PC_DELTA*normal;
                xyz_virtual.push_back(xyz[0]);
                xyz_virtual.push_back(xyz[1]);
                xyz_virtual.push_back(xyz[2]);
                charges_virtual.append(-frac_charge);

                // Negative direction (Towards MM1 atom.)
                xyz = mm2_vec - VIRTUAL_PC_DELTA*normal;
                xyz_virtual.push_back(xyz[0]);
                xyz_virtual.push_back(xyz[1]);
                xyz_virtual.push_back(xyz[2]);
                charges_virtual.append(frac_charge);
            }
        }

        // Store the current number of MM atoms.
        num_mm = charges_mm.size();

        // If there are any virtual point charges, then add to the MM positions
        // and charges.
        if (xyz_virtual.size() > 0)
        {
            xyz_mm.reserve(xyz_mm.size() + xyz_virtual.size());
            xyz_mm.insert(xyz_mm.end(), xyz_virtual.begin(), xyz_virtual.end());
            charges_mm.append(charges_virtual);
        }

        // Update the maximum number of MM atoms that we've seen.
        if (charges_mm.size() > this->max_num_mm)
        {
            this->max_num_mm = charges_mm.size();
        }
        else
        {
            // Resize the charges and positions vectors to the maximum number of MM atoms.
            // This is to try to preserve a static compute graph to avoid re-jitting.
            charges_mm.resize(this->max_num_mm);
            xyz_mm.resize(3*this->max_num_mm);
        }
    }

    // Convert input to Torch tensors.

    // MM charges.
    torch::Tensor charges_mm_torch = torch::from_blob(charges_mm.data(), {charges_mm.size()},
        torch::TensorOptions().dtype(torch::kFloat64))
                              .to(torch::kFloat32).to(device);

    // Atomic numbers.
    torch::Tensor atomic_numbers_torch = torch::from_blob(numbers.data(), {numbers.size()},
        torch::TensorOptions().dtype(torch::kInt32))
                              .to(torch::kInt64).to(device);

    // QM positions.
    torch::Tensor xyz_qm_torch = torch::from_blob(xyz_qm.data(), {numbers.size(), 3},
        torch::TensorOptions().dtype(torch::kFloat32))
                              .to(device);
    xyz_qm_torch.requires_grad_(true);

    // MM positions.
    torch::Tensor xyz_mm_torch = torch::from_blob(xyz_mm.data(), {charges_mm.size(), 3},
        torch::TensorOptions().dtype(torch::kFloat32))
                              .to(device);
    xyz_mm_torch.requires_grad_(true);

    // Create the input vector.
    auto input = std::vector<torch::jit::IValue>{
        atomic_numbers_torch,
        charges_mm_torch,
        xyz_qm_torch,
        xyz_mm_torch
    };

    // Compute the energies.
    auto energies = this->torch_module.forward(input).toTensor();

    // Store the sum of the energy in kJ.
    const auto energy = energies.sum().item<double>() * HARTREE_TO_KJ_MOL;

    // Compute the gradients.
    const auto gradients = torch::autograd::grad({energies.sum()}, {xyz_qm_torch, xyz_mm_torch});

    // Compute the forces, converting from Hatree/Anstrom to kJ/mol/nm.
    const auto forces_qm = -(gradients[0] * HARTREE_TO_KJ_MOL * 10).cpu();
    const auto forces_mm = -(gradients[1] * HARTREE_TO_KJ_MOL * 10).cpu();

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
        // Fall back on the lambda value stored in the TorchQMForce object.
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

    // Flatten the forces to std::vector<float>.
    std::vector<float> forces_qm_flat(
        forces_qm.data_ptr<float>(),
        forces_qm.data_ptr<float>() + forces_qm.numel());
    std::vector<float> forces_mm_flat(
        forces_mm.data_ptr<float>(),
        forces_mm.data_ptr<float>() + forces_mm.numel());

    // First the QM atoms.
    for (int i=0; i<qm_atoms.size(); i++)
    {
        // Get the index of the atom.
        const auto idx = qm_atoms[i];

        // Convert to OpenMM format.
        OpenMM::Vec3 omm_force(
            forces_qm_flat[3*i],
            forces_qm_flat[3*i+1],
            forces_qm_flat[3*i+2]);

        // Update the force vector.
        forces[idx] = lambda * omm_force;
    }

    // Now the MM atoms.
    for (int i=0; i<num_mm; i++)
    {
        // Get the index of the atom.
        const auto idx = idx_mm[i];

        // Convert to OpenMM format.
        OpenMM::Vec3 omm_force(
            forces_mm_flat[3*i],
            forces_mm_flat[3*i+1],
            forces_mm_flat[3*i+2]);

        // Update the force vector.
        forces[idx] = lambda * omm_force;
    }

    // Update the step count.
    this->step_count++;

    // Finally, return the energy.
    return lambda * energy;
#endif
}

/////////
///////// Implementation of TorchQMEngine
/////////

TorchQMEngine::TorchQMEngine() : ConcreteProperty<TorchQMEngine, QMEngine>()
{
    // Register the serialization proxies.
    OpenMM::registerTorchQMSerializationProxies();
}

TorchQMEngine::TorchQMEngine(
    QString module_path,
    SireUnits::Dimension::Length cutoff,
    int neighbour_list_frequency,
    bool is_mechanical,
    double lambda) :
    ConcreteProperty<TorchQMEngine, QMEngine>(),
    module_path(module_path),
    cutoff(cutoff),
    neighbour_list_frequency(neighbour_list_frequency),
    is_mechanical(is_mechanical),
    lambda(lambda)
{
    // Register the serialization proxies.
    OpenMM::registerTorchQMSerializationProxies();

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

TorchQMEngine::TorchQMEngine(const TorchQMEngine &other) :
    module_path(other.module_path),
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

TorchQMEngine &TorchQMEngine::operator=(const TorchQMEngine &other)
{
    this->module_path = other.module_path;
    this->cutoff = other.cutoff;
    this->neighbour_list_frequency = other.neighbour_list_frequency;
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

void TorchQMEngine::setModulePath(QString module_path)
{
    this->module_path = module_path;
}

QString TorchQMEngine::getModulePath() const
{
    return this->module_path;
}

void TorchQMEngine::setLambda(double lambda)
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

double TorchQMEngine::getLambda() const
{
    return this->lambda;
}

void TorchQMEngine::setCutoff(SireUnits::Dimension::Length cutoff)
{
    this->cutoff = cutoff;
}

SireUnits::Dimension::Length TorchQMEngine::getCutoff() const
{
    return this->cutoff;
}

int TorchQMEngine::getNeighbourListFrequency() const
{
    return this->neighbour_list_frequency;
}

void TorchQMEngine::setNeighbourListFrequency(int neighbour_list_frequency)
{
    // Assume anything less than zero means no neighbour list.
    if (neighbour_list_frequency < 0)
    {
        neighbour_list_frequency = 0;
    }
    this->neighbour_list_frequency = neighbour_list_frequency;
}

bool TorchQMEngine::getIsMechanical() const
{
    return this->is_mechanical;
}

void TorchQMEngine::setIsMechanical(bool is_mechanical)
{
    this->is_mechanical = is_mechanical;
}

QVector<int> TorchQMEngine::getAtoms() const
{
    return this->atoms;
}

void TorchQMEngine::setAtoms(QVector<int> atoms)
{
    this->atoms = atoms;
}

boost::tuple<QMap<int, int>, QMap<int, QVector<int>>, QMap<int, double>> TorchQMEngine::getLinkAtoms() const
{
    return boost::make_tuple(this->mm1_to_qm, this->mm1_to_mm2, this->bond_scale_factors);
}

void TorchQMEngine::setLinkAtoms(
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

QVector<int> TorchQMEngine::getMM2Atoms() const
{
    return this->mm2_atoms;
}

QVector<int> TorchQMEngine::getNumbers() const
{
    return this->numbers;
}

void TorchQMEngine::setNumbers(QVector<int> numbers)
{
    this->numbers = numbers;
}

QVector<double> TorchQMEngine::getCharges() const
{
    return this->charges;
}

void TorchQMEngine::setCharges(QVector<double> charges)
{
    this->charges = charges;
}

const char *TorchQMEngine::typeName()
{
    return QMetaType::typeName(qMetaTypeId<TorchQMEngine>());
}

const char *TorchQMEngine::what() const
{
    return TorchQMEngine::typeName();
}

QMForce* TorchQMEngine::createForce() const
{
    return new TorchQMForce(
        this->module_path,
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
