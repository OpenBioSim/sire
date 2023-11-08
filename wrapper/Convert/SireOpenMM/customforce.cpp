
#include "customforce.h"

#ifdef SIRE_USE_CUSTOMCPPFORCE
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCPPForceImpl.h"
#endif

#include "openmmmolecule.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireOpenMM;

#ifdef SIRE_USE_CUSTOMCPPFORCE
class GridForceImpl : public OpenMM::CustomCPPForceImpl
{
public:
    GridForceImpl(const GridForce &owner)
        : OpenMM::CustomCPPForceImpl(owner),
          owner(owner)
    {
    }

    ~GridForceImpl()
    {
    }

    double computeForce(OpenMM::ContextImpl &context,
                        const std::vector<OpenMM::Vec3> &positions,
                        std::vector<OpenMM::Vec3> &forces)
    {
        // Compute the forces and energy here.  Store the forces into the
        // vector and return the energy.

        // we can get the platform using
        // const auto platform = context.getPlatform().getName();
        // so we could have custom code for different platforms if we wanted

        if (fixed_atoms.empty())
            rebuild_fixed_atoms();

        // now need to calculate the coulomb and LJ forces and energy...
        const auto &params = owner.getParticleParameters();

        const int natoms = positions.size();

        if (natoms != params.count())
            throw SireError::incompatible_error(QObject::tr(
                                                    "Incorrect number of atoms in the grid force. Expected %1, got %2")
                                                    .arg(params.count())
                                                    .arg(natoms),
                                                CODELOC);

        const auto &params_data = params.constData();

        // return the potential energy
        return 0;
    }

    const GridForce &getOwner() const
    {
        return owner;
    }

private:
    void rebuild_fixed_atoms()
    {
        fixed_atoms.clear();

        const auto &atoms = owner.getFieldAtoms();

        const int natoms = atoms.count();

        if (natoms == 0)
            return;

        fixed_atoms = std::vector<std::tuple<OpenMM::Vec3, float, float, float>>(natoms);

        const auto &coords = atoms.getCoords().data();
        const auto &charges = atoms.getCharges().data();
        const auto &sigmas = atoms.getSigmas().data();
        const auto &epsilons = atoms.getEpsilons().data();

        for (int i = 0; i < natoms; ++i)
        {
            fixed_atoms.push_back(std::make_tuple(coords[i],
                                                  charges[i],
                                                  sigmas[i],
                                                  epsilons[i]));
        }
    }

    void rebuild_grid()
    {
    }

    /** All of the atoms data */
    std::vector<std::tuple<OpenMM::Vec3, float, float, float>> fixed_atoms;

    /** The coulomb potential grid */
    std::vector<double> coulomb_grid;

    /** Information about the grid dimensions */
    SireVol::AABox grid_box;

    /** The coulomb cutoff (applies from the center of the grid) */
    double coulomb_cutoff;

    /** The grid spacing */
    double grid_spacing;

    /** The number of grid points along x, y and z */
    unsigned int dimx, dimy, dimz;

    const GridForce &owner;
};
#endif

GridForce::GridForce() : OpenMM::Force()
{
#ifndef SIRE_USE_CUSTOMCPPFORCE
    throw SireError::unsupported(QObject::tr(
                                     "Unable to create a GridForce because OpenMM::CustomCPPForceImpl "
                                     "is not available. You need to use OpenMM 8.1 or later."),
                                 CODELOC);
#endif
}

GridForce::~GridForce()
{
}

OpenMM::ForceImpl *GridForce::createImpl() const
{
#ifdef SIRE_USE_CUSTOMCPPFORCE
    return new GridForceImpl(*this);
#else
    throw SireError::unsupported(QObject::tr(
                                     "Unable to create a GridForce because OpenMM::CustomCPPForceImpl "
                                     "is not available. You need to use OpenMM 8.1 or later."),
                                 CODELOC);
    return 0;
#endif
}

void GridForce::addFieldAtoms(const FieldAtoms &atoms)
{
    field_atoms += atoms;
}

const FieldAtoms &GridForce::getFieldAtoms() const
{
    return field_atoms;
}

void GridForce::addParticle(double charge, double sigma, double epsilon)
{
    params.push_back(std::make_tuple(charge, sigma, epsilon));
}

const QVector<std::tuple<double, double, double>> &GridForce::getParticleParameters() const
{
    return params;
}
