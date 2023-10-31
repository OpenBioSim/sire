
#include "customforce.h"

#ifdef SIRE_USE_CUSTOMCPPFORCE
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCPPForceImpl.h"
#endif

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

        // return the potential energy
        return 0;
    }

    const GridForce &getOwner() const
    {
        return owner;
    }

private:
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
