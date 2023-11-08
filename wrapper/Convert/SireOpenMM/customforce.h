#ifndef SIRE_OPENMM_CUSTOMFORCE_H
#define SIRE_OPENMM_CUSTOMFORCE_H

#include "openmm.h"
#include "openmm/Force.h"

#include "openmmmolecule.h"

#include "SireVol/aabox.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    class GridForce : public OpenMM::Force
    {
    public:
        GridForce();
        ~GridForce();

        void addFieldAtoms(const FieldAtoms &atoms);

        const FieldAtoms &getFieldAtoms() const;

        void addParticle(double charge, double sigma, double epsilon);

        const QVector<std::tuple<double, double, double>> &getParticleParameters() const;

    protected:
        OpenMM::ForceImpl *createImpl() const;

    private:
        /** All of the field atoms */
        FieldAtoms field_atoms;

        /** All of the particle parameters */
        QVector<std::tuple<double, double, double>> params;
    };
}

SIRE_END_HEADER

#endif
