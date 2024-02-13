
#include "openmmminimise.h"

/**
 *  This code is heavily inspired / adapted from the LocalEnergyMinimizer
 *  in LocalEnergyMinimizer.cpp in the OpenMM source code (version 8.1 beta).
 *
 *  The copyright notice for that file is copied below.
 *
 */

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2020 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/OpenMMException.h"
#include "openmm/Platform.h"
#include "openmm/VerletIntegrator.h"
#include "lbgfs/lbfgs.h"
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "SireError/errors.h"
#include "SireBase/releasegil.h"

#include "SireBase/progressbar.h"
#include "SireUnits/units.h"

#include <QDebug>

namespace SireOpenMM
{
    class MinimizerData
    {
    public:
        MinimizerData(OpenMM::Context &c, double kval, SireBase::ProgressBar &bar,
                      int max_iterations, int ratchet_frequency, double ratchet_scale)
            : context(&c), bar(&bar), starting_k(kval), k(kval), cpu_integrator(1.0), it(0), max_its(max_iterations),
              last_fx(1e100), raw_energy(0.0), num_ratchets(0), num_restarts(1), start_it(0),
              ratchet_scale(ratchet_scale), ratchet_frequency(ratchet_frequency)
        {
            QString platform_name = QString::fromStdString(context->getPlatform().getName());
            check_large_forces = (platform_name == "CUDA" || platform_name == "OpenCL" || platform_name == "HIP" || platform_name == "Metal");
        }

        ~MinimizerData()
        {
        }

        OpenMM::Context &getCpuContext()
        {
            // Get an alternate context that runs on the CPU and doesn't place any limits
            // on the magnitude of forces.
            if (context == 0)
                throw SireError::program_bug(QObject::tr("MinimizerData: The context has been destroyed"),
                                             CODELOC);

            if (cpu_context.get() == 0)
            {
                OpenMM::Platform *cpu_platform;
                try
                {
                    cpu_platform = &OpenMM::Platform::getPlatformByName("CPU");

                    // The CPU context sometimes fails to initialize because it flags
                    // a different number of exclusions / exceptions, despite this
                    // working for all other platforms...
                    cpu_context.reset(new OpenMM::Context(context->getSystem(), cpu_integrator, *cpu_platform));
                }
                catch (...)
                {
                    cpu_platform = &OpenMM::Platform::getPlatformByName("Reference");

                    // In those cases, we will fall back to the Reference platform,
                    // which does work
                    cpu_context.reset(new OpenMM::Context(context->getSystem(), cpu_integrator, *cpu_platform));
                }
                cpu_context->setState(context->getState(OpenMM::State::Positions | OpenMM::State::Velocities | OpenMM::State::Parameters));
            }

            return *cpu_context;
        }

        bool checkLargeForces() const
        {
            return check_large_forces;
        }

        OpenMM::Context &getContext()
        {
            if (context == 0)
                throw SireError::program_bug(
                    QObject::tr("MinimizerData: The context has been destroyed"),
                    CODELOC);

            return *context;
        }

        SireBase::ProgressBar &getProgressBar()
        {
            if (bar == 0)
                throw SireError::program_bug(
                    QObject::tr("MinimizerData: The progress bar has been destroyed"),
                    CODELOC);

            return *bar;
        }

        int nRatchets() const
        {
            return num_ratchets;
        }

        int nRestarts() const
        {
            return num_restarts;
        }

        double getK() const
        {
            return k;
        }

        void scaleK(double factor)
        {
            k *= factor;
        }

        void ratchet()
        {
            num_ratchets++;
            scaleK(2.0);
        }

        void aboutToStart()
        {
            start_it = it;
        }

        qint64 getIteration() const
        {
            return it;
        }

        qint64 getMaxIterations() const
        {
            return max_its;
        }

        void resetIteration()
        {
            it = 0;
            last_fx = 1e10;
            num_restarts += 1;
            k = starting_k;
            num_ratchets = 0;

            for (int i = 0; i < num_restarts - 1; ++i)
            {
                this->ratchet();
            }
        }

        void setIteration(int val)
        {
            it = start_it + val;
        }

        bool hasConverged(double val)
        {
            if (std::abs(val - last_fx) < 1e-3)
            {
                last_fx = val;
                return true;
            }
            else if ((ratchet_frequency > 0) and (it - start_it > ratchet_frequency))
            {
                last_fx = val;
                return true;
            }
            else
            {
                last_fx = val;
                return false;
            }
        }

        void setRawEnergy(double val)
        {
            raw_energy = val;
        }

        double getRawEnergy() const
        {
            return raw_energy;
        }

    private:
        /** This is a pointer to the context being minimised */
        OpenMM::Context *context;

        /** An integrator and context in case we need to create a
         *  CPU or Reference context to calculate forces when they
         *  become too large for the GPU-accelerated contexts
         */
        OpenMM::VerletIntegrator cpu_integrator;
        std::shared_ptr<OpenMM::Context> cpu_context;

        /** The progress bar */
        SireBase::ProgressBar *bar;

        /** The current tolerance */
        double k;

        /** The starting value of the tolerance */
        double starting_k;

        /** The last energy value */
        double last_fx;

        /** The raw energy excluding the fake constraints */
        double raw_energy;

        /** The amount to scale k by at each ratchet */
        double ratchet_scale;

        /** The current iteration */
        qint64 it;

        /** The maximum number of iterations */
        qint64 max_its;

        /** The number of ratchets */
        qint64 num_ratchets;

        /** The number of restarts */
        qint64 num_restarts;

        /** The starting iteration of a cycle */
        qint64 start_it;

        /** The frequency of automatically ratchetting */
        qint64 ratchet_frequency;

        /** Whether or not we need to check for large forces on this platform */
        bool check_large_forces;
    };

    /** Calculate the forces and energies for the passed context, given the
     *  positions in 'positions', returning the energy and storing the
     *  forces in 'g'
     */
    static double computeForcesAndEnergy(OpenMM::Context &context,
                                         const std::vector<OpenMM::Vec3> &positions,
                                         lbfgsfloatval_t *g)
    {
        context.setPositions(positions);
        context.computeVirtualSites();

        OpenMM::State state = context.getState(OpenMM::State::Forces | OpenMM::State::Energy,
                                               false, context.getIntegrator().getIntegrationForceGroups());

        const std::vector<OpenMM::Vec3> &forces = state.getForces();
        const OpenMM::System &system = context.getSystem();

        for (int i = 0; i < forces.size(); i++)
        {
            if (system.getParticleMass(i) == 0)
            {
                g[3 * i] = 0.0;
                g[3 * i + 1] = 0.0;
                g[3 * i + 2] = 0.0;
            }
            else
            {
                g[3 * i] = -forces[i][0];
                g[3 * i + 1] = -forces[i][1];
                g[3 * i + 2] = -forces[i][2];
            }
        }

        return state.getPotentialEnergy();
    }

    void update_bar(const MinimizerData &data, SireBase::ProgressBar &bar)
    {
        auto nrg = SireUnits::Dimension::GeneralUnit(data.getRawEnergy() * SireUnits::kJ_per_mol);

        bar.setProgress(data.getIteration(),
                        QString("Minimising %1.%2: %3 : %4").arg(data.nRestarts()).arg(data.nRatchets()).arg(data.getIteration()).arg(nrg.toString()));
    }

    static int progress(void *instance, const lbfgsfloatval_t *x,
                        const lbfgsfloatval_t *g,
                        const lbfgsfloatval_t fx,
                        const lbfgsfloatval_t xnorm,
                        const lbfgsfloatval_t gnorm,
                        const lbfgsfloatval_t step,
                        int n, int k, int ls)
    {
        MinimizerData *data = reinterpret_cast<MinimizerData *>(instance);

        data->setIteration(k + 1);

        auto &bar = data->getProgressBar();

        update_bar(*data, bar);

        if (data->hasConverged(fx))
        {
            return 1;
        }
        else
        {
            // return 0 to keep going, non-zero to stop
            return data->getIteration() >= data->getMaxIterations();
        }
    }

    static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x,
                                    lbfgsfloatval_t *g, const int n,
                                    const lbfgsfloatval_t step)
    {
        MinimizerData *data = reinterpret_cast<MinimizerData *>(instance);

        OpenMM::Context &context = data->getContext();

        const OpenMM::System &system = context.getSystem();

        int num_particles = system.getNumParticles();

        // Compute the force and energy for this configuration.
        std::vector<OpenMM::Vec3> positions(num_particles);

        for (int i = 0; i < num_particles; i++)
            positions[i] = OpenMM::Vec3(x[3 * i], x[3 * i + 1], x[3 * i + 2]);

        double energy = computeForcesAndEnergy(context, positions, g);

        if (data->checkLargeForces())
        {
            // The CUDA, OpenCL and HIP platforms accumulate forces in fixed point, so they
            // can't handle very large forces.  Check for problematic forces (very large,
            // infinite, or NaN) and if necessary recompute them on the CPU.
            for (int i = 0; i < 3 * num_particles; i++)
            {
                if (!(std::fabs(g[i]) < 2e9))
                {
                    energy = computeForcesAndEnergy(data->getCpuContext(), positions, g);
                    break;
                }
            }
        }

        data->setRawEnergy(energy);

        // Add harmonic forces for any constraints.
        int num_constraints = system.getNumConstraints();
        double k = data->getK();

        for (int i = 0; i < num_constraints; ++i)
        {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(i, particle1, particle2, distance);

            OpenMM::Vec3 delta = positions[particle2] - positions[particle1];

            double r2 = delta.dot(delta);
            double r = std::sqrt(r2);
            delta *= 1 / r;
            double dr = r - distance;
            double kdr = k * dr;

            energy += 0.5 * kdr * dr;

            if (system.getParticleMass(particle1) != 0)
            {
                g[3 * particle1] -= kdr * delta[0];
                g[3 * particle1 + 1] -= kdr * delta[1];
                g[3 * particle1 + 2] -= kdr * delta[2];
            }
            if (system.getParticleMass(particle2) != 0)
            {
                g[3 * particle2] += kdr * delta[0];
                g[3 * particle2 + 1] += kdr * delta[1];
                g[3 * particle2 + 2] += kdr * delta[2];
            }
        }

        auto &bar = data->getProgressBar();
        update_bar(*data, bar);

        return energy;
    }

    class LBFGSRaii
    {
    public:
        LBFGSRaii(lbfgsfloatval_t *_x) : x(_x) {}
        ~LBFGSRaii()
        {
            lbfgs_free(x);
        }

    private:
        lbfgsfloatval_t *x;
    };

    void minimise_openmm_context(OpenMM::Context &context,
                                 double tolerance, int max_iterations,
                                 int max_restarts, int max_ratchets,
                                 int ratchet_frequency,
                                 double starting_k, double ratchet_scale)
    {
        if (max_iterations < 0)
        {
            max_iterations = std::numeric_limits<int>::max();
        }

        auto gil = SireBase::release_gil();

        const OpenMM::System &system = context.getSystem();

        int num_particles = system.getNumParticles();

        double constraint_tol = context.getIntegrator().getConstraintTolerance();

        double working_constraint_tol = std::max(1e-4, constraint_tol);

        // this is about 1e6 - it needs to be large to keep constraints in place
        double k = starting_k / working_constraint_tol;

        SireBase::ProgressBar bar("Minimising: initialise", max_iterations);
        bar.setSpeedUnit("steps / s");

        bar = bar.enter();

        MinimizerData data(context, k, bar, max_iterations, ratchet_frequency,
                           ratchet_scale);
        lbfgsfloatval_t *x = lbfgs_malloc(num_particles * 3);

        if (x == 0)
            throw SireError::unavailable_resource(
                QObject::tr("minimise_openmm_context: Failed to allocate memory"),
                CODELOC);

        // this will automatically free x when it goes out of scope
        LBFGSRaii raii(x);

        std::vector<OpenMM::Vec3> initial_pos;
        double energy_before = 1e10;
        bool is_success = true;

        while (data.getIteration() < data.getMaxIterations())
        {
            if (not is_success)
                break;

            try
            {
                energy_before = context.getState(OpenMM::State::Energy).getPotentialEnergy();
                bar.silentTick();

                // Initialize the minimizer.
                lbfgs_parameter_t param;
                lbfgs_parameter_init(&param);

                if (!context.getPlatform().supportsDoublePrecision())
                    param.xtol = 1e-7;

                param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;

                // Make sure the initial configuration satisfies all constraints.
                context.applyConstraints(working_constraint_tol);

                // Record the initial positions and determine a normalization constant for scaling the tolerance.
                initial_pos = context.getState(OpenMM::State::Positions).getPositions();

                double norm = 0.0;
                for (int i = 0; i < num_particles; i++)
                {
                    x[3 * i] = initial_pos[i][0];
                    x[3 * i + 1] = initial_pos[i][1];
                    x[3 * i + 2] = initial_pos[i][2];
                    norm += initial_pos[i].dot(initial_pos[i]);
                }

                norm /= num_particles;
                norm = (norm < 1 ? 1 : sqrt(norm));
                param.epsilon = tolerance / norm;

                double prev_max_error = 1e10;

                // Repeatedly minimize, steadily increasing the strength of the springs until all constraints are satisfied.
                while (data.getIteration() < data.getMaxIterations())
                {
                    param.max_iterations = data.getMaxIterations() - data.getIteration();
                    lbfgsfloatval_t fx; // final energy
                    auto last_it = data.getIteration();

                    // Perform the minimization.
                    data.aboutToStart();
                    lbfgs(num_particles * 3, x, &fx, evaluate, progress, &data, &param);

                    if (std::isinf(fx) || std::isnan(fx))
                    {
                        // we've hit a singularity, so we need to restart
                        throw std::exception();
                    }
                    else if (data.getIteration() == last_it)
                    {
                        // we didn't make any progress, so ratchet up
                        if (data.nRatchets() < max_ratchets)
                        {
                            data.ratchet();
                        }
                        break;
                    }

                    // Check whether all constraints are satisfied.
                    std::vector<OpenMM::Vec3> positions = context.getState(OpenMM::State::Positions).getPositions();

                    int num_constraints = system.getNumConstraints();

                    double max_error = 0.0;

                    for (int i = 0; i < num_constraints; ++i)
                    {
                        int particle1, particle2;
                        double distance;

                        system.getConstraintParameters(i, particle1, particle2, distance);

                        OpenMM::Vec3 delta = positions[particle2] - positions[particle1];

                        double r = std::sqrt(delta.dot(delta));

                        double error = std::fabs(r - distance) / distance;

                        if (error > max_error)
                            max_error = error;
                    }

                    if (max_error <= working_constraint_tol)
                        break; // All constraints are satisfied.

                    if (data.nRatchets() < max_ratchets)
                    {
                        data.ratchet();
                    }

                    if (max_error > 1000 * working_constraint_tol)
                    {
                        // we need to be really sure that this is the right thing to do,
                        // because it could mean that minimisation exits with a structure
                        // that is very close to the un-minimised input structure

                        if (data.nRestarts() >= max_restarts)
                        {
                            // we've tried too many times to get the constraints to be satisfied
                            // so we should just give up
                            is_success = false;
                            break;
                        }

                        // restart from scratch...
                        data.resetIteration();

                        context.setPositions(initial_pos);

                        // We've gotten far enough from a valid state that we might have trouble getting
                        // back, so reset to the original positions.
                        for (int i = 0; i < num_particles; ++i)
                        {
                            x[3 * i] = initial_pos[i][0];
                            x[3 * i + 1] = initial_pos[i][1];
                            x[3 * i + 2] = initial_pos[i][2];
                        }
                    }

                    if (max_error >= prev_max_error)
                        break; // Further tightening the springs doesn't seem to be helping, so just give up.

                    prev_max_error = max_error;
                }
            }
            catch (...)
            {
                // reset to the old coordinates and keep trying
                if (data.nRestarts() > max_restarts)
                {
                    // we've tried too many times to get the constraints to be satisfied
                    // so we should just give up
                    is_success = false;
                    break;
                }

                // restart from scratch...
                data.resetIteration();

                context.setPositions(initial_pos);

                // We've gotten far enough from a valid state that we might have trouble getting
                // back, so reset to the original positions.
                for (int i = 0; i < num_particles; ++i)
                {
                    x[3 * i] = initial_pos[i][0];
                    x[3 * i + 1] = initial_pos[i][1];
                    x[3 * i + 2] = initial_pos[i][2];
                }
            }

            try
            {
                auto energy_after = context.getState(OpenMM::State::Energy).getPotentialEnergy();

                if (std::isinf(energy_after) or std::isnan(energy_after) or energy_after > 1e20)
                {
                    // something has gone wrong, so we need to restart
                    throw std::exception();
                }

                // do a final constraint projection to make sure they are satisfied
                // to the full precision requested by the user.
                context.applyConstraints(working_constraint_tol);

                if (std::abs(energy_after - energy_before) < 100.0)
                {
                    // only 100 kJ/mol difference, so not much more to
                    // be gained by further minimisation (we don't want to
                    // get stuck in a cycle caused by re-application of the
                    // constraints)
                    // (minimisation is only really to remove bad contacts,
                    //  so looking for big changes in energy)
                    break;
                }
            }
            catch (...)
            {
                // reset to the old coordinates and keep trying
                if (data.nRestarts() >= max_restarts)
                {
                    // we've tried too many times to get the constraints to be satisfied
                    // so we should just give up
                    is_success = false;
                    break;
                }

                // restart from scratch...
                data.resetIteration();

                context.setPositions(initial_pos);

                // We've gotten far enough from a valid state that we might have trouble getting
                // back, so reset to the original positions.
                for (int i = 0; i < num_particles; ++i)
                {
                    x[3 * i] = initial_pos[i][0];
                    x[3 * i + 1] = initial_pos[i][1];
                    x[3 * i + 2] = initial_pos[i][2];
                }
            }
        }

        if (is_success)
            bar.success();
        else
        {
            bar.failure("Minimisation failed! Could not satisfy constraints!");
            throw SireError::invalid_state(QObject::tr(
                                               "Despite repeated attempts, the minimiser could not minimise the system "
                                               "while simultaneously satisfying the constraints."),
                                           CODELOC);
        }
    }

} // namespace SireOpenMM
