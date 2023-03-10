/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include "simulation.h"

#include "SireCluster/nodes.h"

#include "SireError/errors.h"

using namespace SireMove;
using namespace SireSystem;
using namespace SireCluster;

/** Null constructor */
Simulation::Simulation()
{
}

/** Internal function used to construct from the promise
    of a running simulation */
Simulation::Simulation(const SireCluster::Promise &promise) : sim_promise(promise)
{
}

/** Copy constructor */
Simulation::Simulation(const Simulation &other) : sim_promise(other.sim_promise)
{
}

/** Destructor */
Simulation::~Simulation()
{
}

/** Copy assignment operator */
Simulation &Simulation::operator=(const Simulation &other)
{
    sim_promise = other.sim_promise;
    return *this;
}

/** Comparison operator */
bool Simulation::operator==(const Simulation &other) const
{
    return sim_promise == other.sim_promise;
}

/** Comparison operator */
bool Simulation::operator!=(const Simulation &other) const
{
    return sim_promise != other.sim_promise;
}

/** Run the simulation contained in the simulation WorkPacket 'simpacket'
    on the node 'node' */
Simulation Simulation::run(Node &node, const SimPacket &simpacket)
{
    return Simulation(node.startJob(simpacket));
}

/** Run a simulation consisting of 'nmoves' moves (in 'moves')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation on the node 'node'. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(Node &node, const System &system, const Moves &moves, int nmoves, bool record_stats)
{
    return Simulation::run(node, SimPacket(system, moves, nmoves, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves (in 'move')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation on the node 'node'. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(Node &node, const System &system, const Move &move, int nmoves, bool record_stats)
{
    return Simulation::run(node, SimPacket(system, SameMoves(move), nmoves, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves (in 'moves')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation on the node 'node'. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(Node &node, const System &system, const Moves &moves, int nmoves, int nmoves_per_chunk,
                           bool record_stats)
{
    return Simulation::run(node, SimPacket(system, moves, nmoves, nmoves_per_chunk, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves (in 'move')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation on the node 'node'. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(Node &node, const System &system, const Move &move, int nmoves, int nmoves_per_chunk,
                           bool record_stats)
{
    return Simulation::run(node, SimPacket(system, SameMoves(move), nmoves, nmoves_per_chunk, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves of the system
    in 'simstore' (which also contains the moves) optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation on the node 'node'. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(Node &node, const SimStore &simstore, int nmoves, bool record_stats)
{
    return Simulation::run(node, SimPacket(simstore, nmoves, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves  of the system
    in 'simstore' (which also contains the moves) optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation on the node 'node'. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(Node &node, const SimStore &simstore, int nmoves, int nmoves_per_chunk, bool record_stats)
{
    return Simulation::run(node, SimPacket(simstore, nmoves, nmoves_per_chunk, record_stats));
}

/** Run the simulation contained in the simulation WorkPacket 'simpacket'
    in the current thread */
Simulation Simulation::run(const SimPacket &simpacket)
{
    Nodes nodes;

    ThisThread this_thread = nodes.borrowThisThread();

    if (nodes.isEmpty())
        throw SireError::unavailable_resource(
            QObject::tr("This thread is unavailable for running a simulation. It is already "
                        "busy doing something else!"),
            CODELOC);

    Node node = nodes.getNode();

    Simulation sim = Simulation::run(node, simpacket);

    sim.wait();

    return sim;
}

/** Run a simulation consisting of 'nmoves' moves (in 'moves')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation in the current thread */
Simulation Simulation::run(const System &system, const Moves &moves, int nmoves, bool record_stats)
{
    return Simulation::run(SimPacket(system, moves, nmoves, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves (in 'move')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation in the current thread */
Simulation Simulation::run(const System &system, const Move &move, int nmoves, bool record_stats)
{
    return Simulation::run(SimPacket(system, SameMoves(move), nmoves, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves (in 'moves')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation in the current thread. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(const System &system, const Moves &moves, int nmoves, int nmoves_per_chunk,
                           bool record_stats)
{
    return Simulation::run(SimPacket(system, moves, nmoves, nmoves_per_chunk, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves (in 'move')
    of the System 'system', optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation in the current thread. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(const System &system, const Move &move, int nmoves, int nmoves_per_chunk, bool record_stats)
{
    return Simulation::run(SimPacket(system, SameMoves(move), nmoves, nmoves_per_chunk, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves of the system
    in 'simstore' (which also contains the moves) optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation in the current thread */
Simulation Simulation::run(const SimStore &simstore, int nmoves, bool record_stats)
{
    return Simulation::run(SimPacket(simstore, nmoves, record_stats));
}

/** Run a simulation consisting of 'nmoves' moves of the system
    in 'simstore' (which also contains the moves) optionally recording simulation
    statistics if 'record_stats' is true. This runs the
    simulation in the current thread. This runs 'nmoves_per_chunk'
    moves in every chunk of the simulation. */
Simulation Simulation::run(const SimStore &simstore, int nmoves, int nmoves_per_chunk, bool record_stats)
{
    return Simulation::run(SimPacket(simstore, nmoves, nmoves_per_chunk, record_stats));
}

/** Abort the simulation */
void Simulation::abort()
{
    sim_promise.abort();
}

/** Stop the simulation */
void Simulation::stop()
{
    sim_promise.stop();
}

/** Wait for the simulation to stop running
    (which can be either because it finished, was stopped,
     was aborted or ended in error) */
void Simulation::wait()
{
    sim_promise.wait();
}

/** Wait for the simulation to stop running, or for 'timeout'
    milliseconds to pass, whichever comes soonest. This returns
    whether or not the simulation has stopped */
bool Simulation::wait(int timeout)
{
    return sim_promise.wait(timeout);
}

/** Return whether or not this simulation is running */
bool Simulation::isRunning()
{
    return sim_promise.isRunning();
}

/** Return whether or not this simulation is in an error state */
bool Simulation::isError()
{
    return sim_promise.isError();
}

/** Return whether or not the simulation has finished
    (completed all of the moves) */
bool Simulation::hasFinished()
{
    if (this->isRunning())
        return false;

    else
    {
        // we aren't running any more - lets see what happened
        if (this->isError() or this->wasAborted())
            return false;

        try
        {
            SimPacket sim = this->result();

            return sim.nCompleted() == sim.nMoves();
        }
        catch (...)
        {
            return false;
        }
    }
}

/** Throw any error associated with this simulation - this does
    nothing if we are not in an error state */
void Simulation::throwError()
{
    sim_promise.throwError();
}

/** Return whether or not the simulation was stopped */
bool Simulation::wasStopped()
{
    return sim_promise.wasStopped();
}

/** Return whether or not the simulation was aborted */
bool Simulation::wasAborted()
{
    return sim_promise.wasAborted();
}

/** Return the progress of the simulation (as a percentage) */
float Simulation::progress()
{
    return sim_promise.progress();
}

/** Return the initial input simulation WorkPacket */
SimPacket Simulation::input()
{
    if (sim_promise.isNull())
        return SimPacket();

    else
    {
        WorkPacket initial_packet = sim_promise.input();

        if (initial_packet.isNull())
        {
            throw SireError::program_bug(QObject::tr("How could we lose the input simulation WorkPacket? How has "
                                                     "it become null?"),
                                         CODELOC);
        }

        if (not initial_packet.isA<SimPacket>())
        {
            throw SireError::program_bug(QObject::tr("How could we lose the input simulation WorkPacket? How has "
                                                     "it turned into a %1?")
                                             .arg(initial_packet.base().what()),
                                         CODELOC);
        }

        return initial_packet.asA<SimPacket>();
    }
}

/** Return the simulation WorkPacket from an intermediate point along
    the simulation. This will throw an error if the simulation is in an
    error state, and the initial packet if the simulation
    was aborted */
SimPacket Simulation::interimResult()
{
    if (sim_promise.isNull())
        return SimPacket();

    else
    {
        WorkPacket interim_packet = sim_promise.interimResult();

        if (interim_packet.wasAborted())
        {
            return this->input();
        }
        else if (interim_packet.isError())
        {
            interim_packet.throwError();
        }

        if (interim_packet.isNull())
        {
            throw SireError::program_bug(QObject::tr("How could we lose the interim simulation WorkPacket? How has "
                                                     "it become null?"),
                                         CODELOC);
        }

        if (not interim_packet.isA<SimPacket>())
        {
            throw SireError::program_bug(QObject::tr("How could we lose the interim simulation WorkPacket? How has "
                                                     "it turned into a %1?")
                                             .arg(interim_packet.base().what()),
                                         CODELOC);
        }

        return interim_packet.asA<SimPacket>();
    }
}

/** Return the final result of the simulation. This blocks until
    the simulation has stopped, and will throw an exception if the
    simulation is in an error state. This returns the initial
    simulation WorkPacket if the simulation was aborted */
SimPacket Simulation::result()
{
    if (sim_promise.isNull())
        return SimPacket();

    else
    {
        WorkPacket result_packet = sim_promise.result();

        if (result_packet.wasAborted())
        {
            return this->input();
        }
        else if (result_packet.isError())
        {
            result_packet.throwError();
        }

        if (result_packet.isNull())
        {
            throw SireError::program_bug(QObject::tr("How could we lose the simulation result WorkPacket? How has "
                                                     "it become null?"),
                                         CODELOC);
        }

        if (not result_packet.isA<SimPacket>())
        {
            throw SireError::program_bug(QObject::tr("How could we lose the simulation result WorkPacket? How has "
                                                     "it turned into a %1?")
                                             .arg(result_packet.base().what()),
                                         CODELOC);
        }

        return result_packet.asA<SimPacket>();
    }
}

/** Return the System in the state it was in before the simulation started */
System Simulation::initialSystem()
{
    return this->input().system();
}

/** Return the Moves in the state they were in before the simulation started */
MovesPtr Simulation::initialMoves()
{
    return this->input().moves();
}

/** Return the current state of the System (updated while the simulation
    is running). This will throw an exception if the system hits an
    error state */
System Simulation::interimSystem()
{
    return this->interimResult().system();
}

/** Return the current state of the moves (updated while the simulation
    is running). This will throw an exception if the system hits an
    error state */
MovesPtr Simulation::interimMoves()
{
    return this->interimResult().moves();
}

/** Return the final state of the system after the simulation. This
    blocks until the simulation has finished and will throw an
    exception if the system hits an error state */
System Simulation::system()
{
    return this->result().system();
}

/** Return the final state of the moves after the simulation. This
    blocks until the simulation has finished and will throw an
    exception if the system hits an error state */
MovesPtr Simulation::moves()
{
    return this->result().moves();
}
