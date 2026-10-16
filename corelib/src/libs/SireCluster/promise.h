/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
  *
  *  This program is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU Lesser General Public License as published
  *  by the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU Lesser General Public License for more details.
  *
  *  You should have received a copy of the GNU Lesser General Public License
  *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  *
  *  For full details of the license please see the LICENSE file
  *  that should have come with this distribution.
  *
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIRECLUSTER_PROMISE_H
#define SIRECLUSTER_PROMISE_H

#include <QByteArray>

#include "SireStream/streamdata.hpp"

SIRE_BEGIN_HEADER

namespace SireCluster
{

    class Node;
    class WorkPacket;

    namespace detail
    {
        class PromisePvt;
    }

    /** This class provides a handle to the (future) result of
        a piece of work that is being conducted on a node
        in the cluster

        @author Christopher Woods
    */
    class SIRECLUSTER_EXPORT Promise
    {

        friend class Node;

    public:
        Promise();

        Promise(const Promise &other);

        ~Promise();

        Promise &operator=(const Promise &other);

        bool operator==(const Promise &other) const;
        bool operator!=(const Promise &other) const;

        void abort();
        void stop();

        void wait();
        bool wait(int timeout);

        bool isNull() const;

        bool isRunning();

        bool isError();
        void throwError();

        bool wasStopped();

        bool wasAborted();

        float progress();

        WorkPacket input();
        WorkPacket interimResult();
        WorkPacket result();

    protected:
        Promise(const Node &node, const WorkPacket &initial_workpacket);

        /** Pointer to the private implementation */
        std::shared_ptr<detail::PromisePvt> d;
    };

} // namespace SireCluster

SIRE_EXPOSE_CLASS(SireCluster::Promise)

SIRE_END_HEADER

#endif
