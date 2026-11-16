#ifndef _HELPERS_SCOPED_GIL_ACQUIRE_HPP_
#define _HELPERS_SCOPED_GIL_ACQUIRE_HPP_

#include "boost/python.hpp"

namespace SireHelpers
{
    /** Acquire the GIL for the lifetime of this object.
     *
     *  Use this before calling into Python from C++ that may run without the
     *  GIL, e.g. within a ScopedGILRelease or on a worker thread. The GIL is
     *  released on every return path, including during stack unwinding.
     */
    class ScopedGILAcquire
    {
    public:
        ScopedGILAcquire() : state(PyGILState_Ensure())
        {
        }

        ~ScopedGILAcquire()
        {
            PyGILState_Release(state);
        }

    private:
        ScopedGILAcquire(const ScopedGILAcquire &);
        ScopedGILAcquire &operator=(const ScopedGILAcquire &);

        PyGILState_STATE state;
    };
}

#endif
