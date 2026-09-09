#ifndef _HELPERS_SCOPED_GIL_RELEASE_HPP_
#define _HELPERS_SCOPED_GIL_RELEASE_HPP_

#include "boost/python.hpp"

namespace SireHelpers
{
    /** Release the GIL for the lifetime of this object.
     *
     *  Use this in a hand-written wrapper function for a call that is hot,
     *  long-running and never re-enters Python. Unlike bp::release_gil_policy,
     *  the GIL is restored during stack unwinding, so this is safe for
     *  functions that throw and for functions with default arguments.
     */
    class ScopedGILRelease
    {
    public:
        ScopedGILRelease() : thread_state(PyEval_SaveThread())
        {
        }

        ~ScopedGILRelease()
        {
            PyEval_RestoreThread(thread_state);
        }

    private:
        ScopedGILRelease(const ScopedGILRelease &);
        ScopedGILRelease &operator=(const ScopedGILRelease &);

        PyThreadState *thread_state;
    };
}

#endif
