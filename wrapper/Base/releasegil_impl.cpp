
#include "boost/python.hpp"

#include "SireBase/releasegil.h"

#include <QDebug>

class ReleaseGIL : public SireBase::detail::ReleaseGILBase
{
public:
    ReleaseGIL() : SireBase::detail::ReleaseGILBase(), thread_state(0)
    {
    }

    ReleaseGIL(PyThreadState *t) : SireBase::detail::ReleaseGILBase(), thread_state(t)
    {
    }

    ~ReleaseGIL()
    {
        if (thread_state)
        {
            // qDebug() << "Re-aquire GIL";
            PyEval_RestoreThread(thread_state);
            thread_state = 0;
        }
    }

    static void register_releasegil()
    {
        SireBase::detail::ReleaseGILBase::registerReleaseGIL(new ReleaseGIL());
    }

protected:
    SireBase::GILHandle releaseGIL() const
    {
        // qDebug() << "Release GIL";
        auto thread_state = PyEval_SaveThread();
        return SireBase::GILHandle(new ReleaseGIL(thread_state));
    }

private:
    PyThreadState *thread_state;
};

void register_releasegil()
{
    ReleaseGIL::register_releasegil();
}
