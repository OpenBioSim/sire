
#include "boost/python.hpp"

#include "SireBase/releasegil.h"

#include <QDebug>

static void print_to_python_stdout(const char *text, bool flush)
{
    // thanks to https://github.com/osqp/osqp/blob/25b6b39247954b74bba17667863fe19e9a8b1ade/include/glob_opts.h#L139
    // and https://stackoverflow.com/questions/35745541/how-to-get-printed-output-from-ctypes-c-functions-into-jupyter-ipython-notebook
    PyGILState_STATE gilstate = PyGILState_Ensure();
    PySys_WriteStdout(text);

    if (flush)
    {
        // copied from the python code
        // https://github.com/python/cpython/blob/f25f2e2e8c8e48490d22b0cdf67f575608701f6f/Python/pylifecycle.c#L1589
        // Thanks to https://stackoverflow.com/questions/69247396/python-c-api-how-to-flush-stdout-and-stderr
        PyObject *fout = PySys_GetObject("__stdout__");
        PyObject *tmp;
        int status = 0;

        if (fout != NULL && fout != Py_None)
        {
            tmp = PyObject_CallMethod(fout, "flush", 0);

            if (tmp == NULL)
            {
                PyErr_WriteUnraisable(fout);
                status = -1;
            }
            else
                Py_DECREF(tmp);
        }
    }

    PyGILState_Release(gilstate);
}

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
    void print(const QString &text, bool flush) const
    {
        QByteArray bytes = text.toUtf8();

        // need to split into chunks of 1000 characters
        // because python can only cope with 1000 characters at once
        while (bytes.length() > 1000)
        {
            auto start = bytes.left(999);
            start.append('\0');
            print_to_python_stdout(start.constData(), false);
            bytes = bytes.right(bytes.length() - 999);
        }

        if (bytes.length() > 0)
            print_to_python_stdout(bytes.constData(), flush);
    }

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
