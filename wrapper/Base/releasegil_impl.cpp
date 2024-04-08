
#include "boost/python.hpp"

#include "SireBase/releasegil.h"
#include "SireBase/console.h"

#include <QMutex>

#include <QDebug>

static bool _is_ipython = false;

bool is_ipython()
{
    return _is_ipython;
}

void set_is_ipython(bool value)
{
    _is_ipython = value;
}

static PyObject *sys_stdout = 0;

static bool sys_stdout_is_atty()
{
    if (is_ipython())
        // ipython uses a weird sys.stdout, which is not a tty, but should
        // be treated as interactive
        return true;

    PyGILState_STATE gilstate = PyGILState_Ensure();

    // copied from the python code
    // https://github.com/python/cpython/blob/f25f2e2e8c8e48490d22b0cdf67f575608701f6f/Python/pylifecycle.c#L1589
    // Thanks to https://stackoverflow.com/questions/69247396/python-c-api-how-to-flush-stdout-and-stderr
    if (sys_stdout == 0)
    {
        sys_stdout = PySys_GetObject("__stdout__");

        if (sys_stdout == 0)
        {
            qWarning() << "CANNOT ACCESS SYS_STDOUT!";
            PyGILState_Release(gilstate);
            return false;
        }
    }

    bool ret = false;

    if (sys_stdout != 0 and sys_stdout != Py_None)
    {
        PyObject *result = PyObject_CallMethod(sys_stdout, "isatty", 0);

        if (result == 0)
        {
            qWarning() << "UNABLE TO DETERMINE IF STDOUT IS A TTY";
            ret = false;
        }
        else
        {
            ret = (result == Py_True);
            Py_DECREF(result);
        }
    }

    PyGILState_Release(gilstate);

    return ret;
}

static void sys_stdout_write(const char *text, bool flush)
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
        if (sys_stdout == 0)
        {
            sys_stdout = PySys_GetObject("__stdout__");

            if (sys_stdout == 0)
            {
                qWarning() << "CANNOT ACCESS SYS_STDOUT!";
                PyGILState_Release(gilstate);
                return;
            }
        }

        if (sys_stdout != 0 and sys_stdout != Py_None)
        {
            PyObject *result = PyObject_CallMethod(sys_stdout, "flush", 0);

            if (result == 0)
            {
                qWarning() << "UNABLE TO FLUSH STDOUT!";
            }
            else
                Py_DECREF(result);
        }
    }

    PyGILState_Release(gilstate);
}

static PyObject *clear_output = 0;

static void ipython_display_clear(bool wait)
{
    PyGILState_STATE gilstate = PyGILState_Ensure();

    if (clear_output == 0)
    {
        PyObject *ipython_name = PyUnicode_FromString("IPython");
        PyObject *ipython = PyImport_GetModule(ipython_name);
        Py_DECREF(ipython_name);

        if (ipython == 0)
        {
            qWarning() << "COULD NOT IMPORT IPYTHON!";
            PyGILState_Release(gilstate);
            return;
        }

        PyObject *display_name = PyUnicode_FromString("display");
        PyObject *ipython_display = PyObject_GetAttr(ipython, display_name);
        Py_DECREF(display_name);
        Py_DECREF(ipython);

        if (ipython_display == 0)
        {
            qWarning() << "COULD NOT IMPORT IPYTHON.DISPLAY";
            PyGILState_Release(gilstate);
            return;
        }

        PyObject *clear_name = PyUnicode_FromString("clear_output");
        clear_output = PyObject_GetAttr(ipython_display, clear_name);
        Py_DECREF(clear_name);
        Py_DECREF(ipython_display);

        if (clear_output == 0)
        {
            qWarning() << "COULD NOT IMPORT IPYTHON.DISPLAY";
            PyGILState_Release(gilstate);
            return;
        }
    }

    PyObject *result;

    if (wait)
    {
        Py_INCREF(Py_True);
        result = PyObject_CallFunctionObjArgs(clear_output, Py_True, 0);
        Py_DECREF(Py_True);
    }
    else
    {
        Py_INCREF(Py_False);
        result = PyObject_CallFunctionObjArgs(clear_output, Py_False, 0);
        Py_DECREF(Py_False);
    }

    if (result == 0)
    {
        qWarning() << "FAILED TO CALL IPYTHON.DISPLAY.CLEAR_OUTPUT";
    }
    else
    {
        Py_DECREF(result);
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
    bool is_ipython() const
    {
        return ::is_ipython();
    }

    void ipython_clear(bool wait) const
    {
        ::ipython_display_clear(wait);
    }

    void move_up(int n) const
    {
        if (n <= 0)
            return;
        else if (n > 3200)
            n = 3200;

        if (::is_ipython())
        {
            // we can only clear
            this->ipython_clear(true);
        }
        else
        {
            ::sys_stdout_write(QString("\x1b[%1A").arg(n).toUtf8().constData(), false);
        }
    }

    bool stdout_is_atty() const
    {
        return ::sys_stdout_is_atty();
    }

    void stdout_write(const QString &text, bool flush) const
    {
        QByteArray bytes = text.toUtf8();

        // need to split into chunks of 1000 characters
        // because python can only cope with 1000 characters at once
        while (bytes.length() > 1000)
        {
            auto start = bytes.left(999);
            start.append('\0');
            ::sys_stdout_write(start.constData(), false);
            bytes = bytes.right(bytes.length() - 999);
        }

        if (bytes.length() > 0)
            ::sys_stdout_write(bytes.constData(), flush);
    }

    SireBase::GILHandle releaseGIL() const
    {
        QMutexLocker lkr(&release_mutex);

        auto handle = current_state.lock();

        if (handle.get() != 0)
            return handle;

        auto thread_state = PyEval_SaveThread();
        handle = SireBase::GILHandle(new ReleaseGIL(thread_state));

        current_state = handle;

        return handle;
    }

private:
    static QMutex release_mutex;
    static std::weak_ptr<SireBase::detail::ReleaseGILBase> current_state;

    PyThreadState *thread_state;
};

class ConsoleImpl : public SireBase::ConsoleBase
{
public:
    ConsoleImpl() : SireBase::ConsoleBase(), pyconsole(0)
    {
    }

    ~ConsoleImpl()
    {
        if (pyconsole)
        {
            PyGILState_STATE gilstate = PyGILState_Ensure();
            Py_DECREF(pyconsole);
            delete pyconsole;
            pyconsole = 0;
            PyGILState_Release(gilstate);
        }
    }

    void debug(const QString &message) const
    {
        const_cast<ConsoleImpl *>(this)->get_pyconsole();

        if (pyconsole)
        {
            PyGILState_STATE gilstate = PyGILState_Ensure();

            PyObject *result = PyObject_CallMethod(pyconsole, "debug", "s", message.toUtf8().constData());

            if (result == 0)
            {
                qDebug() << "UNABLE TO CALL DEBUG";
                qDebug() << message;
            }
            else
            {
                Py_DECREF(result);
            }

            PyGILState_Release(gilstate);
        }
        else
        {
            qDebug() << "NO PYCONSOLE";
            qDebug() << message;
        }
    }

    void warning(const QString &message) const
    {
        const_cast<ConsoleImpl *>(this)->get_pyconsole();

        if (pyconsole)
        {
            PyGILState_STATE gilstate = PyGILState_Ensure();

            PyObject *result = PyObject_CallMethod(pyconsole, "warning", "s", message.toUtf8().constData());

            if (result == 0)
            {
                qWarning() << message;
            }
            else
            {
                Py_DECREF(result);
            }

            PyGILState_Release(gilstate);
        }
        else
        {
            qWarning() << message;
        }
    }

    void error(const QString &message) const
    {
        const_cast<ConsoleImpl *>(this)->get_pyconsole();

        if (pyconsole)
        {
            PyGILState_STATE gilstate = PyGILState_Ensure();

            PyObject *result = PyObject_CallMethod(pyconsole, "error", "s", message.toUtf8().constData());

            if (result == 0)
            {
                qCritical() << message;
            }
            else
            {
                Py_DECREF(result);
            }

            PyGILState_Release(gilstate);
        }
        else
        {
            qCritical() << message;
        }
    }

    void info(const QString &message) const
    {
        const_cast<ConsoleImpl *>(this)->get_pyconsole();

        if (pyconsole)
        {
            PyGILState_STATE gilstate = PyGILState_Ensure();

            PyObject *result = PyObject_CallMethod(pyconsole, "info", "s", message.toUtf8().constData());

            if (result == 0)
            {
                qInfo() << "UNABLE TO CALL INFO";
                qInfo() << message;
            }
            else
            {
                Py_DECREF(result);
            }

            PyGILState_Release(gilstate);
        }
        else
        {
            qInfo() << "NO PYCONSOLE";
            qInfo() << message;
        }
    }

private:
    void get_pyconsole()
    {
        if (pyconsole)
            return;

        PyGILState_STATE gilstate = PyGILState_Ensure();

        PyObject *sire_utils = PyImport_ImportModule("sire.utils");

        if (sire_utils == 0)
        {
            qWarning() << "COULD NOT IMPORT SIRE.UTILS";
            PyGILState_Release(gilstate);
            return;
        }

        pyconsole = PyObject_GetAttrString(sire_utils, "Console");

        Py_DECREF(sire_utils);

        if (pyconsole == 0)
        {
            qWarning() << "Could not import Console from sire.utils";
            PyGILState_Release(gilstate);
            return;
        }
    }

    PyObject *pyconsole;
};

QMutex ReleaseGIL::release_mutex;

std::weak_ptr<SireBase::detail::ReleaseGILBase> ReleaseGIL::current_state;

void register_releasegil()
{
    ReleaseGIL::register_releasegil();

    boost::python::def("set_is_ipython", &set_is_ipython);

    SireBase::Console::setConsole(new ConsoleImpl());
}
