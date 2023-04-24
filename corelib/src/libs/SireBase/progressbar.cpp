/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2023  Christopher Woods
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
  *  You can contact the authors via the website
  *  at http://sire.openbiosim.org
  *
\*********************************************/

#include "progressbar.h"

#include <QDateTime>
#include <QThread>
#include <QAtomicInteger>

#include <boost/noncopyable.hpp>

#include "SireBase/releasegil.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;

using namespace SireStream;

/////////
///////// Functions copied from
/////////

#if defined(_MSC_VER)
#if !defined(NOMINMAX)
#define NOMINMAX
#endif
#include <io.h>
#include <windows.h>
#else
#include <iostream>
#endif

#ifdef _MSC_VER

static inline void move(int x, int y)
{
    auto hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    if (!hStdout)
        return;

    CONSOLE_SCREEN_BUFFER_INFO csbiInfo;
    GetConsoleScreenBufferInfo(hStdout, &csbiInfo);

    COORD cursor;

    cursor.X = csbiInfo.dwCursorPosition.X + x;
    cursor.Y = csbiInfo.dwCursorPosition.Y + y;
    SetConsoleCursorPosition(hStdout, cursor);
}

static inline void move_up(int lines) { move(0, -lines); }
static inline void move_down(int lines) { move(0, -lines); }
static inline void move_right(int cols) { move(cols, 0); }
static inline void move_left(int cols) { move(-cols, 0); }

#else

static inline void move_up(int lines)
{
    std::cout << "\033[" << lines << "A";
}
static inline void move_down(int lines) { std::cout << "\033[" << lines << "B"; }
static inline void move_right(int cols) { std::cout << "\033[" << cols << "C"; }
static inline void move_left(int cols) { std::cout << "\033[" << cols << "D"; }

#endif

/////////
///////// Implementation of SireBase::detail::BarData
/////////

namespace SireBase
{
    namespace detail
    {
        class BarManager;

        class BarData : boost::noncopyable
        {
        public:
            BarData()
                : boost::noncopyable(), total(0), current(0), failed(false)
            {
            }

            BarData(quint32 t)
                : boost::noncopyable(), total(t), current(0), failed(false)
            {
            }

            ~BarData()
            {
            }

            std::shared_ptr<BarData> clone() const
            {
                std::shared_ptr<BarData> ret(new BarData());

                ret->total = total;
                ret->current = current;
                ret->failed = failed;

                return ret;
            }

            void setSuccess();
            void setFailure();

            void tick();
            void tick(const QString &text);

            void setProgress(quint32 c);
            void setProgress(quint32 c, const QString &text);

            std::tuple<QString, bool> toString(qint64 elapsed_ms,
                                               const QString &message,
                                               const QString &speed_unit) const;

            void exit();

            void success();

            void success(const QString &message);

            void failure();

            void failure(const QString &message);

            void setSpeedUnit(const QString &unit);

            std::shared_ptr<BarManager> manager;

            QAtomicInteger<quint32> total;
            QAtomicInteger<quint32> current;
            bool failed;
        };

        /** Internal struct of data */
        struct Bar
        {
            std::weak_ptr<BarData> d;
            QString message;
            QString speed_unit;
            qint64 start_time;
            QString last_string;
        };

        void manager_event_loop();
        void print_bars(const QStringList &bars, int n_move_up);

        /** Internal class that manages all of the progress bars that are visible.
         *  The design is that the progress bars add themselves to a manager, and
         *  the manager takes responsibility for printing them to the screen
         *  in a background thread. This way, updating of the progress bar
         *  shouldn't slow down the main application, plus all progress bars
         *  from multiple threads won't trip over each other
         */
        class BarManager : boost::noncopyable
        {
        public:
            BarManager() : n_move_up(0), has_finished(false)
            {
                start_time = QDateTime::currentMSecsSinceEpoch();
            }

            ~BarManager()
            {
                if (n_move_up <= 0)
                {
                    // we haven't printed anything?
                    return;
                }

                int attempt = 0;

                while (true)
                {
                    if (mutex.tryLock(10))
                    {
                        break;
                    }
                    else if (attempt > 100)
                    {
                        return;
                    }

                    attempt += 1;
                }

                // must print the final bars, so that we
                // get the 100% and completed signal
                QStringList bar_strings;
                try
                {
                    bar_strings = this->_lkr_printBars();
                    mutex.unlock();
                }
                catch (...)
                {
                    mutex.unlock();
                }

                if (not bar_strings.isEmpty())
                {
                    print_bars(bar_strings, n_move_up);
                }
            }

            static std::shared_ptr<BarManager> getManager()
            {
                QMutexLocker lkr(&mutex);

                auto manager = global_manager.lock();

                if (manager.get() == 0)
                {
                    manager.reset(new BarManager());

                    // we only need to create a new thread
                    // if one doesn't already exist
                    if (global_thread != 0)
                    {
                        global_thread->requestInterruption();

                        while (not global_thread->wait(1000))
                        {
                            if (global_thread->isFinished())
                                break;
                        }

                        delete global_thread;
                        global_thread = 0;
                    }

                    global_manager = manager;

                    global_thread = QThread::create(&manager_event_loop);
                    global_thread->start();
                }

                return manager;
            }

            static void eventLoop()
            {
                // don't print anything for the first 250 ms
                for (int i = 0; i < 5; ++i)
                {
                    QThread::msleep(50);

                    if (QThread::currentThread()->isInterruptionRequested())
                    {
                        return;
                    }
                }

                while (true)
                {
                    if (QThread::currentThread()->isInterruptionRequested())
                    {
                        break;
                    }

                    int attempt = 0;

                    while (true)
                    {
                        if (mutex.tryLock(10))
                        {
                            break;
                        }
                        else if (QThread::currentThread()->isInterruptionRequested())
                        {
                            return;
                        }
                        else if (attempt > 200)
                        {
                            qWarning() << "FAILED TO GET A LOCK";
                            return;
                        }

                        attempt += 1;
                    }

                    auto manager = global_manager.lock();

                    if (manager.get() == 0)
                    {
                        mutex.unlock();
                        break;
                    }
                    else if (manager->has_finished)
                    {
                        mutex.unlock();
                        break;
                    }

                    QStringList bar_strings;

                    try
                    {
                        bar_strings = manager->_lkr_printBars();
                    }
                    catch (...)
                    {
                        mutex.unlock();
                        throw;
                    }

                    int n_move_up = manager->n_move_up;
                    manager->n_move_up = bar_strings.count();

                    mutex.unlock();

                    // printing unlocked so we don't deadlock on the GIL
                    print_bars(bar_strings, n_move_up);

                    manager.reset();

                    if (QThread::currentThread()->isInterruptionRequested())
                    {
                        break;
                    }

                    QThread::msleep(sleep_time);
                }
            }

            void add(ProgressBar &bar)
            {
                QMutexLocker lkr(&mutex);

                if (not bar.d.unique())
                    bar.d = bar.d->clone();

                if (has_finished)
                    return;

                bar.d->manager = global_manager.lock();

                Bar b;
                b.d = bar.d;
                b.message = bar.message();
                b.speed_unit = bar.speedUnit();
                b.start_time = QDateTime::currentMSecsSinceEpoch();

                bars.append(b);
            }

            void setMessage(const BarData *d, const QString &message)
            {
                QMutexLocker lkr(&mutex);

                if (has_finished)
                    return;

                for (auto &bar : this->bars)
                {
                    if (bar.d.lock().get() == d)
                    {
                        bar.message = message;
                        return;
                    }
                }
            }

            void setSpeedUnit(const BarData *d, const QString &unit)
            {
                QMutexLocker lkr(&mutex);

                if (has_finished)
                    return;

                for (auto &bar : this->bars)
                {
                    if (bar.d.lock().get() == d)
                    {
                        bar.speed_unit = unit;
                        return;
                    }
                }
            }

            void reprint(const BarData *d)
            {
                QMutexLocker lkr(&mutex);

                if (has_finished)
                    return;

                for (auto &bar : this->bars)
                {
                    if (bar.d.lock().get() == d)
                    {
                        auto current_time = QDateTime::currentMSecsSinceEpoch();
                        auto result = d->toString(current_time - bar.start_time,
                                                  bar.message, bar.speed_unit);
                        bar.last_string = std::get<0>(result);
                        return;
                    }
                }
            }

        private:
            QStringList _lkr_printBars();
            QStringList _lkr_finish();

            static QMutex mutex;
            static std::weak_ptr<BarManager> global_manager;
            static quint64 sleep_time;

            QList<Bar> bars;

            static QThread *global_thread;

            qint64 start_time;

            int n_move_up;

            bool has_finished;
        };

        void manager_event_loop()
        {
            BarManager::eventLoop();
        }

        QMutex BarManager::mutex;

        QThread *BarManager::global_thread(0);

        std::weak_ptr<BarManager> BarManager::global_manager;

        quint64 BarManager::sleep_time(100);

        void BarData::exit()
        {
            if (total != current)
            {
                // we've exited without signalling success or failure
                // Let's assume success
                this->setSuccess();
            }

            manager->reprint(this);
            manager.reset();
        }

        void BarData::success()
        {
            this->setSuccess();
            this->exit();
        }

        void BarData::success(const QString &message)
        {
            if (manager.get() != 0)
            {
                manager->setMessage(this, message);
            }

            this->setSuccess();
            this->exit();
        }

        void BarData::failure()
        {
            this->setFailure();
            this->exit();
        }

        void BarData::failure(const QString &message)
        {
            if (manager.get() != 0)
            {
                manager->setMessage(this, message);
            }

            this->setFailure();
            this->exit();
        }

        void BarData::setSuccess()
        {
            total = current.loadAcquire();
        }

        void BarData::setFailure()
        {
            total = 0;
            current = 0;
        }

        void BarData::tick()
        {
            current += 1;
        }

        void BarData::tick(const QString &text)
        {
            this->tick();

            if (manager.get() != 0)
                manager->setMessage(this, text);
        }

        void BarData::setProgress(quint32 c)
        {
            current = c;
        }

        void BarData::setProgress(quint32 c, const QString &text)
        {
            this->setProgress(c);

            if (manager.get() != 0)
                manager->setMessage(this, text);
        }

    }
}

void SireBase::detail::BarData::setSpeedUnit(const QString &unit)
{
    if (manager.get() != 0)
        manager->setSpeedUnit(this, unit);
}

void ProgressBar::setTheme(QString theme)
{
}

void ProgressBar::setSilent()
{
}

std::tuple<QString, bool> SireBase::detail::BarData::toString(qint64 elapsed,
                                                              const QString &text,
                                                              const QString &speed_unit) const
{
    static int frame_counter = 0;

    QString bar;
    bool finished = false;

    quint32 c = this->current;
    quint32 t = this->total;

    const float secs = 0.001 * elapsed;
    QString speed;

    if (secs > 0)
    {
        if (speed_unit.length() == 0)
            speed = QString("%1 its / s").arg(float(c) / secs, 0, 'F', 1);
        else
            speed = QString("%1 %2").arg(float(c) / secs, 0, 'F', 1).arg(speed_unit);
    }

    if (this->failed)
    {
        bar = QString("Failed : %1 s : %2")
                  .arg(secs, 0, 'F', 1)
                  .arg(speed);
        finished = true;
    }
    else if (t > 0)
    {
        // progress bar
        if (c >= t)
        {
            bar = QString("Completed : %1 s : %2")
                      .arg(secs, 0, 'F', 1)
                      .arg(speed);
            finished = true;
        }
        else
        {
            int percent = int((100.0 * c) / float(t));

            bar = QString("%1% : %2 s : %3")
                      .arg(percent, 2)
                      .arg(secs, 0, 'F', 1)
                      .arg(speed, 0);
        }
    }
    else
    {
        int frame = frame_counter % 9;
        frame_counter += 1;

        QString part = ".  ";

        if (frame >= 3 and frame < 6)
            part = " . ";
        else if (frame >= 6)
            part = "  .";

        bar = QString("%1 : %2 s : %3")
                  .arg(part)
                  .arg(secs, 0, 'F', 1)
                  .arg(speed, 0);
    }

    if (text.length() > 0)
    {
        bar = QString("%1 : %2").arg(text).arg(bar);
    }

    return std::make_tuple(bar, finished);
}

void SireBase::detail::print_bars(const QStringList &bars, int n_move_up)
{
    if (bars.isEmpty())
        return;

    // construct the bars - cap at 80 columns wide
    const int bar_width = 80;

    QString to_print;

    for (const auto &bar : bars)
    {
        if (bar.length() > bar_width)
        {
            to_print += bar.left(bar_width);
        }
        else if (bar.length() < bar_width)
        {
            QString b = bar;
            b.resize(bar_width, ' ');
            to_print += b;
        }
        else
        {
            to_print += bar;
        }
    }

    to_print += "\n";

    if (n_move_up > 0)
        move_up(n_move_up);

    std::cout << to_print.toUtf8().constData();
    std::cout.flush();
}

QStringList SireBase::detail::BarManager::_lkr_printBars()
{
    auto current_time = QDateTime::currentMSecsSinceEpoch();

    if (this->bars.isEmpty())
    {
        // need to wait for some bars to be added...
        auto delta = current_time - start_time;

        if (delta > 20000)
        {
            // we should have added some bars within 20 seconds...
            // We haven't so exit
            has_finished = true;
            n_move_up = 0;
            return QStringList();
        }
    }

    bool all_finished = true;

    QStringList bar_strings;

    for (auto &bar : bars)
    {
        auto b = bar.d.lock();

        if (b.get() != 0)
        {
            auto result = b->toString(current_time - bar.start_time,
                                      bar.message, bar.speed_unit);

            auto bar_string = std::get<0>(result);
            bar_strings.append(bar_string);
            bar.last_string = bar_string;

            if (not std::get<1>(result))
                all_finished = false;
        }
        else
        {
            bar_strings.append(bar.last_string);
        }
    }

    if (all_finished)
    {
        bar_strings += this->_lkr_finish();
    }

    return bar_strings;
}

QStringList SireBase::detail::BarManager::_lkr_finish()
{
    // print something to finish the bars...

    has_finished = true;
    return QStringList();
}

/////////
///////// Implementation of ProgressBar
/////////

static const RegisterMetaType<ProgressBar> r_bar;

QDataStream &operator<<(QDataStream &ds, const ProgressBar &bar)
{
    writeHeader(ds, r_bar, 1);

    SharedDataStream sds(ds);

    sds << bar.message() << bar.speedUnit() << bar.total() << bar.current();

    return ds;
}

QDataStream &operator>>(QDataStream &ds, ProgressBar &bar)
{
    VersionID v = readHeader(ds, r_bar);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        QString message;
        QString speed_unit;
        quint32 total;
        quint32 current;

        sds >> message >> speed_unit >> total >> current;

        bar = ProgressBar(total);
        bar.setSpeedUnit(speed_unit);

        if (message.length() > 0)
        {
            bar.setProgress(current, message);
        }
        else
        {
            bar.setProgress(current);
        }
    }
    else
        throw version_error(v, "1", r_bar, CODELOC);

    return ds;
}

ProgressBar::ProgressBar()
    : ConcreteProperty<ProgressBar, Property>()
{
    d = std::shared_ptr<detail::BarData>(
        new detail::BarData());
}

ProgressBar::ProgressBar(quint32 total)
    : ConcreteProperty<ProgressBar, Property>()
{
    d = std::shared_ptr<detail::BarData>(
        new detail::BarData(total));
}

ProgressBar::ProgressBar(quint32 total, const QString &text)
    : ConcreteProperty<ProgressBar, Property>()
{
    d = std::shared_ptr<detail::BarData>(
        new detail::BarData(total));

    msg = text;
}

ProgressBar::ProgressBar(const QString &text, quint32 total)
    : ConcreteProperty<ProgressBar, Property>()
{
    d = std::shared_ptr<detail::BarData>(
        new detail::BarData(total));

    msg = text;
}

ProgressBar::ProgressBar(const QString &text)
    : ConcreteProperty<ProgressBar, Property>()
{
    d = std::shared_ptr<detail::BarData>(
        new detail::BarData());

    msg = text;
}

ProgressBar::ProgressBar(const ProgressBar &other)
    : ConcreteProperty<ProgressBar, Property>(),
      d(other.d), msg(other.msg), speed_unit(other.speed_unit)
{
}

ProgressBar::~ProgressBar()
{
}

ProgressBar *ProgressBar::clone() const
{
    return new ProgressBar(*this);
}

ProgressBar &ProgressBar::operator=(const ProgressBar &other)
{
    if (this != &other)
    {
        d = other.d;
        msg = other.msg;
        speed_unit = other.speed_unit;
    }

    return *this;
}

bool ProgressBar::operator==(const ProgressBar &other) const
{
    return d.get() == other.d.get();
}

bool ProgressBar::operator!=(const ProgressBar &other) const
{
    return not ProgressBar::operator==(other);
}

const char *ProgressBar::typeName()
{
    return QMetaType::typeName(qMetaTypeId<ProgressBar>());
}

const char *ProgressBar::what() const
{
    return ProgressBar::typeName();
}

void ProgressBar::setSpeedUnit(const QString &unit)
{
    speed_unit = unit;

    if (d.get() != 0)
        d->setSpeedUnit(unit);
}

void ProgressBar::tick()
{
    if (d.get() != 0)
        d->tick();
}

void ProgressBar::tick(const QString &text)
{
    if (d.get() != 0)
        d->tick(text);

    msg = text;
}

void ProgressBar::success()
{
    if (d.get() != 0)
        d->success();
}

void ProgressBar::success(const QString &message)
{
    if (d.get() != 0)
        d->success(message);

    msg = message;
}

void ProgressBar::failure()
{
    if (d.get() != 0)
        d->failure();
}

void ProgressBar::failure(const QString &message)
{
    if (d.get() != 0)
        d->failure(message);

    msg = message;
}

void ProgressBar::setProgress(quint32 value)
{
    if (d.get() != 0)
        d->setProgress(value);
}

void ProgressBar::setProgress(quint32 value, const QString &text)
{
    if (d.get() != 0)
        d->setProgress(value, text);

    msg = text;
}

void ProgressBar::setProgress(const QString &text, quint32 value)
{
    if (d.get() != 0)
        d->setProgress(value, text);

    msg = text;
}

ProgressBar ProgressBar::enter() const
{
    ProgressBar ret(*this);

    detail::BarManager::getManager()->add(ret);

    return ret;
}

quint32 ProgressBar::total() const
{
    if (d.get() == 0)
    {
        return 0;
    }
    else
    {
        return d->total.loadAcquire();
    }
}

quint32 ProgressBar::current() const
{
    if (d.get() == 0)
    {
        return 0;
    }
    else
    {
        return d->current.loadAcquire();
    }
}

QString ProgressBar::speedUnit() const
{
    return speed_unit;
}

QString ProgressBar::message() const
{
    return msg;
}

void ProgressBar::exit()
{
    if (d.get() != 0)
        d->exit();

    d.reset();
}
