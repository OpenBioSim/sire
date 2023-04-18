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

#include "third_party/indicators.hpp" // CONDITIONAL_INCLUDE

#include <QDateTime>

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<ProgressBar> r_bar;

QDataStream &operator<<(QDataStream &ds, const ProgressBar &bar)
{
    writeHeader(ds, r_bar, 1);

    SharedDataStream sds(ds);

    sds << bar.total_value << bar.progress_text << bar.show_time;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, ProgressBar &bar)
{
    VersionID v = readHeader(ds, r_bar);

    if (v == 1)
    {
        bar = ProgressBar();

        SharedDataStream sds(ds);
        sds >> bar.total_value >> bar.progress_text >> bar.show_time;
    }
    else
        throw version_error(v, "1", r_bar, CODELOC);

    return ds;
}

ProgressBar::ProgressBar()
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(0), start_ms(0), show_time(true)
{
}

ProgressBar::ProgressBar(qint64 total, bool time)
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(total), start_ms(0), show_time(time)
{
    if (total_value < 0)
        total_value = 0;
}

ProgressBar::ProgressBar(qint64 total, const QString &text, bool time)
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(total), start_ms(0), show_time(time)
{
    if (total_value < 0)
        total_value = 0;

    progress_text = text.toUtf8();
}

ProgressBar::ProgressBar(const QString &text)
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(0), start_ms(0), show_time(false)
{
    progress_text = text.toUtf8();
}

ProgressBar::ProgressBar(const ProgressBar &other)
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(other.total_value),
      progress_text(other.progress_text),
      progress_ptr(other.progress_ptr),
      spinner_ptr(other.spinner_ptr),
      start_ms(other.start_ms),
      show_time(other.show_time)
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
        total_value = other.total_value;
        progress_text = other.progress_text;
        progress_ptr = other.progress_ptr;
        spinner_ptr = other.spinner_ptr;
        start_ms = other.start_ms;
        show_time = other.show_time;
    }

    return *this;
}

bool ProgressBar::operator==(const ProgressBar &other) const
{
    return total_value == other.total_value and
           progress_text == other.progress_text and
           show_time == other.show_time and
           start_ms == other.start_ms;
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

bool ProgressBar::isActive() const
{
    return this->isProgress() or this->isSpinner();
}

bool ProgressBar::isProgress() const
{
    return progress_ptr.get() != 0;
}

bool ProgressBar::isSpinner() const
{
    return spinner_ptr.get() != 0;
}

void ProgressBar::tick()
{
    if (this->isProgress())
    {
        auto ms = QDateTime::currentMSecsSinceEpoch();

        if (ms < start_ms)
            return;

        progress_ptr->tick();

        // make we don't do more than 10 updates per second
        start_ms = ms + 100;
    }
    else if (this->isSpinner())
    {
        auto ms = QDateTime::currentMSecsSinceEpoch();

        if (ms < start_ms)
            return;

        spinner_ptr->tick();

        // make we don't do more than 10 updates per second
        start_ms = ms + 100;
    }
}

void ProgressBar::setProgress(qint64 value)
{
    if (this->isProgress())
    {
        auto ms = QDateTime::currentMSecsSinceEpoch();

        if (value >= total_value)
        {
            progress_ptr->set_progress(100);
            progress_ptr->mark_as_completed();
            start_ms = ms + 3600000;
            return;
        }

        if (ms < start_ms)
            return;

        int percent = (100 * value) / total_value;

        if (value < 0)
            percent = 0;

        progress_ptr->set_progress(percent);

        // make we don't do more than 10 updates per second
        start_ms = ms + 100;
    }
    else if (this->isSpinner())
    {
        auto ms = QDateTime::currentMSecsSinceEpoch();

        if (ms < start_ms)
            return;

        spinner_ptr->tick();

        // make we don't do more than 10 updates per second
        start_ms = ms + 100;
    }
}

void ProgressBar::setText(const QString &text)
{
    progress_text = text.toUtf8();

    if (this->isActive())
    {
        if (this->isSpinner())
            spinner_ptr->set_option(indicators::option::PostfixText{progress_text.constData()});

        else if (this->isProgress())
            progress_ptr->set_option(indicators::option::PostfixText{progress_text.constData()});
    }
}

bool ProgressBar::isDeterministic() const
{
    return total_value > 0;
}

static int output_fileno = STDOUT_FILENO;

ProgressBar ProgressBar::enter() const
{
    ProgressBar ret(*this);

    if (this->isDeterministic())
    {
        ret.progress_ptr.reset(new indicators::BlockProgressBar(
            indicators::option::PostfixText{this->progress_text.constData()},
            indicators::option::ForegroundColor{indicators::Color::white},
            indicators::option::ShowElapsedTime{show_time},
            indicators::option::ShowRemainingTime{show_time},
            indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}));
    }
    else
    {
        ret.spinner_ptr.reset(new indicators::ProgressSpinner(
            indicators::option::PostfixText{this->progress_text.constData()},
            indicators::option::ForegroundColor{indicators::Color::yellow},
            indicators::option::SpinnerStates{std::vector<std::string>{"⠈", "⠐", "⠠", "⢀", "⡀", "⠄", "⠂", "⠁"}},
            indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}));
    }

    ret.start_ms = QDateTime::currentMSecsSinceEpoch();

    return ret;
}

void ProgressBar::exit()
{
    if (this->isProgress())
    {
        progress_ptr.reset();
    }
    else if (this->isSpinner())
    {
        spinner_ptr.reset();
    }

    start_ms = 0;
}

void ProgressBar::setTheme(const QString &theme)
{
}

void ProgressBar::setSilent()
{
}

void ProgressBar::set_fileno(int fileno)
{
    output_fileno = fileno;
}
