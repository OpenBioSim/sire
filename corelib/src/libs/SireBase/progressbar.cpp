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

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

using namespace indicators::option;
using typename indicators::Color;
using typename indicators::FontStyle;
using typename indicators::ProgressBarPtr;
using typename indicators::SpinnerPtr;

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
      total_value(0), start_ms(0), show_time(true),
      has_displayed(false), has_completed(true)
{
}

ProgressBar::ProgressBar(qint64 total, bool time)
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(total), start_ms(0), show_time(time),
      has_displayed(false), has_completed(true)
{
    if (total_value < 0)
        total_value = 0;
}

ProgressBar::ProgressBar(qint64 total, const QString &text, bool time)
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(total), start_ms(0), show_time(time),
      has_displayed(false), has_completed(true)
{
    if (total_value < 0)
        total_value = 0;

    progress_text = text.toUtf8();
}

ProgressBar::ProgressBar(const QString &text)
    : ConcreteProperty<ProgressBar, Property>(),
      total_value(0), start_ms(0), show_time(false),
      has_displayed(false), has_completed(true)
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
      show_time(other.show_time),
      has_displayed(other.has_displayed),
      has_completed(other.has_completed)
{
}

ProgressBar::~ProgressBar()
{
    this->exit();
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
        has_displayed = other.has_displayed;
        has_completed = other.has_completed;
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

// minimum amount of time to wait between updates
static const int CONSOLE_TICK_DELAY_TIME = 100;
static const int JUPYTER_TICK_DELAY_TIME = 250;

static int TICK_DELAY_TIME = CONSOLE_TICK_DELAY_TIME;

void ProgressBar::tick()
{
    this->tick(QString());
}

void ProgressBar::tick(const QString &text)
{
    const bool update_text = not text.isNull();

    if (this->isProgress())
    {
        if (has_completed)
            return;

        auto ms = QDateTime::currentMSecsSinceEpoch();

        if ((not update_text) and ms < start_ms)
            return;

        if (update_text)
            progress_ptr->set_option(PostfixText{text.toUtf8().constData()});

        if (not has_displayed)
        {
            indicators::show_console_cursor(false);
            has_displayed = true;
        }

        progress_ptr->tick();

        start_ms = ms + TICK_DELAY_TIME;
    }
    else if (this->isSpinner())
    {
        if (has_completed)
            return;

        auto ms = QDateTime::currentMSecsSinceEpoch();

        if ((not update_text) and ms < start_ms)
            return;

        if (update_text)
            spinner_ptr->set_option(PostfixText{text.toUtf8().constData()});

        if (not has_displayed)
        {
            indicators::show_console_cursor(false);
            has_displayed = true;
        }

        spinner_ptr->tick();

        start_ms = ms + TICK_DELAY_TIME;
    }
}

void ProgressBar::setCompleted()
{
    // qDebug() << "setCompleted" << this << has_displayed << has_completed << this->isProgress() << this->isSpinner();

    if (has_completed)
        return;

    if (not has_displayed)
    {
        // no need to show anything as nothing has yet been displayed
        has_completed = true;
        return;
    }

    if (this->isProgress())
    {
        // qDebug() << this << "clear_progress";
        //  calling function I added to actually clear the bar
        //  progress_ptr->clear_at_end();
    }

    if (this->isSpinner())
    {
        // qDebug() << this << "clear_spinner";
        //  calling function I added to actually clear the spinner
        spinner_ptr->reset_at_end();
    }

    has_completed = true;
}

void ProgressBar::setProgress(qint64 value)
{
    this->setProgress(value, QString());
}

void ProgressBar::setProgress(qint64 value, const QString &text)
{
    const bool update_text = not text.isNull();

    if (this->isProgress())
    {
        if (has_completed)
            return;

        auto ms = QDateTime::currentMSecsSinceEpoch();

        if (value >= total_value)
        {
            if (not has_displayed)
            {
                // we will never display
                has_completed = true;
                return;
            }

            if (update_text)
                progress_ptr->set_option(PostfixText{text.toUtf8().constData()});

            progress_ptr->set_progress(100);
            return;
        }

        if ((not update_text) and ms < start_ms)
            return;

        int percent = (100 * value) / total_value;

        if (value < 0)
            percent = 0;

        if (update_text)
            progress_ptr->set_option(PostfixText{text.toUtf8().constData()});

        if (not has_displayed)
        {
            indicators::show_console_cursor(false);
            has_displayed = true;
        }

        progress_ptr->set_progress(percent);

        start_ms = ms + TICK_DELAY_TIME;
    }
    else if (this->isSpinner())
    {
        if (has_completed)
            return;

        auto ms = QDateTime::currentMSecsSinceEpoch();

        if ((not update_text) and ms < start_ms)
            return;

        if (update_text)
            spinner_ptr->set_option(PostfixText{text.toUtf8().constData()});

        if (not has_displayed)
        {
            indicators::show_console_cursor(false);
            has_displayed = true;
        }

        spinner_ptr->tick();

        start_ms = ms + TICK_DELAY_TIME;
    }
}

const char *ProgressBar::text() const
{
    return progress_text.constData();
}

bool ProgressBar::showTime() const
{
    return show_time;
}

bool ProgressBar::isDeterministic() const
{
    return total_value > 0;
}

int ProgressBar::barSize() const
{
    return std::max(20, 50 - progress_text.count());
}

void jupyter_theme(const ProgressBar &p, ProgressBarPtr &bar, SpinnerPtr &spinner)
{
    if (p.isDeterministic())
    {
        bar.reset(new indicators::ProgressBar(
            PostfixText{p.text()},
            BarWidth{p.barSize()},
            Start{"■"},
            Fill{"■"},
            Lead{"▶"},
            Remainder{" "},
            End{" "},
            ShowElapsedTime{p.showTime()},
            ShowRemainingTime{p.showTime()}));
    }
    else
    {
        spinner.reset(new indicators::IndeterminateProgressBar(
            PostfixText{p.text()},
            BarWidth{10},
            Start{" "},
            Fill{" "},
            Lead{"◀■▶"},
            End{" "}));
    }
}

void color_theme(const ProgressBar &p, ProgressBarPtr &bar, SpinnerPtr &spinner)
{
    if (p.isDeterministic())
    {
        bar.reset(new indicators::ProgressBar(
            PostfixText{p.text()},
            BarWidth{p.barSize()},
            Start{"■"},
            Fill{"■"},
            Lead{"▶"},
            Remainder{" "},
            End{" "},
            ForegroundColor{Color::green},
            ShowElapsedTime{p.showTime()},
            ShowRemainingTime{p.showTime()}));
    }
    else
    {
        spinner.reset(new indicators::IndeterminateProgressBar(
            PostfixText{p.text()},
            BarWidth{10},
            Start{" "},
            Fill{" "},
            Lead{"◀■▶"},
            End{" "},
            ForegroundColor{Color::green},
            FontStyles{std::vector<FontStyle>{FontStyle::bold}}));
    }
}

void simple_theme(const ProgressBar &p, ProgressBarPtr &bar, SpinnerPtr &spinner)
{
    if (p.isDeterministic())
    {
        bar.reset(new indicators::ProgressBar(
            PostfixText{p.text()},
            BarWidth{p.barSize()},
            Start{"["},
            Fill{"="},
            Lead{">"},
            Remainder{" "},
            End{"]"},
            ForegroundColor{Color::white},
            ShowElapsedTime{p.showTime()},
            ShowRemainingTime{p.showTime()},
            FontStyles{std::vector<FontStyle>{FontStyle::bold}}));
    }
    else
    {
        spinner.reset(new indicators::IndeterminateProgressBar(
            PostfixText{p.text()},
            BarWidth{10},
            Start{"["},
            Fill{" "},
            Lead{"<=>"},
            End{"]"},
            ForegroundColor{Color::white},
            FontStyles{std::vector<FontStyle>{FontStyle::bold}}));
    }
}

ProgressBar ProgressBar::enter() const
{
    ProgressBar ret(*this);
    ret.has_displayed = false;
    ret.has_completed = false;

    if (current_theme == 1)
        simple_theme(ret, ret.progress_ptr, ret.spinner_ptr);
    else if (current_theme == 2)
        color_theme(ret, ret.progress_ptr, ret.spinner_ptr);
    else if (current_theme == 3)
        jupyter_theme(ret, ret.progress_ptr, ret.spinner_ptr);

    // start this after a delay of 500 ms
    ret.start_ms = QDateTime::currentMSecsSinceEpoch() + 500;

    return ret;
}

void ProgressBar::exit()
{
    // qDebug() << "exit" << this << has_displayed << has_completed << this->isProgress() << this->isSpinner();

    this->setCompleted();

    progress_ptr.reset();
    spinner_ptr.reset();

    if (has_displayed)
    {
        // qDebug() << this << "reshow_console_cursor";
        indicators::show_console_cursor(true);
        has_displayed = false;
    }

    start_ms = 0;
}

// Default to the 'color' theme
int ProgressBar::current_theme(2);

void ProgressBar::setTheme(QString theme)
{
    theme = theme.toLower().simplified();

    TICK_DELAY_TIME = CONSOLE_TICK_DELAY_TIME;

    if (theme == "simple")
        current_theme = 1;
    else if (theme == "color" or theme == "colour")
        current_theme = 2;
    else if (theme == "silent" or theme == "quiet")
        current_theme = 0;
    else if (theme == "jupyter")
    {
        current_theme = 3;
        TICK_DELAY_TIME = JUPYTER_TICK_DELAY_TIME;
    }
    else
    {
        throw SireError::invalid_key(QObject::tr(
                                         "Unrecognised progress bar theme '%1'. Recognised themes are "
                                         "'silent', 'simple' or 'color'")
                                         .arg(theme),
                                     CODELOC);
    }
}

void ProgressBar::setSilent()
{
    current_theme = 0;
}
