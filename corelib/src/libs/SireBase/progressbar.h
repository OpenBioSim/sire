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
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREBASE_PROGRESSBAR_H
#define SIREBASE_PROGRESSBAR_H

#include "property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
    class ProgressBar;
}

SIREBASE_EXPORT QDataStream &operator<<(QDataStream &, const SireBase::ProgressBar &);
SIREBASE_EXPORT QDataStream &operator>>(QDataStream &, SireBase::ProgressBar &);

namespace SireBase
{
    namespace detail
    {
        class BarManager;
        class BarData;
    }

    /** This is a progress bar */
    class SIREBASE_EXPORT ProgressBar : public SireBase::ConcreteProperty<ProgressBar, Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const ProgressBar &);
        friend QDataStream & ::operator>>(QDataStream &, ProgressBar &);

        friend class detail::BarManager;

    public:
        ProgressBar();
        ProgressBar(quint32 total);
        ProgressBar(quint32 total, const QString &text);
        ProgressBar(const QString &text, quint32 total);

        ProgressBar(const QString &text);

        ProgressBar(const ProgressBar &other);

        ~ProgressBar();

        ProgressBar *clone() const;

        ProgressBar &operator=(const ProgressBar &other);

        bool operator==(const ProgressBar &other) const;
        bool operator!=(const ProgressBar &other) const;

        static const char *typeName();
        const char *what() const;

        QString message() const;
        QString speedUnit() const;

        void setSpeedUnit(const QString &unit);

        void tick();
        void tick(const QString &text);
        void silentTick();

        void setProgress(quint32 value);
        void setProgress(quint32 value, const QString &text);
        void setProgress(const QString &text, quint32 value);

        void success();
        void success(const QString &message);

        void failure();
        void failure(const QString &message);

        quint32 total() const;
        quint32 current() const;

        ProgressBar enter() const;
        void exit();

        static void setTheme(QString theme);
        static void setSilent();

    private:
        std::shared_ptr<detail::BarData> d;
        QString msg;
        QString speed_unit;
    };
}

Q_DECLARE_METATYPE(SireBase::ProgressBar)

SIRE_EXPOSE_CLASS(SireBase::ProgressBar)

SIRE_END_HEADER

#endif
