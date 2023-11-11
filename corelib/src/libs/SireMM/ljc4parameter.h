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

#ifndef SIREMM_LJC4PARAMETER_H
#define SIREMM_LJC4PARAMETER_H

#include "ljparameter.h"

#include "SireUnits/generalunit.h"

#include <QHash>

#include <memory>

SIRE_BEGIN_HEADER

namespace SireMM
{
    class LJC4Parameter;
}

QDataStream &operator<<(QDataStream &, const SireMM::LJC4Parameter &);
QDataStream &operator>>(QDataStream &, SireMM::LJC4Parameter &);

namespace SireMM
{
    /** This class holds the C4 parameter for a pair of LJParameters.
     *
     *  Each LJParameter that has a corresponding C4 parameter will
     *  be given a unique ID. This is used to look up the C4 parameter
     *  for the pair.
     */
    class SIREMM_EXPORT LJC4Parameter
    {
        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const LJC4Parameter &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, LJC4Parameter &);

    public:
        LJC4Parameter();
        LJC4Parameter(const LJC4Parameter &other);
        ~LJC4Parameter();

        LJC4Parameter &operator=(const LJC4Parameter &other);

        bool operator==(const LJC4Parameter &other) const;
        bool operator!=(const LJC4Parameter &other) const;

        void addParameter(LJParameter &lj0, const LJParameter &lj1,
                          SireUnits::Dimension::GeneralUnit c4);

        void getParameter(const LJParameter &lj0, const LJParameter &lj1);

        std::shared_ptr<LJC4Parameter> self() const;

    private:
        std::weak_ptr<LJC4Parameter> m_self;

        QVector<QVector<double>> params;
    }
}

SIRE_END_HEADER

#endif
