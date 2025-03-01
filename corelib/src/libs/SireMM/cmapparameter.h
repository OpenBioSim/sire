/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2025  Christopher Woods
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

#ifndef SIREMM_CMAPPARAMETER_H
#define SIREMM_CMAPPARAMETER_H

#include "SireBase/array2d.hpp"

namespace SireMM
{
    class CMAPParameter;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::CMAPParameter &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::CMAPParameter &);

namespace SireMM
{
    SIREMM_EXPORT uint qHash(const CMAPParameter &param);

    class SIREMM_EXPORT CMAPParameter
    {

        friend QDataStream & ::operator<<(QDataStream &, const CMAPParameter &);
        friend QDataStream & ::operator>>(QDataStream &, CMAPParameter &);

        friend uint qHash(const CMAPParameter &param);

    public:
        CMAPParameter();
        CMAPParameter(const SireBase::Array2D<double> &grid);
        CMAPParameter(const CMAPParameter &other);
        ~CMAPParameter();

        CMAPParameter &operator=(const CMAPParameter &other);

        bool operator==(const CMAPParameter &other) const;
        bool operator!=(const CMAPParameter &other) const;

        bool operator<(const CMAPParameter &other) const;
        bool operator<=(const CMAPParameter &other) const;
        bool operator>(const CMAPParameter &other) const;
        bool operator>=(const CMAPParameter &other) const;

        static const char *typeName();
        const char *what() const;

        CMAPParameter *clone() const;

        QString toString() const;

        const SireBase::Array2D<double> &grid() const;

        bool isEmpty() const;
        bool isNull() const;

        QVector<double> values() const;

        int nRows() const;
        int nColumns() const;

    private:
        SireBase::Array2D<double> param;
    };

} // namespace SireMM

Q_DECLARE_METATYPE(SireMM::CMAPParameter);

SIRE_EXPOSE_CLASS(SireMM::CMAPParameter)

SIRE_END_HEADER

#endif