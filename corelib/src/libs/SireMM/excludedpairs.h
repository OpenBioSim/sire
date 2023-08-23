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

#ifndef SIREMM_EXCLUDEDPAIRS_H
#define SIREMM_EXCLUDEDPAIRS_H

#include "SireMol/connectivity.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class ExcludedPairs;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::ExcludedPairs &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::ExcludedPairs &);

namespace SireMM
{

    /** This class holds the excluded atoms pair list for a molecule.
     *  Excluded atom pairs means that the non-bonded coulomb and
     *  LJ energy between these pair of atoms are not calculated.
     */
    class SIREMM_EXPORT ExcludedPairs : public SireBase::ConcreteProperty<ExcludedPairs, SireMol::MolViewProperty>
    {

        friend QDataStream & ::operator<<(QDataStream &, const ExcludedPairs &);
        friend QDataStream & ::operator>>(QDataStream &, ExcludedPairs &);

    public:
        ExcludedPairs();
        ExcludedPairs(const SireMol::MoleculeView &molecule,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());
        ExcludedPairs(const ExcludedPairs &other);

        ~ExcludedPairs();

        static const char *typeName();

        const char *what() const;

        ExcludedPairs &operator=(const ExcludedPairs &other);

        bool operator==(const ExcludedPairs &other) const;
        bool operator!=(const ExcludedPairs &other) const;

        int count() const;

        std::tuple<SireMol::AtomIdx, SireMol::AtomIdx> operator[](int i) const;

        ExcludedPairs *clone() const;

        QString toString() const;

        SireMol::MoleculeInfo info() const;

        bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

        void updateBondMatrix(QVector<QVector<bool>> &bond_matrix) const;

        int nExcludedPairs() const;

        bool areExcluded(const SireMol::AtomID &atom0, const SireMol::AtomID &atom1) const;

        void setExcluded(const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
                         bool are_excluded);

    private:
        int getIndex(qint64 atom0, qint64 atom1) const;

        /** The info object that describes the molecule */
        SireMol::MoleculeInfo minfo;

        /** All of the excluded pairs (atom indicies) */
        QVector<qint64> excl_pairs;
    };
}

Q_DECLARE_METATYPE(SireMM::ExcludedPairs)

SIRE_EXPOSE_CLASS(SireMM::ExcludedPairs)

SIRE_END_HEADER

#endif
