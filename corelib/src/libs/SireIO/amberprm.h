/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#ifndef SIREIO_AMBERPRM_H
#define SIREIO_AMBERPRM_H

#include "moleculeparser.h"

#include "SireBase/sparsematrix.hpp"

#include "SireMM/mmdetail.h"

#include "SireMM/atomljs.h"
#include "SireMM/ljparameter.h"
#include "SireMM/lj1264parameter.h"

#include "SireMaths/vector.h"

#include <QHash>
#include <QSet>

SIRE_BEGIN_HEADER

namespace SireIO
{
    class AmberPrm;
    class AmberRst7;
} // namespace SireIO

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::AmberPrm &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::AmberPrm &);

namespace SireMol
{
    class MolEditor;
    class MoleculeInfoData;
} // namespace SireMol

namespace SireMM
{
    class AmberParams;
}

namespace SireIO
{

    /** This class represents an Amber-format parameter file, currently
        supporting top files produced from Amber7 until Amber16

        The format of this file is described here;

        http://ambermd.org/formats.html

        (specifically the "PARM" parameter/topology file specification)

        @author Christopher Woods
    */
    class SIREIO_EXPORT AmberPrm : public SireBase::ConcreteProperty<AmberPrm, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const AmberPrm &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, AmberPrm &);

    public:
        enum FLAG_TYPE
        {
            UNKNOWN = 0,
            INTEGER = 1,
            FLOAT = 2,
            STRING = 3
        };

        AmberPrm();

        AmberPrm(const QString &filename, const PropertyMap &map = PropertyMap());
        AmberPrm(const QStringList &lines, const PropertyMap &map = PropertyMap());

        AmberPrm(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        AmberPrm(const AmberPrm &other);

        ~AmberPrm();

        AmberPrm &operator=(const AmberPrm &other);

        bool operator==(const AmberPrm &other) const;
        bool operator!=(const AmberPrm &other) const;

        static const char *typeName();

        const char *what() const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;

        QString formatName() const;
        QStringList formatSuffix() const;

        QString formatDescription() const;

        bool isTopology() const;

        static AmberPrm parse(const QString &filename, const PropertyMap &map = PropertyMap());

        SireMol::Molecule getMolecule(int i, const PropertyMap &map = PropertyMap()) const;

        SireMol::Molecule getMolecule(int i, const AmberRst7 &rst, const PropertyMap &map = PropertyMap()) const;

        QVector<QString> linesForFlag(const QString &flag) const;

        QStringList flags() const;

        FLAG_TYPE flagType(const QString &flag) const;

        QVector<qint64> intData(const QString &flag) const;
        QVector<double> floatData(const QString &flag) const;
        QVector<QString> stringData(const QString &flag) const;

        QString title() const;

        int nAtoms() const;

        int nTypes() const;

        int nBonds() const;
        int nBondsWithHydrogen() const;
        int nBondsNoHydrogen() const;

        int nAngles() const;
        int nAnglesWithHydrogen() const;
        int nAnglesNoHydrogen() const;

        int nDihedrals() const;
        int nDihedralsWithHydrogen() const;
        int nDihedralsNoHydrogen() const;

        int nExcluded() const;
        int nResidues() const;

        int nMolecules() const;

        int nAtoms(int molidx) const;

        void assertSane() const;

        SireMM::MMDetail forcefield() const;

        QStringList warnings() const;

    protected:
        SireSystem::System startSystem(const PropertyMap &map) const;

    private:
        void parse(const PropertyMap &map);

        void rebuildAfterReload();
        void rebuildLJParameters();
        void rebuildBADIndicies();
        void rebuildExcludedAtoms();
        void rebuildMolNumToAtomNums();

        SireMM::AmberParams getAmberParams(int imol, const SireMol::MoleculeInfoData &molinfo) const;

        SireMol::MolEditor getMolStructure(int molidx, const SireBase::PropertyName &cutting) const;

        SireMol::MolEditor getMoleculeEditor(int molidx, const PropertyMap &map) const;

        QVector<int> getAtomIndexToMolIndex() const;

        /** Function to process all flags, returning the parsing score */
        double processAllFlags();

        /** A map showing the line number of all flags. This holds
            the start index and number of lines for each flag */
        QHash<QString, QPair<qint64, qint64>> flag_to_line;

        /** The raw int data for the integer flags */
        QHash<QString, QVector<qint64>> int_data;

        /** The raw float data for the float flags */
        QHash<QString, QVector<double>> float_data;

        /** The raw string data for the string flags */
        QHash<QString, QVector<QString>> string_data;

        /** All of the LJ parameters, indexed by atom type */
        QVector<SireMM::LJParameter> lj_data;

        /** All of the LJ exceptions, indexed by LJ parameter ID */
        QHash<quint64, QList<SireMM::LJException>> lj_exceptions;

        /** The indicies of the bonds for each molecule */
        QVector<QVector<int>> bonds_inc_h, bonds_exc_h;

        /** The indicies of the angles for each molecule */
        QVector<QVector<int>> angs_inc_h, angs_exc_h;

        /** The indicies of the dihedrals for each molecule */
        QVector<QVector<int>> dihs_inc_h, dihs_exc_h;

        /** The excluded atoms for each atom of each molecule */
        QVector<QVector<QVector<int>>> excl_atoms;

        /** The AtomNums of each atom in each molecule (indexed by MolNum) */
        QVector<QVector<int>> molnum_to_atomnums;

        /** A copy of the POINTER data to prevent over-lookup */
        QVector<qint64> pointers;

        /** The forcefield for the molecules in this file */
        SireMM::MMDetail ffield;

        /** Warnings when parsing this file */
        QStringList warns;

        /** The combining rules defined in this file */
        QString comb_rules;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::AmberPrm)

SIRE_EXPOSE_CLASS(SireIO::AmberPrm)

SIRE_END_HEADER

#endif
