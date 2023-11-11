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

#ifndef SIREIO_PDBX_H
#define SIREIO_PDBX_H

#include "moleculeparser.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
    class PDBx;
} // namespace SireIO

namespace SireMol
{
    class Atom;
    class MolEditor;
    class MoleculeInfoData;
    class MoleculeView;
    class Residue;
} // namespace SireMol

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::PDBx &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::PDBx &);

namespace SireIO
{
    /** This class holds a parser for reading and writing PDBx/mmcif files */
    class SIREIO_EXPORT PDBx : public SireBase::ConcreteProperty<PDBx, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const PDBx &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, PDBx &);

    public:
        PDBx();
        PDBx(const QString &filename, const PropertyMap &map = PropertyMap());

        PDBx(const QStringList &lines, const PropertyMap &map = PropertyMap());
        PDBx(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        PDBx(const PDBx &other);

        ~PDBx();

        PDBx &operator=(const PDBx &other);

        bool operator==(const PDBx &other) const;
        bool operator!=(const PDBx &other) const;

        static const char *typeName();

        const char *what() const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;
        QVector<QString> toLines() const;

        QString formatName() const;
        QString formatDescription() const;
        QStringList formatSuffix() const;

        bool isTopology() const;
        bool isFrame() const;

        int nAtoms() const;

        int nFrames() const;
        SireMol::Frame getFrame(int i) const;

    protected:
        SireSystem::System startSystem(const PropertyMap &map) const;
        void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

    private:
        void assertSane() const;
        void parseLines(const PropertyMap &map);

        /** Any warnings that were raised when reading the file. */
        QStringList parse_warnings;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::PDBx)

SIRE_EXPOSE_CLASS(SireIO::PDBx)

SIRE_END_HEADER

#endif
