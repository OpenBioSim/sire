/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREIO_DCD_H
#define SIREIO_DCD_H

#include "moleculeparser.h"

#include "SireMol/trajectory.h"

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
    class DCD;

    namespace detail
    {
        class DCDFile;
    }

} // namespace SireIO

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::DCD &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::DCD &);

namespace SireIO
{
    class DCDFile;

    /** This class represents a DCD file reader.

        The format is described here

        https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html

        @author Christopher Woods
    */
    class SIREIO_EXPORT DCD : public SireBase::ConcreteProperty<DCD, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const DCD &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, DCD &);

    public:
        DCD();
        DCD(const QString &filename, const PropertyMap &map = PropertyMap());
        DCD(const QStringList &lines, const PropertyMap &map = PropertyMap());
        DCD(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        DCD(const DCD &other);

        ~DCD();

        DCD &operator=(const DCD &other);

        bool operator==(const DCD &other) const;
        bool operator!=(const DCD &other) const;

        static const char *typeName();

        const char *what() const;

        DCD operator[](int i) const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;

        QString formatName() const;
        QString formatDescription() const;
        QStringList formatSuffix() const;

        static DCD parse(const QString &filename);

        bool isFrame() const;

        int nFrames() const;
        int count() const;
        int size() const;

        SireMol::Frame getFrame(int i) const;

        int nAtoms() const;

        bool isTextFile() const;

        QStringList writeToFile(const QString &filename) const;

    protected:
        void addToSystem(SireSystem::System &system, const PropertyMap &map) const;

    private:
        void parse();

        /** The current frame */
        SireMol::Frame current_frame;

        /** Any warnings that were raised when reading the file */
        QStringList parse_warnings;

        /** The number of frames in this file */
        qint64 nframes;

        /** The current frame index */
        qint64 frame_idx;

        /** Pointer to the underlying file */
        std::shared_ptr<DCDFile> f;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::DCD)

SIRE_EXPOSE_CLASS(SireIO::DCD)

SIRE_END_HEADER

#endif
