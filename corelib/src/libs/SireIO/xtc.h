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

#ifndef SIREIO_XTC_H
#define SIREIO_XTC_H

#include "moleculeparser.h"

#include "SireMol/trajectory.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
    class XTC;
}

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::XTC &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::XTC &);

namespace SireIO
{
    class XTCFile;

    /** This class represents a Gromacs XTC (XDR file) trajectory.
     *  This is a compressed trajectory format
     */
    class SIREIO_EXPORT XTC : public SireBase::ConcreteProperty<XTC, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const XTC &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, XTC &);

    public:
        XTC();
        XTC(const QString &filename, const PropertyMap &map = PropertyMap());
        XTC(const QStringList &lines, const PropertyMap &map = PropertyMap());
        XTC(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        XTC(const XTC &other);

        ~XTC();

        XTC &operator=(const XTC &other);

        bool operator==(const XTC &other) const;
        bool operator!=(const XTC &other) const;

        static const char *typeName();

        const char *what() const;

        XTC operator[](int i) const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;

        QString formatName() const;
        QString formatDescription() const;
        QStringList formatSuffix() const;

        static XTC parse(const QString &filename);

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
        std::shared_ptr<XTCFile> f;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::XTC)

SIRE_EXPOSE_CLASS(SireIO::XTC)

SIRE_END_HEADER

#endif
