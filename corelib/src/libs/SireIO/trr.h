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

#ifndef SIREIO_TRR_H
#define SIREIO_TRR_H

#include "moleculeparser.h"

#include "SireMol/trajectory.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
    class TRR;
}

SIREIO_EXPORT QDataStream &operator<<(QDataStream &, const SireIO::TRR &);
SIREIO_EXPORT QDataStream &operator>>(QDataStream &, SireIO::TRR &);

namespace SireIO
{
    class TRRFile;

    /** This class represents a Gromacs TRR (XDR file) trajectory.
     */
    class SIREIO_EXPORT TRR : public SireBase::ConcreteProperty<TRR, MoleculeParser>
    {

        friend SIREIO_EXPORT QDataStream & ::operator<<(QDataStream &, const TRR &);
        friend SIREIO_EXPORT QDataStream & ::operator>>(QDataStream &, TRR &);

    public:
        TRR();
        TRR(const QString &filename, const PropertyMap &map = PropertyMap());
        TRR(const QStringList &lines, const PropertyMap &map = PropertyMap());
        TRR(const SireSystem::System &system, const PropertyMap &map = PropertyMap());

        TRR(const TRR &other);

        ~TRR();

        TRR &operator=(const TRR &other);

        bool operator==(const TRR &other) const;
        bool operator!=(const TRR &other) const;

        static const char *typeName();

        const char *what() const;

        TRR operator[](int i) const;

        MoleculeParserPtr construct(const QString &filename, const PropertyMap &map) const;

        MoleculeParserPtr construct(const QStringList &lines, const PropertyMap &map) const;

        MoleculeParserPtr construct(const SireSystem::System &system, const PropertyMap &map) const;

        QString toString() const;

        QString formatName() const;
        QString formatDescription() const;
        QStringList formatSuffix() const;

        static TRR parse(const QString &filename);

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
        std::shared_ptr<TRRFile> f;
    };

} // namespace SireIO

Q_DECLARE_METATYPE(SireIO::TRR)

SIRE_EXPOSE_CLASS(SireIO::TRR)

SIRE_END_HEADER

#endif
