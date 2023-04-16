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

#include "SireIO/trr.h"
#include "SireIO/xdrfile.h"

#include "SireSystem/system.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/core.h"
#include "SireMol/mgname.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/molidx.h"
#include "SireMol/trajectory.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/generalunitproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/timeproperty.h"

#include "SireIO/errors.h"
#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMaths;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<TRR> r_traj;
const RegisterParser<TRR> register_traj;

QDataStream &operator<<(QDataStream &ds, const TRR &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);

    sds << static_cast<const MoleculeParser &>(traj);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, TRR &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> static_cast<MoleculeParser &>(traj);
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

/** Constructor */
TRR::TRR()
    : ConcreteProperty<TRR, MoleculeParser>(),
      nframes(0), frame_idx(0)
{
}

/** Return the format name that is used to identify this file format within Sire */
QString TRR::formatName() const
{
    return "TRR";
}

/** Return the suffixes that TRR files will typically have */
QStringList TRR::formatSuffix() const
{
    static const QStringList suffixes = {"trr"};
    return suffixes;
}

/** Return a description of the file format */
QString TRR::formatDescription() const
{
    return QObject::tr("Gromacs TRR (XDR file) coordinate / velocity / force trajectory file");
}

/** This is not a text file */
bool TRR::isTextFile() const
{
    return false;
}

/** Scan the file to work out how many values there are,
 *  and to extrat the title
 */
void TRR::parse()
{
    // open the XDR file and extract the data for the first frame
    TRRFile trr(this->filename());

    trr.open(QIODevice::ReadOnly);

    // read in the first frame - also find out how many frames there
    // are, and potentially build an index?

    trr.close();

    this->setScore(0);
}

/** Construct by parsing the passed file */
TRR::TRR(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<TRR, MoleculeParser>(map),
      nframes(0), frame_idx(0)
{
    // this gets the absolute file path
    this->setFilename(filename);
    this->parse();
}

/** Construct by parsing the data in the passed text lines */
TRR::TRR(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<TRR, MoleculeParser>(lines, map)
{
    throw SireIO::parse_error(QObject::tr("You cannot create a binary Gromacs TRR file from a set of text lines!"),
                              CODELOC);
}

/** Construct by extracting the necessary data from the passed System */
TRR::TRR(const System &system, const PropertyMap &map)
    : ConcreteProperty<TRR, MoleculeParser>(),
      nframes(1), frame_idx(0)
{
    current_frame = MoleculeParser::createFrame(system, map);
}

/** Copy constructor */
TRR::TRR(const TRR &other)
    : ConcreteProperty<TRR, MoleculeParser>(other),
      current_frame(other.current_frame), parse_warnings(other.parse_warnings),
      nframes(other.nframes), frame_idx(other.frame_idx)

{
}

/** Destructor */
TRR::~TRR()
{
}

TRR &TRR::operator=(const TRR &other)
{
    if (this != &other)
    {
        current_frame = other.current_frame;
        parse_warnings = other.parse_warnings;
        nframes = other.nframes;
        frame_idx = other.frame_idx;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool TRR::operator==(const TRR &other) const
{
    return MoleculeParser::operator==(other);
}

bool TRR::operator!=(const TRR &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *TRR::typeName()
{
    return QMetaType::typeName(qMetaTypeId<TRR>());
}

const char *TRR::what() const
{
    return TRR::typeName();
}

bool TRR::isFrame() const
{
    return true;
}

int TRR::nFrames() const
{
    return nframes;
}

Frame TRR::getFrame(int frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    if (frame < 0)
        frame = 0;

    if (frame == frame_idx)
        return current_frame;

    return Frame();
}

QString TRR::toString() const
{
    if (this->nAtoms() == 0)
    {
        return QObject::tr("TRR::null");
    }
    else
    {
        return QObject::tr("TRR( nAtoms() = %1, nFrames() = %2 )")
            .arg(this->nAtoms())
            .arg(this->nFrames());
    }
}

/** Parse from the passed file */
TRR TRR::parse(const QString &filename)
{
    return TRR(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void TRR::addToSystem(System &system, const PropertyMap &map) const
{
    MoleculeParser::copyFromFrame(current_frame, system, map);

    // update the System fileformat property to record that it includes
    // data from this file format
    QString fileformat = this->formatName();

    PropertyName fileformat_property = map["fileformat"];

    try
    {
        QString last_format = system.property(fileformat_property).asA<StringProperty>().value();
        fileformat = QString("%1,%2").arg(last_format, fileformat);
    }
    catch (...)
    {
    }

    if (fileformat_property.hasSource())
    {
        system.setProperty(fileformat_property.source(), StringProperty(fileformat));
    }
    else
    {
        system.setProperty("fileformat", StringProperty(fileformat));
    }
}

/** Return the number of atoms whose coordinates are contained in this restart file */
int TRR::nAtoms() const
{
    return current_frame.nAtoms();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr TRR::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(TRR(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr TRR::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(TRR(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr TRR::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(TRR(system, map));
}
