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

#include "SireIO/dcd.h"
#include "SireIO/fortranfile.h"

#include "SireSystem/system.h"

#include "SireIO/amberformat.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomforces.h"
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

#include "SireBase/booleanproperty.h"
#include "SireBase/generalunitproperty.h"
#include "SireBase/getinstalldir.h"
#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/timeproperty.h"
#include "SireBase/unittest.h"

#include "SireIO/errors.h"

#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDataStream>
#include <QDebug>
#include <QFile>
#include <QFileInfo>

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

//// Thanks to the MDAnalysis parser
//// (https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/lib/formats/include/readdcd.h)
//// which really helped with the reverse engineering of the DCD fileformat

/** This class provides a low-level interface to reading and writing
 *  a DCD file. It is designed to be used only with the
 *  DCD class
 */
class DCDFile
{
public:
    DCDFile();
    DCDFile(const QString &filename);
    ~DCDFile();

    bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly);

    SireMol::Frame readFrame(int i, bool use_parallel = true) const;
    void writeFrame(const SireMol::Frame &frame,
                    bool use_parallel = true);

    int nAtoms() const;
    int nFrames() const;

private:
    void readHeader();

    void _lkr_reset();
    void _lkr_reindexFrames();
    void _lkr_readFrameIntoBuffer(int i);
    void _lkr_writeBufferToFile();

    /** The current frame that has been read into the buffer */
    FortranFile f;

    QStringList title;

    QVector<qint32> fixed_atoms;

    /** The number of atoms in the frame - we assume all
     *  frames have the same number of atoms
     */
    qint64 natoms;

    /** The number of frames in the file */
    qint64 nframes;

    double timestep;

    SireVol::SpacePtr spc;

    qint64 istart;
    qint64 nsavc;
    qint64 nfixed;

    qint64 first_frame_line;

    bool CHARMM_FORMAT;
    bool HAS_EXTRA_BLOCK;
    bool HAS_FOUR_DIMS;
};

//// DCDFile

DCDFile::DCDFile()
    : natoms(0), nframes(0), timestep(0), istart(0), nsavc(0), nfixed(0), first_frame_line(0), CHARMM_FORMAT(false),
      HAS_EXTRA_BLOCK(false), HAS_FOUR_DIMS(false)
{
}

DCDFile::DCDFile(const QString &filename)
    : natoms(0), nframes(0), timestep(0), istart(0), nsavc(0), nfixed(0), first_frame_line(0), CHARMM_FORMAT(false),
      HAS_EXTRA_BLOCK(false), HAS_FOUR_DIMS(false)
{
    fFortranFile(filename);
    this->readHeader();
}

DCDFile::~DCDFile()
{
}

void SireIO::detail::DCDFile::readHeader(FortranFile &file)
{
    auto line = file[0];

    auto typ = line.readChar(4);

    if (typ != "CORD")
        throw SireIO::parse_error(QObject::tr("This does not look like a DCD file, because it does "
                                              "not start with 'CORD'. Starts with %1.")
                                      .arg(typ),
                                  CODELOC);

    auto ints = line.readInt32(9);

    nframes = ints[0];
    istart = ints[1];
    nsavc = ints[2];
    nfixed = ints[8];

    // now need to check the value at buffer[80] as, if this is non-zero,
    // then this is a CHARMM format DCD file
    CHARMM_FORMAT = line.readInt32At(80) != 0;

    // the value at buffer[44] says if there is an extra block
    HAS_EXTRA_BLOCK = line.readInt32At(44) != 0;

    // the value at buffer[48] says whether or not we have four dimensions
    HAS_FOUR_DIMS = line.readInt32At(48) != 0;

    timestep = 0;

    // read the timestep between frames (it is assumed to be in picoseconds)
    if (CHARMM_FORMAT)
    {
        timestep = line.readFloat32At(40);
    }
    else
    {
        timestep = line.readFloat64At(40);
    }

    line = file[1];

    int ntitle = line.readInt32(1)[0];

    for (int i = 0; i < ntitle; ++i)
    {
        QString t = line.readChar(32).simplified().replace(QChar('\0'), "");

        if (not t.isEmpty())
            title.append(t);
    }

    line = file[2];

    natoms = line.readInt32(1)[0];

    int linenum = 3;

    if (nfixed != 0)
    {
        line = file[linenum];
        linenum += 1;

        fixed_atoms = line.readInt32(nfixed);
    }

    first_frame_line = linenum;

    // now read in the space
    spc = this->readSpace(file, 0);

    if (nfixed != 0)
    {
        // we have to read in the first set of coordinates, as these
        // hold the fixed atoms as well as the movable atoms
        if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
        {
            line = file[linenum];
            linenum += 1;

            line.readFloat64(6);
        }

        line = file[linenum];
        linenum += 1;
        auto x = line.readFloat32(natoms);

        line = file[linenum];
        linenum += 1;
        auto y = line.readFloat32(natoms);

        line = file[linenum];
        linenum += 1;
        auto z = line.readFloat32(natoms);

        first_frame = QVector<Vector>(natoms);

        for (int i = 0; i < natoms; ++i)
        {
            first_frame[i] = Vector(x[i], y[i], z[i]);
        }
    }

    // now sanity check the rest of the file
    int num_lines_per_frame = 3;

    if (HAS_FOUR_DIMS)
    {
        num_lines_per_frame += 1;
    }

    if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
    {
        num_lines_per_frame += 1;
    }

    if (nframes != 0)
    {
        if (file.nRecords() != first_frame_line + (num_lines_per_frame * nframes))
        {
            throw SireIO::parse_error(QObject::tr("Wrong number of records in the DCD file. Expect to have %1 "
                                                  "for %2 frames, but actually have %3.")
                                          .arg(first_frame_line + (num_lines_per_frame * nframes))
                                          .arg(nframes)
                                          .arg(file.nRecords()),
                                      CODELOC);
        }
    }
    else
    {
        // we need to calculate nframes
        nframes = (file.nRecords() - first_frame_line) / num_lines_per_frame;
    }
}

double SireIO::detail::DCDFile::getTimeAtFrame(int frame) const
{
    return (istart * timestep) + (frame * timestep);
}

double SireIO::detail::DCDFile::getCurrentTime() const
{
    return getTimeAtFrame(0);
}

void SireIO::detail::DCDFile::setCurrentTime(double time)
{
    if (timestep != 0)
    {
        istart = int(time / timestep);
    }
    else
    {
        timestep = time;
        istart = 1;
    }
}

void SireIO::detail::DCDFile::setSpace(const Space &s)
{
    spc = s;
}

const Space &SireIO::detail::DCDFile::getSpace() const
{
    return *spc;
}

SpacePtr SireIO::detail::DCDFile::readSpace(FortranFile &file, int frame) const
{
    if (frame < 0 or frame >= nframes)
    {
        throw SireIO::parse_error(
            QObject::tr("Trying to access an invalid frame (%1) from a DCD with %2 frames.").arg(frame).arg(nframes),
            CODELOC);
    }

    // get the line number for this frame
    int num_lines_per_frame = 3;

    if (HAS_FOUR_DIMS)
        num_lines_per_frame += 1;

    if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
    {
        num_lines_per_frame += 1;

        int linenum = first_frame_line + (frame * num_lines_per_frame);

        auto line = file[linenum];

        auto boxinfo = line.readFloat64(6);

        // construct from the above boxinfo
        return SpacePtr();
    }

    return SpacePtr();
}

QVector<Vector> SireIO::detail::DCDFile::readCoordinates(FortranFile &file, int frame) const
{
    if (frame < 0 or frame >= nframes)
    {
        throw SireIO::parse_error(
            QObject::tr("Trying to access an invalid frame (%1) from a DCD with %2 frames.").arg(frame).arg(nframes),
            CODELOC);
    }

    // get the line number for this frame
    int num_lines_per_frame = 3;
    int skip_unitcell = 0;

    if (HAS_FOUR_DIMS)
        num_lines_per_frame += 1;

    if (CHARMM_FORMAT and HAS_EXTRA_BLOCK)
    {
        num_lines_per_frame += 1;
        skip_unitcell = 1;
    }

    int linenum = first_frame_line + (frame * num_lines_per_frame) + skip_unitcell;

    if (nfixed == 0)
    {
        // read in the x, y, and z data
        auto line = file[linenum];
        auto x = line.readFloat32(natoms);
        line = file[linenum + 1];
        auto y = line.readFloat32(natoms);
        line = file[linenum + 2];
        auto z = line.readFloat32(natoms);

        QVector<Vector> coords(natoms);

        for (int i = 0; i < natoms; ++i)
        {
            coords[i] = Vector(x[i], y[i], z[i]);
        }

        return coords;
    }
    else if (frame == 0)
    {
        return first_frame;
    }
    else
    {
        // read the coordinates and map them into a copy of first_frame
        QVector<Vector> frame = first_frame;

        return frame;
    }
}

Frame SireIO::detail::DCDFile::readFrame(FortranFile &file, int frame) const
{
    frame = SireID::Index(frame).map(nframes);

    auto space = this->readSpace(file, frame);
    auto coords = this->readCoordinates(file, frame);

    return Frame(coords, space, getTimeAtFrame(frame) * picosecond);
}

QString SireIO::detail::DCDFile::getTitle() const
{
    return title.join("");
}

void SireIO::detail::DCDFile::setTitle(QString t)
{
    // need to split into blocks of 32 characters
    title.clear();

    while (t.length() > 32)
    {
        title.append(t.mid(0, 32));
        t.remove(0, 32);
    }

    if (t.length() > 0)
    {
        title.append(t);
    }
}

double SireIO::detail::DCDFile::getTimeStep() const
{
    return timestep;
}

qint64 SireIO::detail::DCDFile::getFrameStart() const
{
    return istart;
}

qint64 SireIO::detail::DCDFile::getFrameDelta() const
{
    return nsavc;
}

qint64 SireIO::detail::DCDFile::nAtoms() const
{
    return natoms;
}

qint64 SireIO::detail::DCDFile::nFrames() const
{
    return nframes;
}

////////
//////// Implemenetation of DCD
////////

static const RegisterMetaType<DCD> r_dcd;
const RegisterParser<DCD> register_dcd;

QDataStream &operator<<(QDataStream &ds, const DCD &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);

    sds << traj.current_frame << traj.parse_warnings
        << traj.nframes << traj.frame_idx
        << static_cast<const MoleculeParser &>(traj);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, DCD &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> traj.current_frame >> traj.parse_warnings >> traj.nframes >> traj.frame_idx >> static_cast<MoleculeParser &>(traj);
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

/** Constructor */
DCD::DCD()
    : ConcreteProperty<DCD, MoleculeParser>(),
      nframes(0), frame_idx(0)
{
}

/** Return the format name that is used to identify this file format within Sire */
QString DCD::formatName() const
{
    return "DCD";
}

/** Return the suffixes that DCD files will typically have */
QStringList DCD::formatSuffix() const
{
    static const QStringList suffixes = {"DCD"};
    return suffixes;
}

/** Return a description of the file format */
QString DCD::formatDescription() const
{
    return QObject::tr("DCD coordinate/velocity binary trajectory files "
                       "based on charmm / namd / x-plor format.");
}

/** This is not a text file */
bool DCD::isTextFile() const
{
    return false;
}

/** Open the file and read in all the metadata */
void DCD::parse()
{
    f.reset(new DCDFile(this->filename()));

    try
    {
        if (not f->open(QIODevice::ReadOnly))
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Failed to open DCDFile %1")
                                          .arg(this->filename()),
                                      CODELOC);
        }

        nframes = f->nFrames();
        current_frame = f->readFrame(0, this->usesParallel());
        frame_idx = 0;

        this->setScore(f->nFrames() * current_frame.nAtoms());
    }
    catch (...)
    {
        this->setScore(0);
        f.reset();
        throw;
    }
}

/** Construct by parsing the passed file */
DCD::DCD(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<DCD, MoleculeParser>(map),
      nframes(0), frame_idx(0)
{
    // this gets the absolute file path
    this->setFilename(filename);
    this->parse();
}

/** Construct by parsing the data in the passed text lines */
DCD::DCD(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<DCD, MoleculeParser>(lines, map)
{
    throw SireIO::parse_error(QObject::tr("You cannot create a binary Gromacs DCD file from a set of text lines!"),
                              CODELOC);
}

/** Construct by extracting the necessary data from the passed System */
DCD::DCD(const System &system, const PropertyMap &map)
    : ConcreteProperty<DCD, MoleculeParser>(system, map),
      nframes(1), frame_idx(0)
{
    current_frame = MoleculeParser::createFrame(system, map);
}

/** Copy constructor */
DCD::DCD(const DCD &other)
    : ConcreteProperty<DCD, MoleculeParser>(other),
      current_frame(other.current_frame), parse_warnings(other.parse_warnings),
      nframes(other.nframes), frame_idx(other.frame_idx), f(other.f)
{
}

/** Destructor */
DCD::~DCD()
{
}

DCD &DCD::operator=(const DCD &other)
{
    if (this != &other)
    {
        current_frame = other.current_frame;
        parse_warnings = other.parse_warnings;
        nframes = other.nframes;
        frame_idx = other.frame_idx;
        f = other.f;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool DCD::operator==(const DCD &other) const
{
    return MoleculeParser::operator==(other);
}

bool DCD::operator!=(const DCD &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *DCD::typeName()
{
    return QMetaType::typeName(qMetaTypeId<DCD>());
}

const char *DCD::what() const
{
    return DCD::typeName();
}

bool DCD::isFrame() const
{
    return true;
}

int DCD::nFrames() const
{
    return nframes;
}

int DCD::count() const
{
    return this->nFrames();
}

int DCD::size() const
{
    return this->nFrames();
}

Frame DCD::getFrame(int frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    if (frame < 0)
        frame = 0;

    if (frame == frame_idx)
        return current_frame;

    if (f.get() == 0)
    {
        throw SireError::file_error(QObject::tr(
                                        "Somehow we don't have access to the underlying DCD file?"),
                                    CODELOC);
    }

    return f->readFrame(frame, this->usesParallel());
}

DCD DCD::operator[](int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    DCD ret(*this);

    ret.current_frame = this->getFrame(i);
    ret.frame_idx = i;

    return ret;
}

QString DCD::toString() const
{
    if (this->nAtoms() == 0)
    {
        return QObject::tr("DCD::null");
    }
    else
    {
        return QObject::tr("DCD( nAtoms() = %1, nFrames() = %2 )")
            .arg(this->nAtoms())
            .arg(this->nFrames());
    }
}

/** Parse from the passed file */
DCD DCD::parse(const QString &filename)
{
    return DCD(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void DCD::addToSystem(System &system, const PropertyMap &map) const
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
int DCD::nAtoms() const
{
    return current_frame.nAtoms();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr DCD::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(DCD(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr DCD::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(DCD(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr DCD::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(DCD(system, map));
}

/** Write this binary file 'filename' */
QStringList DCD::writeToFile(const QString &filename) const
{
    if (this->nFrames() == 0 or this->nAtoms() == 0)
        return QStringList();

    auto gil = SireBase::release_gil();

    createDirectoryForFile(filename);

    DCDFile outfile(filename);

    if (not outfile.open(QIODevice::WriteOnly))
        throw SireError::file_error(QObject::tr(
                                        "Could not open %1 to write the DCD file.")
                                        .arg(filename),
                                    CODELOC);

    if (this->writingTrajectory())
    {
        const auto frames = this->framesToWrite();

        ProgressBar bar("Save DCD", frames.count());
        bar.setSpeedUnit("frames / s");

        bar = bar.enter();

        for (int i = 0; i < frames.count(); ++i)
        {
            const auto frame = this->createFrame(frames[i]);
            outfile.writeFrame(frame, usesParallel());
            bar.setProgress(i + 1);
        }

        bar.success();
    }
    else
    {
        outfile.writeFrame(current_frame, usesParallel());
    }

    outfile.close();

    return QStringList(filename);
}
