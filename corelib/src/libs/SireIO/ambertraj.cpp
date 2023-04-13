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

#include "SireIO/ambertraj.h"
#include "SireIO/amberformat.h"

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

static const RegisterMetaType<AmberTraj> r_traj;
const RegisterParser<AmberTraj> register_traj;

QDataStream &operator<<(QDataStream &ds, const AmberTraj &traj)
{
    writeHeader(ds, r_traj, 1);

    SharedDataStream sds(ds);

    sds << traj.ttle << traj.natoms << traj.nframes << traj.nvalues << traj.has_box_dims
        << static_cast<const MoleculeParser &>(traj);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AmberTraj &traj)
{
    VersionID v = readHeader(ds, r_traj);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> traj.ttle >> traj.natoms >> traj.nframes >> traj.nvalues >> traj.has_box_dims >> static_cast<MoleculeParser &>(traj);
    }
    else
        throw version_error(v, "1", r_traj, CODELOC);

    return ds;
}

/** Constructor */
AmberTraj::AmberTraj()
    : ConcreteProperty<AmberTraj, MoleculeParser>(),
      natoms(0), nframes(0), nvalues(0), has_box_dims(false)
{
}

/** Return the format name that is used to identify this file format within Sire */
QString AmberTraj::formatName() const
{
    return "TRAJ";
}

/** Return the suffixes that RST7 files will typically have */
QStringList AmberTraj::formatSuffix() const
{
    static const QStringList suffixes = {"traj", "trj", "crd"};
    return suffixes;
}

/** Return a description of the file format */
QString AmberTraj::formatDescription() const
{
    return QObject::tr("Amber trajectory (ascii) coordinate or velocity files "
                       "supported from Amber 7 upwards.");
}

/** Scan the file to work out how many values there are,
 *  and to extrat the title
 */
void AmberTraj::parse()
{
    const auto &l = lines();

    if (l.count() < 2)
        // there is nothing in the file
        return;

    // read in the title FORMAT(20A4)
    ttle = l[0].simplified();

    nvalues = 0;
    natoms = 0;
    nframes = 0;
    has_box_dims = false;

    double linevals[8];

    for (int i = 0; i < 8; ++i)
    {
        linevals[i] = 0;
    }

    int nvals_per_line = 0;
    int nvals_last_line = -1;
    int nvals_per_frame = 0;

    auto end_of_frame = [&]()
    {
        // we must have read in the last line of what can only be
        // coordinate data for this frame
        nvals_per_frame += nvals_per_line;

        if (nvals_per_frame % 3 != 0)
        {
            throw SireIO::parse_error(QObject::tr(
                                          "This does not look like a valid Amber Traj file as the "
                                          "number of values in the frame (%1) is not divisible "
                                          "by three.")
                                          .arg(nvals_per_frame),
                                      CODELOC);
        }

        if (natoms <= 0)
        {
            natoms = nvals_per_frame / 3;
        }
        else if (nvals_per_frame / 3 != natoms)
        {
            throw SireIO::parse_error(QObject::tr(
                                          "This does not look like a valid Amber Traj file as the number "
                                          "of atoms in the first frame (%1) does not equal the number "
                                          "of atoms in a subsequent frame (%2)")
                                          .arg(natoms)
                                          .arg(nvals_per_frame / natoms),
                                      CODELOC);
        }

        nvals_per_frame = 0;

        nframes += 1;
    };

    // we will read the data, using the fact it is laid out in
    // blocks of 10 F8.3 values, with different blocks for the
    // atom coordinates (which must be 3 x atoms long) and
    // then box dimensions (a line of 3 x values long), repeating
    // for the number of frames, to work out the number of frame,
    // whether or not there are box dimensions, and the number
    // of atoms. This is not foolproof - there are edge cases
    // where we will read this as one big frame of a large number
    // of atoms, or we will count an atom as being a periodic
    // box (as we will lean towards interpreting a line with
    // three values as the box dimensions). In either of these
    // cases, we can fix the issue when we convert to a System

    for (int i = 1; i < l.count(); ++i)
    {
        const auto &line = l[i];
        const auto length = line.count();

        // read each line - this is a 10F8.3 format
        nvals_per_line = 0;
        for (int j = 0; j < 10; ++j)
        {
            int pos = j * 8;

            bool ok = true;
            double val = 0;

            if (pos + 8 > length)
            {
                ok = false;
            }
            else
            {
                val = line.midRef(pos, 8).toDouble(&ok);
            }

            if (not ok)
            {
                break;
            }

            linevals[j] = val;
            nvals_per_line += 1;
        }

        nvalues += nvals_per_line;

        if (nvals_per_line == 10)
        {
            // we have read in a full set of what can only be coordinate data
            nvals_per_frame += 10;
        }
        else if (nvals_per_line == 3)
        {
            // we must have either reached the end of a frame, or are
            // reading in the periodic box dimensions
            if (nvals_last_line == 3)
            {
                // this is the edge case where the number of coordinates
                // happens to have 3 space on the last line of the frame.
                // This line is the periodic box, and the previous line
                // is the end of the frame
                nvals_per_line = 0;
                end_of_frame();
                nvals_per_line = 3;
                has_box_dims = true;

                if (natoms == 1)
                    throw SireIO::parse_error(QObject::tr(
                                                  "This parser does not support the reading of single-atom "
                                                  "trajectories in Amber Traj format."),
                                              CODELOC);
            }
            else if (nvals_last_line == 10)
            {
                // this is the edge case where the number of coordinates
                // happens to fill the whole line (10 values) and this
                // line is the periodic box
                nvals_per_line = 0;
                end_of_frame();
                nvals_per_line = 3;
                has_box_dims = true;
            }
            else
            {
                // we must have already finished a frame
                if (nvals_per_frame != 0)
                    throw SireIO::parse_error(QObject::tr(
                                                  "This does not look like a valid Amber Traj file as the "
                                                  "number of consecutive lines without 10 values looks wrong."),
                                              CODELOC);

                has_box_dims = true;
            }
        }
        else
        {
            end_of_frame();
        }

        nvals_last_line = nvals_per_line;
    }

    this->setScore(natoms * nframes);
}

/** Construct by parsing the passed file */
AmberTraj::AmberTraj(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<AmberTraj, MoleculeParser>(filename, map),
      natoms(0), nframes(0), nvalues(0), has_box_dims(false)
{
    this->parse();
}

/** Construct by parsing the data in the passed text lines */
AmberTraj::AmberTraj(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<AmberTraj, MoleculeParser>(lines, map),
      natoms(0), nframes(0), nvalues(0), has_box_dims(false)
{
    this->parse();
}

static QStringList toLines(const QVector<QVector<Vector>> &all_coords,
                           QStringList *errors,
                           bool uses_parallel)
{
    // now find the start index of each molecule
    QVector<qint64> start_idx;
    start_idx.reserve(all_coords.count());

    qint64 last_idx = 0;

    for (const auto &molcoords : all_coords)
    {
        start_idx.append(last_idx);
        last_idx += 3 * molcoords.count();
    }

    const qint64 nats = last_idx / 3;
    QVector<double> coords(nats * 3, 0.0);

    if (uses_parallel)
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, all_coords.count()), [&](const tbb::blocked_range<int> &r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                const qint64 idx = start_idx.constData()[i];
                const auto molcoords = all_coords.constData()[i];

                for (int j = 0; j < molcoords.count(); ++j)
                {
                    const Vector &atomcoords = molcoords[j];
                    coords[idx + 3 * j + 0] = atomcoords.x();
                    coords[idx + 3 * j + 1] = atomcoords.y();
                    coords[idx + 3 * j + 2] = atomcoords.z();
                }
            } });
    }
    else
    {
        for (int i = 0; i < all_coords.count(); ++i)
        {
            const qint64 idx = start_idx.constData()[i];
            const auto molcoords = all_coords.constData()[i];

            for (int j = 0; j < molcoords.count(); ++j)
            {
                const Vector &atomcoords = molcoords[j];
                coords[idx + 3 * j + 0] = atomcoords.x();
                coords[idx + 3 * j + 1] = atomcoords.y();
                coords[idx + 3 * j + 2] = atomcoords.z();
            }
        }
    }

    QStringList lines = writeFloatData(coords, AmberFormat(AmberPrm::FLOAT, 10, 8, 3), errors, false, 'f');

    return lines;
}

static QVector<Vector> getCoordinates(const Molecule &mol, const PropertyName &coords_property)
{
    if (not mol.hasProperty(coords_property))
    {
        return QVector<Vector>();
    }

    QVector<Vector> coords(mol.nAtoms());

    const auto molcoords = mol.property(coords_property).asA<AtomCoords>();

    const auto molinfo = mol.info();

    for (int i = 0; i < mol.nAtoms(); ++i)
    {
        // coords are already in angstroms :-)
        coords[i] = molcoords.at(molinfo.cgAtomIdx(AtomIdx(i)));
    }

    return coords;
}

/** Construct by extracting the necessary data from the passed System */
AmberTraj::AmberTraj(const System &system, const PropertyMap &map)
    : ConcreteProperty<AmberTraj, MoleculeParser>(),
      natoms(0), nframes(0), nvalues(0), has_box_dims(false)
{
    // get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();

    if (molnums.isEmpty())
    {
        // no molecules in the system
        this->operator=(AmberTraj());
        return;
    }

    // get the coordinates (and velocities if available) for each molecule in the system
    QVector<QVector<Vector>> all_coords(molnums.count());

    const auto coords_property = map["coordinates"];

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, molnums.count()), [&](const tbb::blocked_range<int> r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                const auto mol = system[molnums[i]].molecule();
                all_coords[i] = ::getCoordinates(mol, coords_property);
            } });
    }
    else
    {
        for (int i = 0; i < molnums.count(); ++i)
        {
            const auto mol = system[molnums[i]].molecule();

            all_coords[i] = ::getCoordinates(mol, coords_property);
        }
    }

    QStringList errors;

    // extract the space of the system
    SpacePtr space;

    try
    {
        space = system.property(map["space"]).asA<Space>();
    }
    catch (...)
    {
    }

    // now convert these into text lines that can be written as the file
    QStringList lines = ::toLines(all_coords, &errors, usesParallel());

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error(
            QObject::tr("Errors converting the system to a Amber Traj format...\n%1").arg(errors.join("\n")), CODELOC);
    }

    // we don't need the coords and vels data any more, so free the memory
    all_coords.clear();

    lines.prepend(system.name().value());

    // finally add on the box dimensions and angles
    if (space.read().isA<PeriodicBox>())
    {
        Vector dims = space.read().asA<PeriodicBox>().dimensions();

        QVector<double> boxdims(3);
        boxdims[0] = dims.x();
        boxdims[1] = dims.y();
        boxdims[2] = dims.z();

        lines += writeFloatData(boxdims, AmberFormat(AmberPrm::FLOAT, 3, 8, 3), &errors, false, 'f');
    }
    else if (space.read().isA<Cartesian>())
    {
        // write an infinite box as three zero-length box dimensions
        QVector<double> boxdims(3, 0.0);
        lines += writeFloatData(boxdims, AmberFormat(AmberPrm::FLOAT, 3, 8, 3), &errors, false, 'f');
    }
    else
    {
        throw SireIO::parse_error(QObject::tr(
                                      "The Amber Traj format only supports the infinite cartesian or standard periodic box "
                                      "spaces. It doesn't support writing a system in the space %1")
                                      .arg(space.read().toString()),
                                  CODELOC);
    }

    if (not errors.isEmpty())
    {
        throw SireIO::parse_error(
            QObject::tr("Errors converting the system to a Amber Rst7 format...\n%1").arg(errors.join("\n")), CODELOC);
    }

    // now generate this object by re-reading these lines
    AmberTraj parsed(lines, map);

    this->operator=(parsed);
}

/** Copy constructor */
AmberTraj::AmberTraj(const AmberTraj &other)
    : ConcreteProperty<AmberTraj, MoleculeParser>(other),
      ttle(other.ttle), natoms(other.natoms), nframes(other.nframes),
      nvalues(other.nvalues), has_box_dims(other.has_box_dims)
{
}

/** Destructor */
AmberTraj::~AmberTraj()
{
}

AmberTraj &AmberTraj::operator=(const AmberTraj &other)
{
    if (this != &other)
    {
        ttle = other.ttle;
        natoms = other.natoms;
        nframes = other.nframes;
        nvalues = other.nvalues;
        has_box_dims = other.has_box_dims;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool AmberTraj::operator==(const AmberTraj &other) const
{
    return MoleculeParser::operator==(other);
}

bool AmberTraj::operator!=(const AmberTraj &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *AmberTraj::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberTraj>());
}

const char *AmberTraj::what() const
{
    return AmberTraj::typeName();
}

bool AmberTraj::isFrame() const
{
    return true;
}

int AmberTraj::nFrames() const
{
    return nframes;
}

Frame AmberTraj::getFrame(int frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    if (natoms <= 0)
        return Frame();

    if (frame < 0)
        frame = 0;

    // sanity check what we inferred when parsing
    qint32 nvalues_per_frame = natoms * 3;

    if (has_box_dims)
    {
        nvalues_per_frame += 3;
    }

    if (nframes * nvalues_per_frame != nvalues)
    {
        throw SireError::incompatible_error(QObject::tr(
                                                "The number of values in this trajectory (%1) is not "
                                                "compatible with the specified number of atoms (%2) and "
                                                "number of frames (%3).")
                                                .arg(nvalues)
                                                .arg(natoms)
                                                .arg(nframes),
                                            CODELOC);
    }

    qint32 start_idx = nvalues_per_frame * frame;
    qint32 end_idx = start_idx + nvalues_per_frame;

    if (end_idx > nvalues)
        throw SireError::invalid_index(QObject::tr(
                                           "It is not possible to read frame %1 as the total number of "
                                           "frames is only %2.")
                                           .arg(frame)
                                           .arg(nframes),
                                       CODELOC);

    const auto &l = lines();

    // this is a very simple format - just lots of doubles written
    // in a 10F8.3 format (not necessarily 10 numbers per line)
    // We will read them until we get to the specified start value,
    // and will then try to interpret them based on what we need

    int idx = 0;

    QVector<double> values;
    values.reserve(end_idx - start_idx);

    for (int i = 1; i < l.count(); ++i)
    {
        const auto &line = l[i];
        const auto length = line.count();

        bool line_ok = false;

        for (int j = 0; j < 10; ++j)
        {
            int pos = j * 8;

            if (pos + 8 > length)
                // nothing left on this line to read
                break;

            bool ok = false;
            double val = line.midRef(pos, 8).toDouble(&ok);

            if (not ok)
                // assume the rest of the line is corrupted
                break;

            if (idx >= start_idx and idx < end_idx)
                values.append(val);

            idx += 1;
            if (idx >= end_idx)
            {
                line_ok = false;
                break;
            }

            // we have read at least one value from this line
            line_ok = true;
        }

        if (not line_ok)
            break;
    }

    // ok, we now have the raw values - interpret them as coordinates
    auto coords = QVector<Vector>(natoms);
    auto coords_data = coords.data();

    const double *values_data = values.constData();

    for (int i = 0; i < natoms; ++i)
    {
        coords_data[i] = Vector(values_data[(3 * i) + 0],
                                values_data[(3 * i) + 1],
                                values_data[(3 * i) + 2]);
    }

    if (has_box_dims)
    {
        Vector box_dims(values_data[(3 * natoms) + 0],
                        values_data[(3 * natoms) + 1],
                        values_data[(3 * natoms) + 2]);

        if (box_dims.isZero())
        {
            // this is the infinite cartesian space
            return Frame(coords, Cartesian(), SireUnits::Dimension::Time(0));
        }
        else
        {
            return Frame(coords, PeriodicBox(box_dims), SireUnits::Dimension::Time(0));
        }
    }
    else
    {
        return Frame(coords, Cartesian(), SireUnits::Dimension::Time(0));
    }
}

QString AmberTraj::toString() const
{
    if (this->nAtoms() == 0)
    {
        return QObject::tr("AmberTraj::null");
    }
    else
    {
        return QObject::tr("AmberTraj( title() = %1, nAtoms() = %2, nFrames() = %3 )")
            .arg(this->title())
            .arg(this->nAtoms())
            .arg(this->nFrames());
    }
}

/** Parse from the passed file */
AmberTraj AmberTraj::parse(const QString &filename)
{
    return AmberTraj(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void AmberTraj::addToSystem(System &system, const PropertyMap &map) const
{
    // first, we are going to work with the group of all molecules, which should
    // be called "all". We have to assume that the molecules are ordered in "all"
    // in the same order as they are in this restart file, with the data
    // in MolIdx/AtomIdx order (this should be the default for all parsers!)
    MoleculeGroup allmols = system[MGName("all")];

    const int nmols = allmols.nMolecules();

    QVector<int> atom_pointers(nmols + 1, -1);

    int natoms = 0;

    for (int i = 0; i < nmols; ++i)
    {
        atom_pointers[i] = natoms;
        const int nats = allmols[MolIdx(i)].data().info().nAtoms();
        natoms += nats;
    }

    atom_pointers[nmols] = natoms;

    Frame frame = this->getFrame(0);

    if (natoms != frame.nAtoms())
    {
        // disagreement caused by us inferring the wrong number of atoms
        // when reading the trajectory. We will need to infer the right
        // number at a higher level...
        throw SireIO::parse_error(QObject::tr("Incompatibility between the files, as this trajectory file contains data "
                                              "for %1 atom(s), while the other file(s) have created a system with "
                                              "%2 atom(s)")
                                      .arg(frame.nAtoms())
                                      .arg(natoms),
                                  CODELOC);
    }

    // next, copy the coordinates into the molecules
    QVector<Molecule> mols(nmols);
    Molecule *mols_array = mols.data();
    const Vector *coords_array = frame.coordinates().constData();

    const PropertyName coords_property = map["coordinates"];

    auto add_coords = [&](int i)
    {
        const int atom_start_idx = atom_pointers.constData()[i];
        auto mol = system[MolIdx(i)].molecule();
        const auto molinfo = mol.data().info();

        // create space for the coordinates
        auto coords = QVector<QVector<Vector>>(molinfo.nCutGroups());

        for (int j = 0; j < molinfo.nCutGroups(); ++j)
        {
            coords[j] = QVector<Vector>(molinfo.nAtoms(CGIdx(j)));
        }

        for (int j = 0; j < mol.nAtoms(); ++j)
        {
            auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(j));

            const int atom_idx = atom_start_idx + j;

            coords[cgatomidx.cutGroup()][cgatomidx.atom()] = coords_array[atom_idx];
        }

        mols_array[i] = mol.edit().setProperty(coords_property, AtomCoords(CoordGroupArray(coords))).commit();
    };

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                add_coords(i);
            } });
    }
    else
    {
        for (int i = 0; i < nmols; ++i)
        {
            add_coords(i);
        }
    }

    system.update(Molecules(mols));

    PropertyName space_property = map["space"];
    if (space_property.hasValue())
    {
        system.setProperty("space", space_property.value());
    }
    else
    {
        system.setProperty(space_property.source(), frame.space());
    }

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

/** Return the title of the file */
QString AmberTraj::title() const
{
    return ttle;
}

/** Return the number of atoms whose coordinates are contained in this restart file */
int AmberTraj::nAtoms() const
{
    return natoms;
}

/** Return the parsed coordinate data */
QVector<SireMaths::Vector> AmberTraj::coordinates() const
{
    return this->getFrame(0).coordinates();
}

/** Return the parsed simulation box */
SireVol::SpacePtr AmberTraj::space() const
{
    return this->getFrame(0).space();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr AmberTraj::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberTraj(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr AmberTraj::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberTraj(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr AmberTraj::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberTraj(system, map));
}
