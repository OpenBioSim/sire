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

#include "moleculeparser.h"
#include "filetrajectory.h"
#include "filetrajectoryparser.h"
#include "supplementary.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireBase/booleanproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/releasegil.h"

#include "SireFF/ffdetail.h"
#include "SireMM/mmdetail.h"

#include "SireMol/core.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/trajectory.h"
#include "SireMol/mgname.h"
#include "SireMol/mgnum.h"
#include "SireMol/molidx.h"
#include "SireMol/trajectoryaligner.h"

#include "SireBase/timeproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/propertylist.h"
#include "SireBase/releasegil.h"
#include "SireBase/progressbar.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QMutex>
#include <QTextStream>

using namespace SireIO;
using namespace SireSystem;
using namespace SireMaths;
using namespace SireMol;
using namespace SireFF;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits;

//////////////
////////////// Implementation of ParserFactory and ParserFactoryHelper
//////////////

namespace SireIO
{
    namespace detail
    {
        /** Null constructor */
        ParserFactoryHelper::ParserFactoryHelper()
        {
        }

        /** Copy constructor */
        ParserFactoryHelper::ParserFactoryHelper(const ParserFactoryHelper &other) : parser(other.parser)
        {
        }

        /** Destructor */
        ParserFactoryHelper::~ParserFactoryHelper()
        {
        }

        ParserFactoryHelper &ParserFactoryHelper::operator=(const ParserFactoryHelper &other)
        {
            parser = other.parser;
            return *this;
        }

        bool ParserFactoryHelper::operator<(const ParserFactoryHelper &other) const
        {
            if (isValid())
            {
                if (other.isValid())
                {
                    return parser->formatName() < other.parser->formatName();
                }
                else
                {
                    return true;
                }
            }
            else
            {
                return not other.isValid();
            }
        }

        bool ParserFactoryHelper::operator==(const ParserFactoryHelper &other) const
        {
            if (isValid())
            {
                if (other.isValid())
                {
                    return parser->formatName() == other.parser->formatName();
                }
                else
                    return false;
            }
            else
            {
                return not other.isValid();
            }
        }

        bool ParserFactoryHelper::operator>(const ParserFactoryHelper &other) const
        {
            return not(operator==(other) or operator<(other));
        }

        bool ParserFactoryHelper::operator!=(const ParserFactoryHelper &other) const
        {
            return not operator==(other);
        }

        bool ParserFactoryHelper::operator<=(const ParserFactoryHelper &other) const
        {
            return operator==(other) or operator<(other);
        }

        bool ParserFactoryHelper::operator>=(const ParserFactoryHelper &other) const
        {
            return not operator<(other);
        }

        /** Return whether or not this helper is valid */
        bool ParserFactoryHelper::isValid() const
        {
            return parser.get() != 0;
        }

        /** Return whether or not this is the supplementary parser */
        bool ParserFactoryHelper::isSupplementary() const
        {
            return parser.get() != 0 and parser->isA<Supplementary>();
        }

        /** Return the unique ID name of the parser in the program */
        QString ParserFactoryHelper::formatName() const
        {
            if (isValid())
            {
                return parser->formatName();
            }
            else
            {
                return QString();
            }
        }

        /** Return the description of the parser */
        QString ParserFactoryHelper::formatDescription() const
        {
            if (isValid())
            {
                return parser->formatDescription();
            }
            else
            {
                return QString();
            }
        }

        QString ParserFactoryHelper::toString() const
        {
            return QString("Parser( %1 : %2 )").arg(formatName()).arg(formatDescription());
        }

        /** Return all of the suffixes recognised by this parser, in their order
            of preference */
        QStringList ParserFactoryHelper::suffixes() const
        {
            if (isValid())
            {
                return parser->formatSuffix();
            }
            else
            {
                return QStringList();
            }
        }

        /** Return the preferred suffix for the parser */
        QString ParserFactoryHelper::preferredSuffix() const
        {
            const auto s = this->suffixes();

            if (not s.isEmpty())
            {
                return s[0];
            }
            else
            {
                return QString();
            }
        }

        /** Use this factory helper to construct a new parser that parses
            the file called 'filename' */
        MoleculeParserPtr ParserFactoryHelper::construct(const QString &filename, const PropertyMap &map) const
        {
            if (isValid())
            {
                return parser->construct(filename, map);
            }
            else
                return MoleculeParserPtr();
        }

        /** Use this factory helper to construct a new parser that parses
            the data in the passed lines of text */
        MoleculeParserPtr ParserFactoryHelper::construct(const QStringList &lines, const PropertyMap &map) const
        {
            if (isValid())
            {
                return parser->construct(lines, map);
            }
            else
                return MoleculeParserPtr();
        }

        /** Use this factory helper to construct a new parser from the information
            contained in the passed system */
        MoleculeParserPtr ParserFactoryHelper::construct(const SireSystem::System &system, const PropertyMap &map) const
        {
            if (isValid())
            {
                return parser->construct(system, map);
            }
            else
                return MoleculeParserPtr();
        }

        /** The parser factory */
        class ParserFactory
        {
        public:
            ParserFactory()
            {
            }

            ~ParserFactory()
            {
            }

            void registerParser(const ParserFactoryHelper &helper)
            {
                if (not helper.isValid())
                {
                    return;
                }

                QMutexLocker lkr(&mutex);

                helpers_by_id.insert(helper.formatName().toLower(), helper);

                for (const auto &suffix : helper.suffixes())
                {
                    helpers_by_suffix.insert(suffix.toLower(), helper);
                }
            }

            QList<ParserFactoryHelper> getFactories(const QStringList &parser_names)
            {
                QMutexLocker lkr(&mutex);
                QList<ParserFactoryHelper> helpers;
                QStringList missing;

                for (const auto &name : parser_names)
                {
                    const auto lower_name = name.toLower();
                    auto helper = helpers_by_id.value(lower_name);

                    if (not helper.isValid())
                    {
                        // search for the helped in the secondary suffixes...
                        for (auto h : helpers_by_id.values())
                        {
                            for (auto suffix : h.suffixes())
                            {
                                if (name.toLower() == suffix.toLower())
                                {
                                    helper = h;
                                    break;
                                }
                            }

                            if (helper.isValid())
                                break;
                        }
                    }

                    helpers.append(helper);

                    if (not helpers.last().isValid())
                    {
                        missing.append(name);
                    }
                }

                if (not missing.isEmpty())
                {
                    lkr.unlock();
                    throw SireError::io_error(QObject::tr("Cannot find parsers that support the following formats: %1.\n"
                                                          "Supported parsers are:\n%2")
                                                  .arg(missing.join(", "))
                                                  .arg(this->supportedFormats()),
                                              CODELOC);
                }

                return helpers;
            }

            QList<ParserFactoryHelper> factoriesForSuffix(const QString &suffix, bool disable_supplementary)
            {
                QMutexLocker lkr(&mutex);
                auto helpers = helpers_by_suffix.values(suffix);
                std::sort(helpers.begin(), helpers.end());

                if (disable_supplementary)
                {
                    QMutableListIterator<ParserFactoryHelper> it(helpers);

                    while (it.hasNext())
                    {
                        const auto &value = it.next();

                        if (value.isSupplementary())
                        {
                            it.remove();
                        }
                    }
                }

                return helpers;
            }

            QList<ParserFactoryHelper> factoriesExcludingSuffix(const QString &suffix, bool disable_supplementary)
            {
                QMutexLocker lkr(&mutex);

                if (suffix.isEmpty())
                {
                    auto helpers = helpers_by_id.values();
                    std::sort(helpers.begin(), helpers.end());

                    if (disable_supplementary)
                    {
                        QMutableListIterator<ParserFactoryHelper> it(helpers);

                        while (it.hasNext())
                        {
                            const auto &value = it.next();

                            if (value.isSupplementary())
                            {
                                it.remove();
                            }
                        }
                    }

                    return helpers;
                }

                QList<ParserFactoryHelper> helpers;

                for (const auto &helper : helpers_by_id)
                {
                    if (not helper.suffixes().contains(suffix))
                    {
                        helpers.append(helper);
                    }
                }

                std::sort(helpers.begin(), helpers.end());

                if (disable_supplementary)
                {
                    QMutableListIterator<ParserFactoryHelper> it(helpers);

                    while (it.hasNext())
                    {
                        const auto &value = it.next();

                        if (value.isSupplementary())
                        {
                            it.remove();
                        }
                    }
                }

                return helpers;
            }

            ParserFactoryHelper factory(const QString &name)
            {
                QMutexLocker lkr(&mutex);
                return helpers_by_id.value(name);
            }

            QString supportedFormats()
            {
                QMutexLocker lkr(&mutex);

                auto keys = helpers_by_id.keys();
                std::sort(keys.begin(), keys.end());

                QStringList lines;

                for (const auto &key : keys)
                {
                    const auto parser = helpers_by_id.value(key);

                    if (not parser.isSupplementary())
                    {
                        lines.append(QObject::tr("## Parser %1 ##").arg(key));
                        lines.append(QObject::tr("Supports files: %1").arg(parser.suffixes().join(", ")));
                        lines.append(parser.formatDescription());
                        lines += QString("#").repeated(13 + key.length()) + "\n";
                    }
                }

                return lines.join("\n");
            }

        private:
            /** Mutex to serialise access to the factory */
            QMutex mutex;

            /** All of the factory helpers arranged by the suffix of
                file that they support */
            QMultiHash<QString, ParserFactoryHelper> helpers_by_suffix;

            /** All of the factor helpers arranged by their unique ID */
            QHash<QString, ParserFactoryHelper> helpers_by_id;
        };

    } // end of namespace detail
} // end of namespace SireIO

Q_GLOBAL_STATIC(SireIO::detail::ParserFactory, getParserFactory);

/** This registers a ParserFactoryHelper with the ParserFactory for the
    specified parser */
SireIO::detail::ParserFactoryHelper::ParserFactoryHelper(MoleculeParser *p)
{
    parser.reset(p);
    getParserFactory()->registerParser(*this);
}

//////////////
////////////// Implementation of MoleculeParser
//////////////

static const RegisterMetaType<MoleculeParser> r_parser(MAGIC_ONLY, MoleculeParser::typeName());

QDataStream &operator<<(QDataStream &ds, const MoleculeParser &parser)
{
    writeHeader(ds, r_parser, 4);

    SharedDataStream sds(ds);
    sds << parser.fname << parser.lnes
        << parser.saved_system
        << parser.loaded_order
        << parser.frames_to_write
        << parser.propmap
        << parser.scr << parser.run_parallel << static_cast<const Property &>(parser);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MoleculeParser &parser)
{
    VersionID v = readHeader(ds, r_parser);

    parser.loaded_order.clear();

    if (v == 4)
    {
        SharedDataStream sds(ds);
        sds >> parser.fname >> parser.lnes >> parser.saved_system >> parser.loaded_order >> parser.frames_to_write >> parser.propmap >> parser.scr >> parser.run_parallel >> static_cast<Property &>(parser);
    }
    else if (v == 3)
    {
        SharedDataStream sds(ds);
        sds >> parser.fname >> parser.lnes >> parser.saved_system >> parser.frames_to_write >> parser.propmap >> parser.scr >> parser.run_parallel >> static_cast<Property &>(parser);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);
        sds >> parser.fname >> parser.lnes >> parser.scr >> parser.run_parallel >> static_cast<Property &>(parser);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> parser.lnes >> parser.scr >> parser.run_parallel >> static_cast<Property &>(parser);

        parser.fname = QString();
    }
    else
        throw version_error(v, "1, 2, 3, 4", r_parser, CODELOC);

    return ds;
}

/** Constructor */
MoleculeParser::MoleculeParser(const PropertyMap &map) : Property(), propmap(map), scr(0), run_parallel(true)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }
}

void MoleculeParser::setParsedSystem(const System &system, const PropertyMap &map)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (map["frames_to_write"].hasValue())
    {
        const qint32 num_frames = system.nFrames();

        auto frames = map["frames_to_write"].value().asAnArray();

        for (int i = 0; i < frames.count(); ++i)
        {
            qint32 frame_idx = frames[i].asAnInteger();

            // make sure this is a valid frame
            if (not Index(frame_idx).canMap(num_frames))
            {
                throw SireError::invalid_index(QObject::tr(
                                                   "Cannot save trajectory frame '%1' as the number of "
                                                   "frames is only '%2'")
                                                   .arg(frame_idx)
                                                   .arg(num_frames),
                                               CODELOC);
            }

            frames_to_write.append(frame_idx);
        }
    }

    saved_system = system;
    propmap = map;
}

/** Constructor */
MoleculeParser::MoleculeParser(const System &system, const PropertyMap &map)
    : Property(), saved_system(system), propmap(map), scr(0), run_parallel(true)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (map["frames_to_write"].hasValue())
    {
        const qint32 num_frames = system.nFrames();

        auto frames = map["frames_to_write"].value().asAnArray();

        for (int i = 0; i < frames.count(); ++i)
        {
            qint32 frame_idx = frames[i].asAnInteger();

            // make sure this is a valid frame
            if (not Index(frame_idx).canMap(num_frames))
            {
                throw SireError::invalid_index(QObject::tr(
                                                   "Cannot save trajectory frame '%1' as the number of "
                                                   "frames is only '%2'")
                                                   .arg(frame_idx)
                                                   .arg(num_frames),
                                               CODELOC);
            }

            frames_to_write.append(frame_idx);
        }
    }
}

/** Internal function that provides a file cache */
class FileContentsCache
{
public:
    FileContentsCache()
    {
    }

    ~FileContentsCache()
    {
    }

    QVector<QString> read(QString filename)
    {
        QMutexLocker lkr(&cache_mutex);
        auto it = cache.constFind(filename);

        if (it != cache.constEnd())
        {
            return it.value();
        }
        else
        {
            return QVector<QString>();
        }
    }

    void save(QString filename, QVector<QString> filecontents)
    {
        QMutexLocker lkr(&cache_mutex);
        cache[filename] = filecontents;
    }

    void clear()
    {
        QMutexLocker lkr(&cache_mutex);
        cache.clear();
    }

private:
    QMutex cache_mutex;
    QHash<QString, QVector<QString>> cache;
};

Q_GLOBAL_STATIC(FileContentsCache, getFileCache);

/** Internal function that can be used by the parsers to read the contents
    of a text file into memory. This uses a cache to ensure that every file
    is read only once */
QVector<QString> MoleculeParser::readTextFile(QString filename)
{
    filename = QFileInfo(filename).absoluteFilePath();

    QVector<QString> lines = getFileCache()->read(filename);

    if (not lines.isEmpty())
        return lines;

    QFile file(filename);

    if (not file.open(QIODevice::ReadOnly | QIODevice::Unbuffered))
    {
        throw SireError::file_error(file, CODELOC);
    }

    QTextStream ts(&file);

    QStringList l;

    while (not ts.atEnd())
    {
        l.append(ts.readLine());
    }

    file.close();

    lines = l.toVector();

    if (not lines.isEmpty())
        getFileCache()->save(filename, lines);

    return lines;
}

/** Construct the parser, parsing in all of the lines in the file
    with passed filename */
MoleculeParser::MoleculeParser(const QString &filename, const PropertyMap &map) : Property(), propmap(map), scr(0), run_parallel(true)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    fname = QFileInfo(filename).absoluteFilePath();
    lnes = readTextFile(fname);
}

/** Construct the parser, parsing in all of the passed text lines */
MoleculeParser::MoleculeParser(const QStringList &lines, const PropertyMap &map)
    : Property(), propmap(map), scr(0), run_parallel(true)
{
    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (not lines.isEmpty())
    {
        lnes = lines.toVector();
    }
}

/** Copy constructor */
MoleculeParser::MoleculeParser(const MoleculeParser &other)
    : Property(other), fname(other.fname), lnes(other.lnes),
      saved_system(other.saved_system),
      loaded_order(other.loaded_order), frames_to_write(other.frames_to_write),
      propmap(other.propmap), scr(other.scr), run_parallel(other.run_parallel)
{
}

/** Destructor */
MoleculeParser::~MoleculeParser()
{
}

/** Return whether or not this parser runs in parallel - this will depend
 *  on whether parallel was enabled for parsing, we have more than
 *  one thread available, and the number of items (n) makes it worth
 *  parsing in parallel. Pass in n<=0 if you don't want to check the
 *  last condition.
 */
bool MoleculeParser::usesParallel(int n) const
{
    if (n > 0 and n < 8)
    {
        // don't parallelise for 8 items
        return false;
    }

    // how many threads are available
    if (get_max_num_threads() <= 1)
        return false;

    return run_parallel;
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
        coords[i] = molcoords.at(molinfo.cgAtomIdx(AtomIdx(i)));
    }

    return coords;
}

static QVector<Velocity3D> getVelocities(const Molecule &mol, const PropertyName &vels_property)
{
    if (not mol.hasProperty(vels_property))
    {
        return QVector<Velocity3D>();
    }

    try
    {
        const auto molvels = mol.property(vels_property).asA<AtomVelocities>();
        const auto molinfo = mol.info();

        QVector<Velocity3D> vels(mol.nAtoms());

        for (int i = 0; i < mol.nAtoms(); ++i)
        {
            vels[i] = molvels.at(molinfo.cgAtomIdx(AtomIdx(i)));
        }

        return vels;
    }
    catch (...)
    {
        return QVector<Velocity3D>();
    }
}

static QVector<Force3D> getForces(const Molecule &mol, const PropertyName &forces_property)
{
    if (not mol.hasProperty(forces_property))
    {
        return QVector<Force3D>();
    }

    try
    {
        const auto molforces = mol.property(forces_property).asA<AtomForces>();
        const auto molinfo = mol.info();

        QVector<Force3D> forces(mol.nAtoms());

        for (int i = 0; i < mol.nAtoms(); ++i)
        {
            forces[i] = molforces.at(molinfo.cgAtomIdx(AtomIdx(i)));
        }

        return forces;
    }
    catch (...)
    {
        return QVector<Force3D>();
    }
}

template <class T>
static bool hasData(const QVector<QVector<T>> &array)
{
    const auto array_data = array.constData();

    for (int i = 0; i < array.count(); ++i)
    {
        if (not array_data[i].isEmpty())
            return true;
    }

    return false;
}

/** Internal function used to collapse an array of arrays of type T into
    a single array of type T */
template <class T>
static QVector<T> collapse(const QVector<QVector<T>> &arrays)
{
    int nvals = 0;

    for (const auto &array : arrays)
    {
        nvals += array.count();
    }

    if (nvals == 0)
    {
        return QVector<T>();
    }

    QVector<T> values;
    values.reserve(nvals);

    for (const auto &array : arrays)
    {
        values += array;
    }

    return values;
}

static SireUnits::Dimension::Time get_time_from_system(const System &system,
                                                       const PropertyName &time_property)
{
    SireUnits::Dimension::Time time(0);

    if (system.containsProperty(time_property))
    {
        try
        {
            const Property &prop = system.property(time_property);

            if (prop.isA<TimeProperty>())
                time = prop.asA<TimeProperty>().value();
            else
                time = prop.asA<GeneralUnitProperty>().toUnit<SireUnits::Dimension::Time>();
        }
        catch (...)
        {
        }
    }

    return time;
}

/** Convenience function that converts the passed System into a Frame */
Frame MoleculeParser::createFrame(const System &system,
                                  const PropertyMap &map) const
{
    // get the MolNums of each molecule in the System - this returns the
    // numbers in MolIdx order
    const QVector<MolNum> molnums = system.getMoleculeNumbers().toVector();
    const auto molnums_data = molnums.constData();

    if (molnums.isEmpty())
    {
        // no molecules in the system
        return Frame();
    }

    QVector<Vector> coords;
    QVector<Velocity3D> vels;
    QVector<Force3D> frcs;

    // get the coordinates (and velocities if available) for each molecule in the system
    {
        QVector<QVector<Vector>> all_coords(molnums.count());
        QVector<QVector<Velocity3D>> all_vels(molnums.count());
        QVector<QVector<Force3D>> all_forces(molnums.count());

        auto all_coords_data = all_coords.data();
        auto all_vels_data = all_vels.data();
        auto all_forces_data = all_forces.data();

        const auto coords_property = map["coordinates"];
        const auto vels_property = map["velocity"];
        const auto forces_property = map["force"];

        if (should_run_in_parallel(molnums.count(), map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, molnums.count()), [&](const tbb::blocked_range<int> r)
                              {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    const auto mol = system[molnums_data[i]].molecule();

                    tbb::parallel_invoke([&]() { all_coords_data[i] = ::getCoordinates(mol, coords_property); },
                                         [&]() { all_vels_data[i] = ::getVelocities(mol, vels_property); },
                                         [&]() { all_forces_data[i] = ::getForces(mol, forces_property); });
                } });
        }
        else
        {
            for (int i = 0; i < molnums.count(); ++i)
            {
                const auto mol = system[molnums_data[i]].molecule();

                all_coords_data[i] = ::getCoordinates(mol, coords_property);
                all_vels_data[i] = ::getVelocities(mol, vels_property);
                all_forces_data[i] = ::getForces(mol, forces_property);
            }
        }

        if (::hasData(all_coords))
        {
            // check the edge case that some molecules had coordinates defined, but some didn't
            bool none_have_coords = true;
            bool all_have_coords = true;
            bool some_have_coords = false;

            for (const auto &c : all_coords)
            {
                if (c.isEmpty())
                {
                    all_have_coords = false;
                }
                else
                {
                    none_have_coords = false;
                }

                if (all_have_coords == false and none_have_coords == false)
                {
                    // ok, only some have coordinates...
                    some_have_coords = true;
                    break;
                }
            }

            if (some_have_coords)
            {
                // we need to populate the empty velocities with values of 0
                for (int i = 0; i < molnums.count(); ++i)
                {
                    auto &coords_i = all_coords_data[i];

                    if (coords_i.isEmpty())
                    {
                        const int natoms = system[molnums[i]].atoms().count();
                        coords_i = QVector<Vector>(natoms, Vector(0));
                    }
                }
            }

            coords = collapse(all_coords);
        }

        if (::hasData(all_vels))
        {
            // check the edge case that some molecules had velocities defined, but some didn't
            bool none_have_velocities = true;
            bool all_have_velocities = true;
            bool some_have_velocities = false;

            for (const auto &vels : all_vels)
            {
                if (vels.isEmpty())
                {
                    all_have_velocities = false;
                }
                else
                {
                    none_have_velocities = false;
                }

                if (all_have_velocities == false and none_have_velocities == false)
                {
                    // ok, only some have velocities...
                    some_have_velocities = true;
                    break;
                }
            }

            if (some_have_velocities)
            {
                // we need to populate the empty velocities with values of 0
                for (int i = 0; i < molnums.count(); ++i)
                {
                    auto &vels_i = all_vels_data[i];

                    if (vels_i.isEmpty())
                    {
                        const int natoms = system[molnums[i]].atoms().count();
                        vels_i = QVector<Velocity3D>(natoms, Velocity3D(0));
                    }
                }
            }

            vels = collapse(all_vels);
        }

        if (::hasData(all_forces))
        {
            // check the edge case that some molecules had forces defined, but some didn't
            bool none_have_forces = true;
            bool all_have_forces = true;
            bool some_have_forces = false;

            for (const auto &frcs : all_forces)
            {
                if (frcs.isEmpty())
                {
                    all_have_forces = false;
                }
                else
                {
                    none_have_forces = false;
                }

                if (all_have_forces == false and none_have_forces == false)
                {
                    // ok, only some have forces...
                    some_have_forces = true;
                    break;
                }
            }

            if (some_have_forces)
            {
                // we need to populate the empty velocities with values of 0
                for (int i = 0; i < molnums.count(); ++i)
                {
                    auto &frcs_i = all_forces_data[i];

                    if (frcs_i.isEmpty())
                    {
                        const int natoms = system[molnums[i]].atoms().count();
                        frcs_i = QVector<Force3D>(natoms, Force3D(0));
                    }
                }
            }

            frcs = collapse(all_forces);
        }
    }

    SireVol::SpacePtr space;

    const auto space_property = map["space"];
    if (system.containsProperty(space_property))
    {
        space = system.property(space_property).asA<SireVol::Space>();
    }

    auto time = get_time_from_system(system, map["time"]);

    return SireMol::Frame(coords, vels, frcs, space.read(), time);
}

void MoleculeParser::copyFromFrame(const Frame &frame, System &system, const PropertyMap &map) const
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

    const Vector *coords_array = 0;

    if (frame.hasCoordinates())
        coords_array = frame.coordinates().constData();

    const PropertyName coords_property = map["coordinates"];

    const Velocity3D *vels_array = 0;

    if (frame.hasVelocities())
        vels_array = frame.velocities().constData();

    const PropertyName vels_property = map["velocity"];

    const Force3D *frcs_array = 0;

    if (frame.hasForces())
        frcs_array = frame.forces().constData();

    const PropertyName frcs_property = map["force"];

    bool should_make_whole = false;

    const Space &space = frame.space();

    if (space.isPeriodic() and map.specified("make_whole"))
    {
        should_make_whole = map["make_whole"].value().asABoolean();
    }

    auto add_data = [&](int i)
    {
        const int atom_start_idx = atom_pointers.constData()[i];
        auto mol = system[MolIdx(i)].molecule().edit();
        const auto molinfo = mol.data().info();

        if (coords_array != 0)
        {
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

            if (should_make_whole)
                coords = space.makeWhole(coords);

            mol.setProperty(coords_property, AtomCoords(CoordGroupArray(coords)));
        }

        if (vels_array != 0)
        {
            auto vels = AtomVelocities(molinfo);

            for (int j = 0; j < mol.nAtoms(); ++j)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(j));

                const int atom_idx = atom_start_idx + j;

                vels.set(cgatomidx, vels_array[atom_idx]);
            }

            mol.setProperty(vels_property, vels).commit();
        }

        if (frcs_array != 0)
        {
            auto forces = AtomForces(molinfo);

            for (int j = 0; j < mol.nAtoms(); ++j)
            {
                auto cgatomidx = molinfo.cgAtomIdx(AtomIdx(j));

                const int atom_idx = atom_start_idx + j;

                forces.set(cgatomidx, frcs_array[atom_idx]);
            }

            mol.setProperty(frcs_property, forces).commit();
        }

        mols_array[i] = mol.commit();
    };

    if (usesParallel())
    {
        tbb::parallel_for(tbb::blocked_range<int>(0, nmols), [&](tbb::blocked_range<int> r)
                          {
            for (int i = r.begin(); i < r.end(); ++i)
            {
                add_data(i);
            } });
    }
    else
    {
        for (int i = 0; i < nmols; ++i)
        {
            add_data(i);
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

    PropertyName time_property = map["time"];
    if (time_property.hasValue())
    {
        system.setProperty("time", time_property.value());
    }
    else
    {
        system.setProperty(time_property.source(), GeneralUnitProperty(frame.time()));
    }
}

/** Remove any comment lines (those that start with 'comment_flag')
 *  from the file. This should make parsing easier
 */
void MoleculeParser::removeCommentLines(const QString &comment_flag)
{
    bool has_comments = false;

    for (int i = 0; i < lnes.count(); ++i)
    {
        const auto &line = lnes.constData()[i];

        if (line.startsWith(comment_flag))
        {
            has_comments = true;
            break;
        }
    }

    if (has_comments)
    {
        QMutableVectorIterator<QString> it(lnes);

        while (it.hasNext())
        {
            const auto &line = it.next();
            if (line.startsWith(comment_flag))
            {
                it.remove();
            }
        }
    }
}

/** Function used by derived classes to set the lines */
void MoleculeParser::setLines(const QVector<QString> &lines)
{
    lnes = lines;
}

/** Functions used by derived classes to set the filename */
void MoleculeParser::setFilename(const QString &filename)
{
    fname = QFileInfo(filename).absoluteFilePath();
}

const char *MoleculeParser::typeName()
{
    return "SireIO::MoleculeParser";
}

MoleculeParser &MoleculeParser::operator=(const MoleculeParser &other)
{
    if (this != &other)
    {
        fname = other.fname;
        lnes = other.lnes;
        saved_system = other.saved_system;
        loaded_order = other.loaded_order;
        frames_to_write = other.frames_to_write;
        propmap = other.propmap;
        scr = other.scr;
        run_parallel = other.run_parallel;
        Property::operator=(other);
    }

    return *this;
}

bool MoleculeParser::operator==(const MoleculeParser &other) const
{
    return fname == other.fname and lnes == other.lnes and scr == other.scr and run_parallel == other.run_parallel and
           saved_system == other.saved_system and frames_to_write == other.frames_to_write and propmap == other.propmap and
           Property::operator==(other);
}

bool MoleculeParser::operator!=(const MoleculeParser &other) const
{
    return not MoleculeParser::operator==(other);
}

/** Return the name of the file that was parsed */
QString MoleculeParser::filename() const
{
    return fname;
}

/** Return whether or not this parser is broken */
bool MoleculeParser::isBroken() const
{
    return false;
}

/** Return whether or not this parser is a topology parser */
bool MoleculeParser::isTopology() const
{
    return false;
}

/** Return whether or not this parser is a frame parser */
bool MoleculeParser::isFrame() const
{
    return false;
}

/** Return whether or not this parser is a supplementary parser */
bool MoleculeParser::isSupplementary() const
{
    return false;
}

/** Return the error report, if this parser is broken. If it isn't,
 *  then an empty string is returned. */
QString MoleculeParser::errorReport() const
{
    return QString();
}

/** Return any warnings that were generated when loading data
 *  using this parser
 */
QStringList MoleculeParser::warnings() const
{
    return QStringList();
}

/** Return whether there were any warnings when loading the file
 *  using this parser */
bool MoleculeParser::hasWarnings() const
{
    return not this->warnings().isEmpty();
}

/** Internal function to set the atom order that was used to load
 *  the atoms from the topology file - this is only set if the
 *  atoms are not in their expected order
 */
void MoleculeParser::setLoadedOrder(const QVector<qint64> &order)
{
    loaded_order = order;
    this->reorderLoadedFrame();
}

/** Internal function to get the order in which atoms were loaded.
 *  This is empty if the atoms were loaded in the expected order
 */
const QVector<qint64> &MoleculeParser::getLoadedOrder() const
{
    return loaded_order;
}

/** Internal function used to reorder the passed frame based on the
 *  loaded order (if this differs to the expected loaded order)
 */
Frame MoleculeParser::reorderFrame(const Frame &frame) const
{
    if (loaded_order.isEmpty())
        return frame;
    else
        return frame.reorder(loaded_order);
}

/** Return the number of trajectory frames contained in this parser.
 *  Trajectory frames contain coordinates and/or velocities and/or
 *  forces data. It is possible for a parser to have zero frames,
 *  e.g. if it only contains topology information.
 */
int MoleculeParser::nFrames() const
{
    return 0;
}

/** Implement this function to be signalled when you need to
 *  reorder the loaded frame
 */
void MoleculeParser::reorderLoadedFrame()
{
}

/** Return the ith trajectory frame from this parser. Note that
 *  some parsers may have to re-read the file, so this may fail
 *  if the filename changes since the last time this parser
 *  was used
 */
SireMol::Frame MoleculeParser::getFrame(int i) const
{
    // this will raise an exception as we will only
    // be calling this function if the parser has zero frames
    i = SireID::Index(i).map(0);

    return SireMol::Frame();
}

/** Enable code to parse files in parallel */
void MoleculeParser::enableParallel()
{
    run_parallel = true;
}

/** Disable code to parse files in parallel - parsing will happen in serial */
void MoleculeParser::disableParallel()
{
    run_parallel = false;
}

/** Set whether or not to parse files in parallel or serial */
void MoleculeParser::setUseParallel(bool on)
{
    run_parallel = on;
}

/** Extract and return a FFDetail forcefield that is compatible with all of the
    molecules in this system, using the passed property map to find the property.
    Note that this will raise an incompatible_error exception if there is no
    forcefield that adequately covers all of the molecules */
PropertyPtr MoleculeParser::getForceField(const System &system, const PropertyMap &map) const
{
    const auto ffprop = map["forcefield"];

    // make sure that the user is not telling us a specific forcefield to use
    if (ffprop.hasValue())
    {
        if (not ffprop.value().isA<FFDetail>())
            throw SireError::incompatible_error(
                QObject::tr("Cannot convert the passed mapped value '%1' to an object of type FFDetail. "
                            "If you want to specify the focefield it must be an object derived from FFDetail.")
                    .arg(ffprop.value().toString()),
                CODELOC);

        return ffprop.value();
    }

    const auto molnums = system.molNums().toVector();

    if (molnums.isEmpty())
        return SireMM::MMDetail();

    PropertyPtr ffield;
    QStringList errors;

    for (int i = 0; i < molnums.count(); ++i)
    {
        const auto mol = system[molnums[i]].molecule();
        PropertyPtr molff;

        if (mol.hasProperty(ffprop))
        {
            const auto &p = mol.property(ffprop);

            if (p.isA<FFDetail>())
                molff = p;
        }

        if (molff.isNull())
        {
            errors.append(QObject::tr("Molecule '%1' does not have a valid 'forcefield' "
                                      "property. Please make sure that it has a property called '%2' that is "
                                      "derived from FFDetail")
                              .arg(mol.toString())
                              .arg(ffprop.source()));
        }
        else if (ffield.isNull())
        {
            // this is the first valid forcefield
            ffield = molff;
        }
        else if (not ffield.read().asA<FFDetail>().isCompatibleWith(molff.read().asA<FFDetail>()))
        {
            // incompatible forcefields!
            errors.append(QObject::tr("The forcefield for molecule '%1' is not compatible "
                                      "with that for other molecules.\n%2\nversus\n%3.")
                              .arg(mol.toString())
                              .arg(molff.read().toString())
                              .arg(ffield.read().toString()));
        }
    }

    if (not errors.isEmpty())
        throw SireError::incompatible_error(
            QObject::tr("There were some problems "
                        "extracting a valid forcefield object from all of the molecules.\n\n%1")
                .arg(errors.join("\n")),
            CODELOC);

    // we should have a valid FFDetail now...
    return ffield.read().asA<FFDetail>();
}

/** Create a directory for the passed file (if it doesn't already exist) */
void MoleculeParser::createDirectoryForFile(const QString &filename) const
{
    const auto fileinfo = QFileInfo(filename);

    const auto basename = fileinfo.dir().filePath(fileinfo.baseName());
    const auto basedir = fileinfo.absoluteDir();

    if (basedir.exists())
        return;

    // need to use a mutex to minimise chance of a race condition

    static QMutex create_dir_mutex;

    QMutexLocker lkr(&create_dir_mutex);

    if (not basedir.exists())
    {
        // create this directory
        if (not basedir.mkpath("."))
        {
            throw SireError::file_error(QObject::tr(
                                            "Could not create the directory '%1' into which the file "
                                            "'%2' will be written. Make sure there is enough space "
                                            "and you have the right permissions.")
                                            .arg(basedir.absolutePath())
                                            .arg(fileinfo.fileName()),
                                        CODELOC);
        }
    }
}

void _check_and_remove_frame_dir(const QFileInfo &fileinfo)
{
    // we will only remove directories that only contain 'frame_XXX' files
    auto dir = QDir(fileinfo.absoluteFilePath());

    bool ok_to_delete = true;

    for (const auto &subfile : dir.entryInfoList(QDir::Files))
    {
        if (not subfile.fileName().startsWith("frame_"))
        {
            ok_to_delete = false;
            break;
        }
    }

    if (not ok_to_delete)
        return;

    // now try to remove the files...
    for (const auto &subfile : dir.entryInfoList(QDir::Files))
    {
        if (subfile.fileName().startsWith("frame_"))
        {
            dir.remove(subfile.fileName());
        }
    }

    dir = QDir(fileinfo.absolutePath());
    dir.rmdir(fileinfo.fileName());
}

/** Write the parsed data back to the file called 'filename'. This will
    overwrite the file if it exists already, so be careful! Note that
    this will write this to multiple files if trajectory writing
    is requested and this is a frame parser
*/
QStringList MoleculeParser::writeToFile(const QString &filename) const
{
    if (lnes.isEmpty())
        return QStringList();

    if (not this->isTextFile())
        throw SireError::program_bug(
            QObject::tr("Dear programmer - please override the MoleculeParser::writeToFile function "
                        "to work with your binary file format. Text writing is not supported "
                        "for the parser %1.")
                .arg(this->what()),
            CODELOC);

    auto gil = SireBase::release_gil();

    QStringList written_files;

    // Check whether list of trajectory frames have been specified via the map.
    // If so, then the name will be prefixed with "frames_names:"
    QStringList frame_names;
    if (filename.contains("frames_names:"))
    {
        // Split the part of the string following "frames_names:" by commas to get
        // individual frame names.
        const QString frames_names_str = filename.section("frames_names:", 1);
        frame_names = frames_names_str.split(",", Qt::SkipEmptyParts);

        // Now get the absolute path. (bit before "frames_names:")
        const QString abs_path = filename.section("frames_names:", 0, 0);

        // Append the absolute path to each frame name.
        for (auto &frame_name : frame_names)
        {
            frame_name = QDir(abs_path).filePath(frame_name.trimmed());
            frame_name = QFileInfo(frame_name).absoluteFilePath();
        }

        if (frame_names.size() != frames_to_write.size())
        {
            throw SireError::program_bug(
                QObject::tr("The number of frame names provided (%1) does not match the number of frames to write (%2).")
                    .arg(frame_names.size())
                    .arg(frames_to_write.size()),
                CODELOC);
        }
    }
    else
    {
        createDirectoryForFile(filename);
    }

    if (this->writingTrajectory() and this->isFrame())
    {
        System s = this->saved_system;

        // construct a copy of our property map, but without
        // including the list of frames to write
        auto d = this->propmap.toDict();
        d.remove("frames_to_write");
        const PropertyMap m(d);

        const int nframes = frames_to_write.count();

        ProgressBar bar(QString("Save %1").arg(this->formatName()), nframes);
        bar.setSpeedUnit("frames / s");

        bar = bar.enter();

        // find the largest frame so that we can pad with leading zeroes
        const int largest_frame = *std::max_element(frames_to_write.constBegin(),
                                                    frames_to_write.constEnd());

        const int padding = QString::number(largest_frame).length();

        QDir framedir;
        QFileInfo fileinfo;
        bool write_to_frame_dir = false;

        if (frame_names.isEmpty())
        {
            fileinfo = QFileInfo(filename);
            write_to_frame_dir = true;

            // in this case, we are going to create a directory named after the
            // filename, into which all the frames will be written
            fileinfo = QFileInfo(filename);

            if (fileinfo.exists())
            {
                // by default, we support overwriting of files - so remove this if it is a file
                if (fileinfo.isDir())
                {
                    _check_and_remove_frame_dir(fileinfo);

                    // rebuild so it is not cached
                    fileinfo = QFileInfo(filename);

                    if (fileinfo.exists())
                        throw SireError::file_error(QObject::tr(
                                                        "Could not write the trajectory for file '%1' as there "
                                                        "is already a directory with this name that contains "
                                                        "files that don't appear to have been created by sire.")
                                                        .arg(filename),
                                                    CODELOC);
                }
                else
                {
                    auto dir = fileinfo.absoluteDir();

                    if (not dir.remove(fileinfo.fileName()))
                        throw SireError::file_error(QObject::tr(
                                                        "Could not write the trajectory for file '%1' as "
                                                        "we don't have permission to remove the existing "
                                                        "file with this name.")
                                                        .arg(filename),
                                                    CODELOC);
                }
            }

            framedir = QDir(fileinfo.absoluteFilePath());

            if (not framedir.mkpath("."))
                throw SireError::file_error(QObject::tr(
                                                "Could not create the directory into which to write the "
                                                "trajectory for '%1'. Check that there is enough space "
                                                "and you have the correct permissions.")
                                                .arg(framedir.absolutePath()),
                                            CODELOC);
        }

        const auto suffix = fileinfo.completeSuffix();

        const auto time_property = m["time"];

        if (this->usesParallel())
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, nframes), [&](tbb::blocked_range<int> r)
                              {
                System thread_s(s);

                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto frame = frames_to_write[i];

                    // load the specified frame
                    thread_s.loadFrame(frame, m);

                    auto time = get_time_from_system(thread_s, time_property).to(picosecond);

                    if (time < 0)
                        time = 0;

                    // construct a copy of this parser for this frame
                    auto parser = this->construct(thread_s, m);

                    QString frame_filename;

                    if (write_to_frame_dir)
                    {
                        // now write it to the file, numbered by the frame number and time
                        QDir framedir(fileinfo.absoluteFilePath());

                        frame_filename = framedir.filePath(
                            "frame_" +
                            QString::number(i).rightJustified(padding, '0') +
                            "_" +
                            QString::number(time).replace(".", "-") +
                            "." + suffix);
                    }
                    else
                    {
                        // write to the specified frame name
                        frame_filename = frame_names[i];
                        createDirectoryForFile(frame_filename);
                    }

                    parser.read().writeToFile(frame_filename);

                    bar.tick();
                } });
        }
        else
        {
            for (int i = 0; i < nframes; ++i)
            {
                const auto frame = frames_to_write[i];

                // load the specified frame
                s.loadFrame(frame, m);

                auto time = get_time_from_system(s, time_property).to(picosecond);

                // construct a copy of this parser for this frame
                auto parser = this->construct(s, m);

                QString frame_filename;

                if (write_to_frame_dir)
                {
                    // now write it to the file, numbered by the frame number
                    QDir framedir(fileinfo.absoluteFilePath());

                    frame_filename = framedir.filePath(
                        "frame_" +
                        QString::number(i).rightJustified(padding, '0') +
                        "_" +
                        QString::number(time).replace(".", "-") +
                        "." + suffix);
                }
                else
                {
                    // write to the specified frame name
                    frame_filename = frame_names[i];
                    createDirectoryForFile(frame_filename);
                }

                parser.read().writeToFile(frame_filename);

                bar.tick();
            }
        }

        bar.success();

        // only return the directory name, as we can handle
        // reading all the frames contained therein
        if (write_to_frame_dir)
        {
            written_files.append(framedir.absolutePath());
        }
        else
        {
            written_files.append(frame_names);
        }
    }
    else
    {
        QFile f(filename);

        if (not f.open(QIODevice::WriteOnly | QIODevice::Text))
        {
            throw SireError::file_error(f, CODELOC);
        }

        QTextStream ts(&f);

        for (const QString &line : lnes)
        {
            ts << line << '\n';
        }

        f.close();

        written_files.append(filename);
    }

    return written_files;
}

/** Internal function that actually tries to parse the supplied file with name
    'filename'. This will try to find a parser based on suffix, but if that fails,
    it will try all parsers. It will choose the parser that doesn't raise an error
    that scores highest. */
MoleculeParserPtr MoleculeParser::_pvt_parse(const QString &filename, const PropertyMap &map)
{
    bool disable_supplementary = false;

    if (map.specified("DISABLE_SUPPLEMENTARY"))
    {
        disable_supplementary = map["DISABLE_SUPPLEMENTARY"].value().asA<BooleanProperty>().value();
    }

    QFileInfo info(filename);

    if (not info.isReadable())
    {
        throw SireError::file_error(QObject::tr("There is nothing readable called '%1'.").arg(filename), CODELOC);
    }

    if (info.isDir())
    {
        // this is a directory containing (potentially) multiple
        // files representing multiple frames
        return MoleculeParserPtr(FileTrajectoryParser(filename, map));
    }

    // try to find the right parser based on the suffix
    QString suffix = info.suffix().toLower();

    QStringList errors;
    QStringList suffix_errors;
    QStringList recognised_suffixes;
    QMap<float, MoleculeParserPtr> parsers;

    if (not suffix.isEmpty())
    {
        for (auto factory : getParserFactory()->factoriesForSuffix(suffix, disable_supplementary))
        {
            recognised_suffixes.append(factory.preferredSuffix());

            try
            {
                const auto parser = factory.construct(info.absoluteFilePath(), map);

                if (parser.read().score() <= 0)
                {
                    suffix_errors.append(QObject::tr("*-- Failed to parse '%1' with parser '%2'.\n"
                                                     "The file is not recognised as being of the required format.")
                                             .arg(filename)
                                             .arg(factory.formatName()));
                }
                else
                {
                    parsers.insert(parser.read().score(), parser);
                }
            }
            catch (const SireError::exception &e)
            {
                suffix_errors.append(QObject::tr("*-- Failed to parse '%1' with parser '%2'.\n%3")
                                         .arg(filename)
                                         .arg(factory.formatName())
                                         .arg(e.error()));
            }
        }

        if (not parsers.isEmpty())
        {
            // return the parser with the highest score
            return parsers.last();
        }
    }

    // none of the tested parsers worked, so let's now try all of the parsers
    for (auto factory : getParserFactory()->factoriesExcludingSuffix(suffix, disable_supplementary))
    {
        try
        {
            const auto parser = factory.construct(filename, map);

            if (parser.read().score() <= 0)
            {
                errors.append(QObject::tr("Failed to parse '%1' with parser '%2' "
                                          "as this file is not recognised as being of the required format.")
                                  .arg(filename)
                                  .arg(factory.formatName()));
            }
            else
            {
                parsers.insert(parser.read().score(), parser);
            }
        }
        catch (const SireError::exception &e)
        {
            errors.append(QObject::tr("Failed to parse '%1' with parser '%2'\n%3")
                              .arg(filename)
                              .arg(factory.formatName())
                              .arg(e.error()));
        }
    }

    if (not parsers.isEmpty())
    {
        return parsers.last();
    }
    else
    {
        if (not recognised_suffixes.isEmpty())
        {
            return MoleculeParserPtr(BrokenParser(filename, recognised_suffixes.join(","), suffix_errors));
        }
        else
        {
            return MoleculeParserPtr(BrokenParser(filename, suffix_errors + errors));
        }

        return MoleculeParserPtr();
    }
}

/** Parse the passed system, returning the resulting Parser. You must
 *  specify the parser that you want to use
 */
MoleculeParserPtr MoleculeParser::parse(const System &system, const QString &format, const PropertyMap &map)
{
    const auto factories = getParserFactory()->getFactories({format});

    if (factories.count() == 0)
        return MoleculeParserPtr(new BrokenParser());

    return factories[0].construct(system, map);
}

/** Parse the passed file, returning the resulting Parser. This employs a lot
    of magic to automatically work out the format of the file and whether or
    not this is parseable by Sire... This raises an exception if the file
    cannot be recognised, or if there is an error in parsing. */
MoleculeParserPtr MoleculeParser::parse(const QString &filename, const PropertyMap &map)
{
    MoleculeParserPtr parser;

    try
    {
        parser = MoleculeParser::_pvt_parse(filename, map);
        getFileCache()->clear();
    }
    catch (...)
    {
        getFileCache()->clear();
        throw;
    }

    return parser;
}

/** Parse the passed set of files, returning the resulting Parsers */
QList<MoleculeParserPtr> MoleculeParser::parse(const QStringList &filenames, const PropertyMap &map)
{
    QList<MoleculeParserPtr> result;

    try
    {
        if (filenames.count() == 1)
        {
            result.append(MoleculeParser::_pvt_parse(filenames[0], map));
        }
        else
        {
            QVector<MoleculeParserPtr> parsers(filenames.count());

            bool run_parallel = true;

            if (map["parallel"].hasValue())
            {
                run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
            }

            if (run_parallel)
            {
                // parse the files in parallel - we use a grain size of 1
                // as each file can be pretty big, and there won't be many of them
                tbb::parallel_for(
                    tbb::blocked_range<int>(0, filenames.count(), 1),
                    [&](tbb::blocked_range<int> r)
                    {
                        for (int i = r.begin(); i < r.end(); ++i)
                        {
                            parsers[i] = MoleculeParser::_pvt_parse(filenames[i], map);
                        }
                    },
                    tbb::simple_partitioner());
            }
            else
            {
                for (int i = 0; i < filenames.count(); ++i)
                {
                    parsers[i] = MoleculeParser::_pvt_parse(filenames[i], map);
                }
            }

            result = parsers.toList();
        }

        getFileCache()->clear();
    }
    catch (...)
    {
        getFileCache()->clear();
        throw;
    }

    return result;
}

/** Read the passed file called 'filename', returning the System contained therein */
System MoleculeParser::read(const QString &filename, const PropertyMap &map)
{
    MoleculeParserPtr parser = MoleculeParser::parse(filename, map);
    auto system = parser.read().toSystem(map);

    if (system.name().isEmpty())
    {
        system.setName(QFileInfo(filename).baseName());
    }

    return system;
}

/** Read the two passed files, returning the System contained therein. The two
    files must refer to the same System, i.e. they could be a parameter + coordinate file */
System MoleculeParser::read(const QString &file1, const QString &file2, const PropertyMap &map)
{
    MoleculeParserPtr parser1, parser2;

    bool run_parallel = true;

    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (run_parallel)
    {
        tbb::parallel_invoke([&]()
                             { parser1 = MoleculeParser::parse(file1, map); },
                             [&]()
                             { parser2 = MoleculeParser::parse(file2, map); });
    }
    else
    {
        parser1 = MoleculeParser::parse(file1, map);
        parser2 = MoleculeParser::parse(file2, map);
    }

    auto system = parser1.read().toSystem(parser2.read(), map);

    if (system.name().isEmpty())
    {
        auto p1 = QFileInfo(file1).baseName();
        auto p2 = QFileInfo(file2).baseName();

        if (p1 == p2)
        {
            system.setName(p1);
        }
        else
        {
            system.setName(QString("%1:%2").arg(p1).arg(p2));
        }
    }

    return system;
}

/** Read the files with passed filenames, returning the System contained therein.
    Note that all of the files must be connected to the same system
    (i.e. it could be the Amber Parm and Rst file) */
System MoleculeParser::read(const QStringList &filenames, const PropertyMap &map)
{
    QList<MoleculeParserPtr> parsers = MoleculeParser::parse(filenames, map);

    if (parsers.isEmpty())
        return System();

    MoleculeParserPtr parser = parsers.takeFirst();

    auto system = parser.read().toSystem(parsers, map);

    if (system.name().isEmpty())
    {
        QSet<QString> parts;
        for (const auto &filename : filenames)
        {
            parts.insert(QFileInfo(filename).baseName());
        }

        // QStringList names(parts.constBegin(), parts.constEnd());
        QStringList names = parts.values();

        system.setName(names.join(":"));
    }

    return system;
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &filename, const PropertyMap &map)
{
    return MoleculeParser::read(filename, map);
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &file1, const QString &file2, const PropertyMap &map)
{
    return MoleculeParser::read(file1, file2, map);
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QStringList &filenames, const PropertyMap &map)
{
    return MoleculeParser::read(filenames, map);
}

/** Parse the passed file, returning the resulting Parser. This employs a lot
    of magic to automatically work out the format of the file and whether or
    not this is parseable by Sire... This raises an exception if the file
    cannot be recognised, or if there is an error in parsing. */
MoleculeParserPtr MoleculeParser::parse(const QString &filename)
{
    return parse(filename, PropertyMap());
}

/** Parse the passed set of files, returning the resulting Parsers */
QList<MoleculeParserPtr> MoleculeParser::parse(const QStringList &filenames)
{
    return parse(filenames, PropertyMap());
}

/** Read the passed file called 'filename', returning the System contained therein */
System MoleculeParser::read(const QString &filename)
{
    return read(filename, PropertyMap());
}

/** Read the two passed files, returning the System contained therein. The two
    files must refer to the same System, i.e. they could be a parameter + coordinate file */
System MoleculeParser::read(const QString &file1, const QString &file2)
{
    return read(file1, file2, PropertyMap());
}

/** Read the files with passed filenames, returning the System contained therein.
    Note that all of the files must be connected to the same system
    (i.e. it could be the Amber Parm and Rst file) */
System MoleculeParser::read(const QStringList &filenames)
{
    return read(filenames, PropertyMap());
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &filename)
{
    return load(filename, PropertyMap());
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QString &file1, const QString &file2)
{
    return load(file1, file2, PropertyMap());
}

/** Synonym for MoleculeParser::read */
System MoleculeParser::load(const QStringList &filenames)
{
    return load(filenames, PropertyMap());
}

/** Return the suffix (or suffixes) given to files that support this format.
    The first suffix is the preferred on to use */
QStringList MoleculeParser::formatSuffix() const
{
    // just use a lower-case version of the format name
    return QStringList(this->formatName().toLower());
}

/** This returns a human readable set of lines describing the formats supported
    by MoleculeParser. Each line is formatted as "extension : description" where
    extension is the unique extension of the file used by MoleculeParser, and
    description is a description of the file format */
QString MoleculeParser::supportedFormats()
{
    return getParserFactory()->supportedFormats();
}

QStringList pvt_write(System system, const QStringList &filenames, const QStringList &fileformats,
                      const PropertyMap &map)
{
    if (filenames.count() != fileformats.count())
    {
        throw SireError::program_bug(QObject::tr("Disagreement of the number of files... %1 vs %2")
                                         .arg(filenames.count())
                                         .arg(fileformats.count()),
                                     CODELOC);
    }

    // release the GIL here so that progress bars can be displayed
    auto handle = SireBase::release_gil();

    // check that the files can be written
    QVector<QFileInfo> fileinfos(filenames.count());

    QStringList errors;

    for (int i = 0; i < filenames.count(); ++i)
    {
        fileinfos[i] = QFileInfo(filenames[i]);

        if (fileinfos[i].exists())
        {
            if (fileinfos[i].isDir())
            {
                _check_and_remove_frame_dir(fileinfos[i]);

                // rebuild this, so this isn't cached
                fileinfos[i] = QFileInfo(filenames[i]);

                if (fileinfos[i].exists())
                    errors.append(QObject::tr("The file %1 is actually a directory, and not writable!")
                                      .arg(fileinfos[i].absoluteFilePath()));
            }
            else if (not fileinfos[i].isWritable())
            {
                errors.append(
                    QObject::tr("The file %1 exists and is not writable!").arg(fileinfos[i].absoluteFilePath()));
            }
        }
    }

    if (not errors.isEmpty())
    {
        throw SireError::io_error(
            QObject::tr("Cannot write the files as the following errors occurred:\n%1").arg(errors.join("\n\n")),
            CODELOC);
    }

    // now get all of the parsers
    const auto factories = getParserFactory()->getFactories(fileformats);

    QVector<QStringList> written_files(filenames.count());

    // should we write the files in parallel?
    bool run_parallel = true;

    if (map["parallel"].hasValue())
    {
        run_parallel = map["parallel"].value().asA<BooleanProperty>().value();
    }

    if (run_parallel)
    {
        tbb::spin_mutex error_mutex;

        tbb::parallel_for(
            tbb::blocked_range<int>(0, filenames.count(), 1),
            [&](const tbb::blocked_range<int> &r)
            {
                for (int i = r.begin(); i < r.end(); ++i)
                {
                    const auto filename = fileinfos[i].absoluteFilePath();

                    try
                    {
                        written_files[i] = factories[i].construct(system, map).read().writeToFile(filename);
                    }
                    catch (const SireError::exception &e)
                    {
                        tbb::spin_mutex::scoped_lock locker(error_mutex);
                        errors.append(QObject::tr("Failed to write the file '%1' using the parser "
                                                  "for fileformat '%2'. Errors reported by the parser are:\n%3")
                                          .arg(filename)
                                          .arg(factories[i].formatName())
                                          .arg(e.error()));
                    }
                }
            },
            tbb::simple_partitioner());
    }
    else
    {
        for (int i = 0; i < filenames.count(); ++i)
        {
            const auto filename = fileinfos[i].absoluteFilePath();

            try
            {
                written_files[i] = factories[i].construct(system, map).read().writeToFile(filename);
            }
            catch (const SireError::exception &e)
            {
                errors.append(QObject::tr("Failed to write the file '%1' using the parser "
                                          "for fileformat '%2'. Errors reported by the parser are:\n%3")
                                  .arg(filename)
                                  .arg(factories[i].formatName())
                                  .arg(e.error()));
            }
        }
    }

    if (not errors.isEmpty())
    {
        throw SireError::io_error(QObject::tr("Cannot write the (perhaps some of the ) files "
                                              "as the following errors occurred:\n%1")
                                      .arg(errors.join("\n\n")),
                                  CODELOC);
    }

    QStringList w;

    for (const auto &wf : written_files)
    {
        w += wf;
    }

    return w;
}

/** Save the passed System to the file called 'filename'. First, the 'fileformat'
    property is looked at in 'map'. This is used to set the format(s) of
    the files that are written (comma-separated list).

    If this does not exist, then the extension of the
    file is used to work out which format to use. If no extension is given,
    then the System will be queried to find out its preferred format (normally
    the format it was loaded with), via its 'fileformat' property
    (again, comma separated list).

    If their preferred format results in multiple files, then
    multiple files will be written. This returns the full pathnames to
    all of the files that are written
*/
QStringList MoleculeParser::write(const System &system, const QString &filename, const PropertyMap &map)
{
    if (filename.isEmpty() or QFileInfo(filename).baseName().isEmpty())
    {
        throw SireError::io_error(
            QObject::tr("You must supply a valid filename. This '%1' is not sufficient.").arg(filename), CODELOC);
    }

    // build a list of filenames with their associated fileformats
    QStringList filenames;
    QStringList fileformats;

    if (map.specified("fileformat"))
    {
        const auto format_property = map["fileformat"];

        if (format_property.hasSource())
        {
            fileformats = format_property.source().split(",");
        }
        else
        {
            try
            {
                fileformats = format_property.value().asA<StringProperty>().toString().split(",");
            }
            catch (...)
            {
                fileformats.append(format_property.value().asA<MoleculeParser>().formatName());
            }
        }

        const auto fileinfo = QFileInfo(filename);
        QString basename = fileinfo.absoluteDir().absoluteFilePath(fileinfo.completeBaseName());

        for (const auto &format : fileformats)
        {
            filenames.append(QString("%1.%2").arg(basename, format.toLower()));
        }
    }
    else
    {
        QString extension = QFileInfo(filename).completeSuffix();

        if (extension.isEmpty())
        {
            // we need to find the format from the system
            try
            {
                fileformats = system.property(map["fileformat"]).asA<StringProperty>().toString().split(",");
            }
            catch (...)
            {
                throw SireError::io_error(
                    QObject::tr("Cannot work out the fileformat to use to write the System to "
                                "file '%1'. You need to either supply the format using the "
                                "'fileformat' property in the passed map, add this to the System "
                                "as its 'fileformat' property, or pass a filename with an extension "
                                "whose fileformat can be determined. Supported fileformats are;\n%2")
                        .arg(filename)
                        .arg(MoleculeParser::supportedFormats()),
                    CODELOC);
            }

            for (const auto &format : fileformats)
            {
                filenames.append(QString("%1.%2").arg(filename).arg(format.toLower()));
            }
        }
        else
        {
            filenames.append(filename);
            fileformats.append(extension.toUpper());
        }
    }

    // Check for a frame_names property in the map, which is used to control the naming of
    // specific frames when writing trajectories.
    if (map.specified("frame_names"))
    {
        const auto frame_names_property = map["frame_names"];

        QStringList frame_names;

        if (frame_names_property.hasSource())
        {
            frame_names = frame_names_property.source().split(",");
        }
        else
        {
            try
            {
                frame_names = frame_names_property.value().asA<StringProperty>().toString().split(",");
            }
            catch (...)
            {
                throw SireError::incompatible_error(
                    QObject::tr("The 'frame_names' property must be a StringProperty containing "
                                "a comma-separated list of filenames to use for each frame when "
                                "writing trajectories."),
                    CODELOC);
            }
        }

        // add the special prefix to the first filename
        filenames[0] = "frames_names:" + frame_names.join(",");
    }

    // now we have a list of filenames and associated formats, actually
    // write the files
    return ::pvt_write(system, filenames, fileformats, map);
}

/** Extension of MoleculeParser::write which allows many filenames.
    The same rules to locate the fileformats are now used, except that now only
    the number of files written must match the number of filenames */
QStringList MoleculeParser::write(const System &system, const QStringList &files, const PropertyMap &map)
{
    if (files.isEmpty())
    {
        throw SireError::io_error(QObject::tr("You must supply a valid filename. An empty list is not sufficient!"),
                                  CODELOC);
    }
    else if (files.count() == 1)
    {
        return MoleculeParser::write(system, files[0], map);
    }

    // build a list of filenames with their associated fileformats
    QStringList filenames;
    QStringList fileformats;

    if (map.specified("fileformat"))
    {
        const auto format_property = map["fileformat"];

        if (format_property.hasSource())
        {
            fileformats = format_property.source().split(",");
        }
        else
        {
            try
            {
                fileformats = format_property.value().asA<StringProperty>().toString().split(",");
            }
            catch (...)
            {
                fileformats.append(format_property.value().asA<MoleculeParser>().formatName());
            }
        }

        if (files.count() != fileformats.count())
        {
            throw SireError::io_error(QObject::tr("You must match up the number of filenames to fileformats when "
                                                  "specifying both the filenames [%1] and fileformats [%2].")
                                          .arg(filenames.join(","))
                                          .arg(fileformats.join(",")),
                                      CODELOC);
        }

        for (int i = 0; i < fileformats.count(); ++i)
        {
            const QString filename = files[i];

            const auto fileinfo = QFileInfo(filename);

            if (filename.isEmpty() or fileinfo.completeBaseName().isEmpty())
            {
                throw SireError::io_error(
                    QObject::tr("You must supply a valid filename. This '%1' is not sufficient.").arg(filename),
                    CODELOC);
            }

            QString basename = fileinfo.absoluteDir().absoluteFilePath(fileinfo.completeBaseName());

            filenames.append(QString("%1.%2").arg(basename, fileformats[i].toLower()));
        }
    }
    else
    {
        // we may need to find the format from the system
        try
        {
            fileformats = system.property(map["fileformat"]).asA<StringProperty>().toString().split(",");
        }
        catch (...)
        {
        }

        for (int i = 0; i < files.count(); ++i)
        {
            const auto filename = files[i];

            QString extension = QFileInfo(filename).completeSuffix();

            if (extension.isEmpty())
            {
                if (i >= fileformats.count())
                {
                    throw SireError::io_error(
                        QObject::tr("Cannot work out the fileformat to use to write the System to "
                                    "file '%1'. You need to either supply the format using the "
                                    "'fileformat' property in the passed map, add this to the System "
                                    "as its 'fileformat' property, or pass a filename with an extension "
                                    "whose fileformat can be determined. Supported fileformats are;\n%2")
                            .arg(filename)
                            .arg(MoleculeParser::supportedFormats()),
                        CODELOC);
                }
                else
                {
                    filenames.append(QString("%1.%2").arg(filename).arg(fileformats[i].toLower()));
                }
            }
            else if (i >= fileformats.count())
            {
                filenames.append(filename);
                fileformats.append(extension.toUpper());
            }
            else
            {
                filenames.append(filename);
                fileformats[i] = extension.toUpper();
            }
        }
    }

    // now we have a list of filenames and associated formats, actually
    // write the files
    return ::pvt_write(system, filenames, fileformats, map);
}

/** Extension of MoleculeParser::write which allows you to specify two filenames.
    The same rules to locate the fileformats are now used, except now only two
    files are permitted to be written */
QStringList MoleculeParser::write(const System &system, const QString &file1, const QString &file2,
                                  const PropertyMap &map)
{
    QStringList filenames;
    filenames.append(file1);
    filenames.append(file2);
    return MoleculeParser::write(system, filenames, map);
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QString &filename, const PropertyMap &map)
{
    return MoleculeParser::write(system, filename, map);
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QString &file1, const QString &file2,
                                 const PropertyMap &map)
{
    return MoleculeParser::write(system, file1, file2, map);
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QStringList &filenames, const PropertyMap &map)
{
    return MoleculeParser::write(system, filenames, map);
}

/** Return the System that is constructed from the data in this parser */
System MoleculeParser::toSystem(const PropertyMap &map) const
{
    QList<MoleculeParserPtr> others;
    return this->toSystem(others, map);
}

/** Return the System that is constructed from the data in the two
    passed parsers (i.e. representing a topology and a coordinate file) */
System MoleculeParser::toSystem(const MoleculeParser &other, const PropertyMap &map) const
{
    QList<MoleculeParserPtr> others;
    others.append(other);
    return this->toSystem(others, map);
}

/** Return the System that is constructed from the information in the passed
    parsers. This will parse the information in order, meaning that data contained
    in earlier parsers may be overwritten by data from later parsers */
System MoleculeParser::toSystem(const QList<MoleculeParserPtr> &others, const PropertyMap &map) const
{
    auto parsers = this->sortParsers(others, map);

    if (parsers.count() == 0)
        return System();

    if (parsers.value("broken").count() > 0)
    {
        QStringList filenames;
        QStringList parser_errors;

        for (const auto &parser : parsers["broken"])
        {
            parser_errors.append("\n\n");
            parser_errors.append(parser.read().errorReport());
            filenames.append(parser.read().filename());
        }

        throw SireIO::parse_error(QObject::tr(
                                      "Unable to load the file: %1. Errors reported by all parsers are below.%2")
                                      .arg(filenames.join(", "))
                                      .arg(parser_errors.join("\n")),
                                  CODELOC);
    }

    if (parsers.value("topology").count() == 0)
    {
        throw SireError::program_bug(QObject::tr("Should only be here if we already have a topology!"), CODELOC);
    }

    auto topology = parsers["topology"][0];

    // Instantiate an empty system.
    System system;

    bool show_warnings = false;
    bool has_warnings = false;

    if (map["show_warnings"].hasValue())
    {
        show_warnings = map["show_warnings"].value().asA<BooleanProperty>().value();
    }

    if (parsers.value("supplementary").count() > 0)
    {
        QVector<QString> supplementary_lines;

        for (const auto &parser : parsers["supplementary"])
        {
            supplementary_lines += parser.read().lines();
        }

        try
        {
            system = topology.read().startSystem(supplementary_lines, map);
        }
        catch (const SireError::exception &)
        {
            // we couldn't load this supplementary with this topology parser.
            // We need to try other topology parsers, and accept the first one
            // that works.
            auto tops = parsers["topology"];

            bool ok = false;

            for (int i = 1; i < tops.count(); ++i)
            {
                try
                {
                    topology = tops[i];
                    system = topology.read().startSystem(supplementary_lines, map);
                    ok = true;
                }
                catch (...)
                {
                }
            }

            if (not ok)
                throw;
        }
    }
    else
    {
        system = topology.read().startSystem(map);
    }

    if (topology.read().hasWarnings())
    {
        if (show_warnings)
        {
            QTextStream cout(stdout, QIODevice::WriteOnly);

            cout << QObject::tr("\nWARNINGS encountered when parsing the topology:\n");
            cout << topology.read().warnings().join("\n");
            cout << "====\n\n";
        }

        has_warnings = true;
    }

    if (system.containsProperty(QString("loaded_atom_order")))
    {
        // get the atom order from the topology parser
        const auto loaded_order = system.property("loaded_atom_order").asA<IntegerArrayProperty>().value();
        system.removeProperty("loaded_atom_order");

        const int nats = system.nAtoms();

        if (loaded_order.count() != nats)
        {
            throw SireError::program_bug(QObject::tr("The loaded atom order does not match the number of atoms in the system!"),
                                         CODELOC);
        }

        bool in_expected_order = true;

        for (int i = 0; i < loaded_order.count(); ++i)
        {
            if (loaded_order[i] != i)
            {
                in_expected_order = false;

                if (loaded_order[i] < 0 or loaded_order[i] >= nats)
                {
                    throw SireError::program_bug(QObject::tr("The loaded atom order is not valid!"), CODELOC);
                }
            }
        }

        if (not in_expected_order)
        {
            // tell each of the frame parsers that they will need to reorder
            // their frames for this system
            for (auto &parser : parsers["frame"])
            {
                parser.edit().setLoadedOrder(loaded_order);
            }
        }
    }

    bool ignore_topology_frame = false;

    if (map.specified("ignore_topology_frame"))
    {
        ignore_topology_frame = map["ignore_topology_frame"].value().asABoolean();
    }

    if (parsers.value("frame").count() > 0)
    {
        auto frames = parsers["frame"];

        if (ignore_topology_frame)
        {
            frames.removeAll(topology);

            // we need to add frame information from the first file
            if (not frames.isEmpty())
            {
                frames[0].read().addToSystem(system, map);

                if (frames[0].read().hasWarnings())
                {
                    if (show_warnings)
                    {
                        QTextStream cout(stdout, QIODevice::WriteOnly);

                        cout << QObject::tr("\nWARNINGS encountered when adding addition system data:\n");
                        cout << frames[0].read().warnings().join("\n");
                        cout << "====\n\n";
                    }

                    has_warnings = true;
                }
            }
        }
        else if (not topology.read().isFrame())
        {
            // we need to add frame information from the first file
            frames[0].read().addToSystem(system, map);

            if (frames[0].read().hasWarnings())
            {
                if (show_warnings)
                {
                    QTextStream cout(stdout, QIODevice::WriteOnly);

                    cout << QObject::tr("\nWARNINGS encountered when adding addition system data:\n");
                    cout << frames[0].read().warnings().join("\n");
                    cout << "====\n\n";
                }

                has_warnings = true;
            }
        }

        // if there is more than one frame, then we need to store the
        // trajectory too
        int nframes = 0;

        for (const auto &frame : frames)
        {
            nframes += frame.read().nFrames();
        }

        if (nframes > 1)
        {
            QList<SireMol::TrajectoryDataPtr> trajectories;

            for (const auto &frame : frames)
            {
                trajectories.append(SireMol::TrajectoryDataPtr(new FileTrajectory(frame)));

                if (frame.read().hasWarnings())
                {
                    if (show_warnings)
                    {
                        QTextStream cout(stdout, QIODevice::WriteOnly);

                        cout << QObject::tr("\nWARNINGS encountered when adding a trajectory frame:\n");
                        cout << frame.read().warnings().join("\n");
                        cout << "====\n\n";
                    }

                    has_warnings = true;
                }
            }

            // comment out as we don't use the System trajectory property,
            // and it is confusing if it isn't updated...
            // system.setProperty("trajectory", SireMol::Trajectory(trajectories));

            // we now have to assume that the trajectories all had the atomic
            // data in the same order and that this matches the atomidx order
            // in the system...
            int start_atom = 0;

            for (int i = 0; i < system.nMolecules(); ++i)
            {
                auto mol = system[i].molecule();

                int natoms = mol.nAtoms();

                SireMol::Trajectory traj(trajectories, start_atom, natoms);

                mol = mol.edit().setProperty("trajectory", traj).commit();

                system.update(mol);

                start_atom += natoms;
            }
        }
    }

    if (has_warnings and show_warnings)
    {
        QTextStream cout(stdout, QIODevice::WriteOnly);

        cout << QObject::tr("\nSilence these warnings by passing "
                            "`show_warnings=False` when loading.\n");
    }

    return system;
}

/** Save the passed System to the file called 'filename'. First, the 'fileformat'
    property is looked at in 'map'. This is used to set the format(s) of
    the files that are written (comma-separated list).

    If this does not exist, then the extension of the
    file is used to work out which format to use. If no extension is given,
    then the System will be queried to find out its preferred format (normally
    the format it was loaded with), via its 'fileformat' property
    (again, comma separated list).

    If their preferred format results in multiple files, then
    multiple files will be written. This returns the full pathnames to
    all of the files that are written
*/
QStringList MoleculeParser::write(const System &system, const QString &filename)
{
    return write(system, filename, PropertyMap());
}

/** Extension of MoleculeParser::write which allows many filenames.
    The same rules to locate the fileformats are now used, except that now only
    the number of files written must match the number of filenames */
QStringList MoleculeParser::write(const System &system, const QStringList &files)
{
    return write(system, files, PropertyMap());
}

/** Extension of MoleculeParser::write which allows you to specify two filenames.
    The same rules to locate the fileformats are now used, except now only two
    files are permitted to be written */
QStringList MoleculeParser::write(const System &system, const QString &file1, const QString &file2)
{
    return write(system, file1, file2, PropertyMap());
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QString &filename)
{
    return save(system, filename, PropertyMap());
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QString &file1, const QString &file2)
{
    return save(system, file1, file2, PropertyMap());
}

/** Synonym of MoleculeParser::write */
QStringList MoleculeParser::save(const System &system, const QStringList &filenames)
{
    return save(system, filenames, PropertyMap());
}

/** Return the System that is constructed from the data in this parser */
System MoleculeParser::toSystem() const
{
    return toSystem(PropertyMap());
}

/** Return the System that is constructed from the data in the two
    passed parsers (i.e. representing a topology and a coordinate file) */
System MoleculeParser::toSystem(const MoleculeParser &other) const
{
    return toSystem(other, PropertyMap());
}

/** Return the System that is constructed from the information in the passed
    parsers. This will parse the information in order, meaning that data contained
    in earlier parsers may be overwritten by data from later parsers */
System MoleculeParser::toSystem(const QList<MoleculeParserPtr> &others) const
{
    return toSystem(others, PropertyMap());
}

/** Start creating a new System using the information contained in this parser,
    using the (optional) property map to name the properties */
System MoleculeParser::startSystem(const PropertyMap &map) const
{
    throw SireError::io_error(QObject::tr("There is not enough information in this parser (%1) to start "
                                          "the creation of a new System. You need to use a more detailed input file.")
                                  .arg(this->toString()),
                              CODELOC);
}

/** Start creating a new System using the information contained in this parser
    and the supplementary records contained in 'lines', using the (optional)
    property map to name the properties */
System MoleculeParser::startSystem(const QVector<QString> &lines, const PropertyMap &map) const
{
    throw SireError::io_error(QObject::tr("There is not enough information in this parser (%1) to start "
                                          "the creation of a new System. You need to use a more detailed input file.")
                                  .arg(this->toString()),
                              CODELOC);
}

/** Continue adding data to the passed System using the information contained in
    this parser, using the (optional) property map to name the properties */
void MoleculeParser::addToSystem(System &system, const PropertyMap &map) const
{
    throw SireError::io_error(QObject::tr("This parser (%1) cannot be used to add additional information to a "
                                          "System. It can only be used to create a new System from scratch.")
                                  .arg(this->toString()),
                              CODELOC);
}

/** Sort the parsers into different categories (identified by name)
 *
 *  topology : contains molecular topology (and optionally also frame information)
 *             so can be used to construct a System, and, if it has frame information,
 *             also to specify the coordinates, velocities etc
 *
 *  frame : contains frame information (coordinates and/or velocities and/or
 *          forces, with optionally a space and time)
 *
 *  supplementary : supplementary files used to give extra information.
 *                  These don't have standard molecular information
 *
 *  broken : contains all of the BrokenParser objects for files that
 *           could not be parsed
 */
QHash<QString, QList<MoleculeParserPtr>> MoleculeParser::sortParsers(const QList<MoleculeParserPtr> &parsers,
                                                                     const PropertyMap &map) const
{
    QHash<QString, QList<MoleculeParserPtr>> ret;

    // The topology parsers - we should end up with one...
    QList<MoleculeParserPtr> topology;

    // The frame parsers - we can have as many as specified...
    QList<MoleculeParserPtr> frame;

    // The broken parsers...
    QList<MoleculeParserPtr> broken;

    // The supplementary parsers
    QList<MoleculeParserPtr> supplementary;

    if (this->isBroken())
    {
        broken.append(*this);
    }
    else if (this->isSupplementary())
    {
        supplementary.append(*this);
    }
    else
    {
        if (this->isTopology())
        {
            topology.append(*this);
        }

        if (this->isFrame())
        {
            frame.append(*this);
        }
    }

    for (auto parser : parsers)
    {
        // This is a lead parser.
        if (parser.read().isBroken())
        {
            broken.append(parser);
        }
        else if (parser.read().isSupplementary())
        {
            supplementary.append(parser);
        }
        else
        {
            if (parser.read().isTopology())
            {
                topology.append(parser);
            }

            if (parser.read().isFrame())
            {
                frame.append(parser);
            }
        }
    }

    if (broken.count() > 0)
    {
        // everything is broken!
        ret["broken"] = broken;
        return ret;
    }

    // No topology parsers - we can't create a system...
    if (topology.count() == 0)
    {
        if (supplementary.count() > 0)
        {
            // likely a Supplementary is hiding a broken parser - reparse it...
            for (const auto &parser : supplementary)
            {
                PropertyMap m2(map);
                m2.set("DISABLE_SUPPLEMENTARY", BooleanProperty(true));

                auto p = _pvt_parse(parser.read().filename(), m2);

                if (p.read().isBroken())
                {
                    broken.append(p);
                }
            }

            if (broken.count() > 0)
            {
                ret["broken"] = broken;
                return ret;
            }
        }

        throw SireIO::parse_error(QObject::tr("Unable to load any molecules from the files as none "
                                              "contain the necessary molecular information to create a "
                                              "system. Only coordinate or trajectory information "
                                              "has been loaded. Structure or topology information, "
                                              "e.g. as would be found in a topology file, is missing."),
                                  CODELOC);
    }

    // If there are topology parsers. We want to use the first one that
    // can only be used as a topology parser (e.g. if the user has loaded
    // a prmtop and a PDB file in the wrong order)
    QList<int> top_or_frame;

    if (topology.count() > 1)
    {
        int topology_only_idx = -1;

        for (int i = 0; i < topology.count(); ++i)
        {
            if (not topology[i].read().isFrame())
            {
                if (topology_only_idx != -1)
                {
                    throw SireIO::parse_error(
                        QObject::tr("Cannot construct a System from multiple topology-only parsers "
                                    "if none can follow!"),
                        CODELOC);
                }

                topology_only_idx = i;
            }
            else
            {
                top_or_frame.append(i);
            }
        }

        if (topology_only_idx != -1)
        {
            topology = {topology[topology_only_idx]};
        }
        else
        {
            auto tmp = topology;
            tmp.clear();

            for (const int idx : top_or_frame)
            {
                tmp.append(topology[idx]);
            }

            topology = tmp;
        }
    }

    ret["topology"] = topology;
    ret["frame"] = frame;
    ret["supplementary"] = supplementary;

    return ret;
}

bool MoleculeParser::writingTrajectory() const
{
    return not frames_to_write.isEmpty();
}

QList<qint32> MoleculeParser::framesToWrite() const
{
    return frames_to_write;
}

SireMol::Frame MoleculeParser::createFrame(qint32 frame_index) const
{
    PropertyMap map(propmap);

    frame_index = Index(frame_index).map(saved_system.nFrames());

    if (map["frame_aligner"].hasValue())
    {
        map.set("transform",
                map["frame_aligner"].value().asA<TrajectoryAligner>()[frame_index]);
    }

    System local_system(saved_system);

    local_system.loadFrame(frame_index, map);

    return this->createFrame(local_system, map);
}

const PropertyMap &MoleculeParser::propertyMap() const
{
    return propmap;
}

QString MoleculeParser::saveTitle() const
{
    return saved_system.name().value();
}

Q_GLOBAL_STATIC(NullParser, nullParser)

const NullParser &MoleculeParser::null()
{
    return *(nullParser());
}

//////////////
////////////// Implementation of NullParser
//////////////

static const RegisterMetaType<NullParser> r_null;

QDataStream &operator<<(QDataStream &ds, const NullParser &parser)
{
    writeHeader(ds, r_null, 1);
    ds << static_cast<const MoleculeParser &>(parser);
    return ds;
}

QDataStream &operator>>(QDataStream &ds, NullParser &parser)
{
    VersionID v = readHeader(ds, r_null);

    if (v == 1)
    {
        ds >> static_cast<MoleculeParser &>(parser);
    }
    else
        throw version_error(v, "1", r_null, CODELOC);

    return ds;
}

NullParser::NullParser() : ConcreteProperty<NullParser, MoleculeParser>()
{
}

NullParser::NullParser(const NullParser &other) : ConcreteProperty<NullParser, MoleculeParser>(other)
{
}

NullParser::~NullParser()
{
}

NullParser &NullParser::operator=(const NullParser &other)
{
    MoleculeParser::operator=(other);
    return *this;
}

bool NullParser::operator==(const NullParser &other) const
{
    return MoleculeParser::operator==(other);
}

bool NullParser::operator!=(const NullParser &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *NullParser::typeName()
{
    return QMetaType::typeName(qMetaTypeId<NullParser>());
}

QString NullParser::formatName() const
{
    return "NULL";
}

/** Return a description of the file format */
QString NullParser::formatDescription() const
{
    return QObject::tr("Null parser that should not be used for any real parsing.");
}

int NullParser::nAtoms() const
{
    return 0;
}

System NullParser::toSystem(const PropertyMap &) const
{
    return System();
}

System NullParser::toSystem(const MoleculeParser &other, const PropertyMap &) const
{
    if (not other.isA<NullParser>())
        throw SireError::incompatible_error(
            QObject::tr("Null parsers cannot be combined with other parsers (%1)").arg(other.toString()), CODELOC);

    return System();
}

System NullParser::toSystem(const QList<MoleculeParserPtr> &others, const PropertyMap &) const
{
    for (auto other : others)
    {
        if (not other.isNull())
        {
            throw SireError::incompatible_error(
                QObject::tr("Null parsers cannot be combined with other parsers (%1)").arg(other.read().toString()),
                CODELOC);
        }
    }

    return System();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr NullParser::construct(const QString &filename, const PropertyMap &map) const
{
    throw SireError::program_bug(QObject::tr("The NullParser should not be used for an real file IO!"), CODELOC);

    return MoleculeParserPtr();
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr NullParser::construct(const QStringList &lines, const PropertyMap &map) const
{
    throw SireError::program_bug(QObject::tr("The NullParser should not be used for an real file IO!"), CODELOC);

    return MoleculeParserPtr();
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr NullParser::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    throw SireError::program_bug(QObject::tr("The NullParser should not be used for an real file IO!"), CODELOC);

    return MoleculeParserPtr();
}

//////////////
////////////// Implementation of BrokenParser
//////////////

static const RegisterMetaType<BrokenParser> r_broken;

QDataStream &operator<<(QDataStream &ds, const BrokenParser &parser)
{
    writeHeader(ds, r_broken, 1);

    SharedDataStream sds(ds);

    sds << parser.error_report << parser.suffix << static_cast<const MoleculeParser &>(parser);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BrokenParser &parser)
{
    VersionID v = readHeader(ds, r_broken);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> parser.error_report >> parser.suffix >> static_cast<MoleculeParser &>(parser);
    }
    else
        throw version_error(v, "1", r_broken, CODELOC);

    return ds;
}

BrokenParser::BrokenParser() : ConcreteProperty<BrokenParser, MoleculeParser>()
{
}

BrokenParser::BrokenParser(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<BrokenParser, MoleculeParser>(filename, map)
{
}

BrokenParser::BrokenParser(const QString &filename, const QString &s, const QStringList &errors)
    : ConcreteProperty<BrokenParser, MoleculeParser>(filename, PropertyMap())
{
    suffix = s;
    error_report = errors;
}

BrokenParser::BrokenParser(const QString &filename, const QStringList &errors)
    : ConcreteProperty<BrokenParser, MoleculeParser>(filename, PropertyMap())
{
    error_report = errors;
}

BrokenParser::BrokenParser(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<BrokenParser, MoleculeParser>(lines, map)
{
}

BrokenParser::BrokenParser(const System &, const PropertyMap &map) : ConcreteProperty<BrokenParser, MoleculeParser>(map)
{
}

BrokenParser::BrokenParser(const BrokenParser &other)
    : ConcreteProperty<BrokenParser, MoleculeParser>(other), error_report(other.error_report), suffix(other.suffix)
{
}

BrokenParser::~BrokenParser()
{
}

BrokenParser &BrokenParser::operator=(const BrokenParser &other)
{
    if (this != &other)
    {
        error_report = other.error_report;
        suffix = other.suffix;
        MoleculeParser::operator=(other);
    }

    return *this;
}

bool BrokenParser::operator==(const BrokenParser &other) const
{
    return MoleculeParser::operator==(other);
}

bool BrokenParser::operator!=(const BrokenParser &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *BrokenParser::typeName()
{
    return QMetaType::typeName(qMetaTypeId<BrokenParser>());
}

QString BrokenParser::formatName() const
{
    return "BROKEN";
}

/** Return a description of the file format */
QString BrokenParser::formatDescription() const
{
    return QObject::tr("Broken parser used to report an unparseable file.");
}

int BrokenParser::nAtoms() const
{
    return 0;
}

bool BrokenParser::isBroken() const
{
    return true;
}

QString BrokenParser::errorReport() const
{
    if (suffix.isEmpty())
    {
        return QObject::tr("== %1 ==\n\nThis file was not recognised by any of the file parsers!\n\n"
                           "%2\n")
            .arg(this->filename())
            .arg(error_report.join("\n\n"));
    }
    else
    {
        return QObject::tr("== %1 ==\n\nThis file could not be parsed by any of the file parsers! "
                           "It was recognised as a file of type %2, but all parsers failed "
                           "to parse this file. The errors from the parsers associated "
                           "with the suffix %2 are printed below:\n\n"
                           "%3\n")
            .arg(this->filename())
            .arg(suffix)
            .arg(error_report.join("\n\n"));
    }
}

System BrokenParser::toSystem(const PropertyMap &) const
{
    return System();
}

System BrokenParser::toSystem(const MoleculeParser &other, const PropertyMap &) const
{
    if (not other.isA<BrokenParser>())
        throw SireError::incompatible_error(
            QObject::tr("Broken parsers cannot be combined with other parsers (%1)").arg(other.toString()), CODELOC);

    return System();
}

System BrokenParser::toSystem(const QList<MoleculeParserPtr> &others, const PropertyMap &) const
{
    for (const auto &other : others)
    {
        if (not other->isA<BrokenParser>())
        {
            throw SireError::incompatible_error(
                QObject::tr("Broken parsers cannot be combined with other parsers (%1)").arg(other.read().toString()),
                CODELOC);
        }
    }

    return System();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr BrokenParser::construct(const QString &filename, const PropertyMap &map) const
{
    return MoleculeParserPtr(BrokenParser(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr BrokenParser::construct(const QStringList &lines, const PropertyMap &map) const
{
    return MoleculeParserPtr(BrokenParser(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr BrokenParser::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    return MoleculeParserPtr(BrokenParser(system, map));
}
