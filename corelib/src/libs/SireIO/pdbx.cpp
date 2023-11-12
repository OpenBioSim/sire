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

#include "SireIO/pdbx.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/core.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/trajectory.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include <iostream>

using namespace SireBase;
using namespace SireIO;
using namespace SireMol;
using namespace SireSystem;
using namespace SireStream;

const RegisterParser<PDBx> register_pdbx;

static const RegisterMetaType<PDBx> r_pdbx;

QDataStream &operator<<(QDataStream &ds, const PDBx &pdbx)
{
    writeHeader(ds, r_pdbx, 1);

    SharedDataStream sds(ds);

    sds << pdbx.parsed_system << static_cast<const MoleculeParser &>(pdbx);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PDBx &pdbx)
{
    VersionID v = readHeader(ds, r_pdbx);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> pdbx.parsed_system >> static_cast<MoleculeParser &>(pdbx);
    }
    else
        throw version_error(v, "1", r_pdbx, CODELOC);

    return ds;
}

static PDBxReaderFunction reader_function;
static PDBxWriterFunction writer_function;

void SireIO::register_pdbx_loader_functions(const PDBxWriterFunction &writer,
                                            const PDBxReaderFunction &reader)
{
    writer_function = writer;
    reader_function = reader;
}

/** Constructor */
PDBx::PDBx() : ConcreteProperty<PDBx, MoleculeParser>()
{
}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
PDBx::PDBx(const QString &filename, const PropertyMap &map) : ConcreteProperty<PDBx, MoleculeParser>(filename, map)
{
    // the file has been read into memory and is available via
    // the MoleculeParser::lines() function.

    // a parameter has also been read in MoleculeParser to say whether
    // we are allowed to use multiple cores to parse the file, e.g.
    // MoleculeParser::usesParallel() will be true

    // parse the data in the parse function
    this->parseLines(map);

    // now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
PDBx::PDBx(const QStringList &lines, const PropertyMap &map) : ConcreteProperty<PDBx, MoleculeParser>(lines, map)
{
    // the file has been read into memory and is available via
    // the MoleculeParser::lines() function.

    // a parameter has also been read in MoleculeParser to say whether
    // we are allowed to use multiple cores to parse the file, e.g.
    // MoleculeParser::usesParallel() will be true

    // parse the data in the parse function
    this->parseLines(map);

    // now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
PDBx::PDBx(const SireSystem::System &system, const PropertyMap &map) : ConcreteProperty<PDBx, MoleculeParser>(system, map)
{
    if (writer_function.empty())
        throw SireError::unsupported(
            "No PDBx writer function has been registered. You need to "
            "install a library to write PDBx/mmCIF files, e.g. gemmi. "
            "Do this by running 'mamba install -c conda-forge gemmi' "
            "and then re-running this script.",
            CODELOC);

    MoleculeParser::setLines(writer_function(system, map).toVector());
    parsed_system = system;
}

/** Copy constructor */
PDBx::PDBx(const PDBx &other)
    : ConcreteProperty<PDBx, MoleculeParser>(other),
      parsed_system(other.parsed_system), parse_warnings(other.parse_warnings)
{
}

/** Destructor */
PDBx::~PDBx()
{
}

/** Copy assignment operator */
PDBx &PDBx::operator=(const PDBx &other)
{
    if (this != &other)
    {
        parsed_system = other.parsed_system;
        parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool PDBx::operator==(const PDBx &other) const
{
    return parsed_system == other.parsed_system and
           MoleculeParser::operator==(other);
}

/** Comparison operator */
bool PDBx::operator!=(const PDBx &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char *PDBx::typeName()
{
    return QMetaType::typeName(qMetaTypeId<PDBx>());
}

/** Return the C++ name for this class */
const char *PDBx::what() const
{
    return PDBx::typeName();
}

bool PDBx::isTopology() const
{
    return true;
}

bool PDBx::isFrame() const
{
    return true;
}

int PDBx::nFrames() const
{
    return 1;
}

Frame PDBx::getFrame(int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    throw SireError::unsupported(QObject::tr(
                                     "You cannot extra frame data from a PDBx/mmCIF file."),
                                 CODELOC);

    return Frame();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr PDBx::construct(const QString &filename, const PropertyMap &map) const
{
    return PDBx(filename, map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr PDBx::construct(const QStringList &lines, const PropertyMap &map) const
{
    return PDBx(lines, map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr PDBx::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    return PDBx(system, map);
}

/** Return a string representation of this parser */
QString PDBx::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("PDBx::null");
    else
    {
        return QObject::tr("PDBx( nAtoms() = %1 )")
            .arg(nAtoms());
    }
}

/** Convert the the parsed data to a collection of PDBx record lines. */
QVector<QString> PDBx::toLines() const
{
    if (writer_function.empty())
        throw SireError::unsupported(
            "No PDBx writer function has been registered. You need to "
            "install a library to write PDBx/mmCIF files, e.g. gemmi. "
            "Do this by running 'mamba install -c conda-forge gemmi' "
            "and then re-running this script.",
            CODELOC);

    return writer_function(this->parsed_system, this->propertyMap()).toVector();
}

/** Return the format name that is used to identify this file format within Sire */
QString PDBx::formatName() const
{
    return "PDBx";
}

/** Return a description of the file format */
QString PDBx::formatDescription() const
{
    return QObject::tr("Protein Data Bank PDBx/mmCIF format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList PDBx::formatSuffix() const
{
    static const QStringList suffixes = {"pdbx", "cif"};
    return suffixes;
}

/** Return the total number of atoms in all molecules. */
int PDBx::nAtoms() const
{
    return 0;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void PDBx::assertSane() const
{
    // check state, raise SireError::program_bug if we are in an invalid state
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void PDBx::parseLines(const PropertyMap &map)
{
    if (reader_function.empty())
        throw SireError::unsupported(
            "No PDBx reader function has been registered. You need to "
            "install a library to read PDBx/mmCIF files, e.g. gemmi. "
            "Do this by running 'mamba install -c conda-forge gemmi' "
            "and then re-running this script.",
            CODELOC);

    // only parse if the file contains a line with
    // '_audit_conform.dict_name       mmcif_pdbx.dic'

    bool is_conformant = false;

    for (const auto &line : lines())
    {
        if (line.contains("_audit_conform.dict_name") and line.contains("mmcif_pdbx.dic"))
        {
            is_conformant = true;
            break;
        }
    }

    if (not is_conformant)
        throw SireError::unsupported(
            QObject::tr("The file does not appear to be a PDBx/mmCIF file. "
                        "It does not contain the line '_audit_conform.dict_name mmcif_pdbx.dic'"),
            CODELOC);

    parsed_system = reader_function(lines().toList(), map);

    parsed_system.setProperty(map["fileformat"].source(), StringProperty(this->formatName()));

    this->setScore(parsed_system.nAtoms());
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map' */
System PDBx::startSystem(const PropertyMap &map) const
{
    return parsed_system;
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void PDBx::addToSystem(System &system, const PropertyMap &map) const
{
    // you should loop through each molecule in the system and work out
    // which ones are described in the file, and then add data from the file
    // to those molecules.
    throw SireError::unsupported(QObject::tr(
                                     "You cannot add data from a PDBx/mmCIF file to an existing system."),
                                 CODELOC);
}
