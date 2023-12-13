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

#include "forcefieldinfo.h"

#include "SireSystem/system.h"

#include "SireMol/molecule.h"
#include "SireMol/selectormol.h"
#include "SireMol/core.h"

#include "SireVol/cartesian.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"

using namespace SireSystem;
using namespace SireMol;
using namespace SireMM;
using namespace SireFF;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits::Dimension;
using namespace SireUnits;

static const RegisterMetaType<ForceFieldInfo> r_ffinfo;

ForceFieldInfo::CUTOFF_TYPE ForceFieldInfo::string_to_cutoff_type(QString s)
{
    s = s.toUpper().simplified().replace(" ", "_");

    if (s == "PARTICLE_MESH_EWALD")
        s = "PME";
    else if (s.isEmpty() or s == "NONE")
        s = "NO_CUTOFF";
    else if (s == "RF")
        s = "REACTION_FIELD";
    else if (s == "SHIFT")
        s = "SHIFT_ELECTROSTATICS";
    else if (s == "DEFAULT")
        s = "CUTOFF";

    static std::map<QString, ForceFieldInfo::CUTOFF_TYPE> types = {
        {"NO_CUTOFF", ForceFieldInfo::NO_CUTOFF},
        {"CUTOFF", ForceFieldInfo::CUTOFF},
        {"EWALD", ForceFieldInfo::EWALD},
        {"PME", ForceFieldInfo::PME},
        {"REACTION_FIELD", ForceFieldInfo::REACTION_FIELD},
        {"SHIFT_ELECTROSTATICS", ForceFieldInfo::SHIFT_ELECTROSTATICS}};

    auto it = types.find(s);

    if (it == types.end())
    {
        throw SireError::invalid_key(QObject::tr(
                                         "Unsupported cutoff type '%1'. Supported types are (case-insensitive): "
                                         "* PME / PARTICLE MESH EWALD / PARTICLE_MESH_EWALD * or "
                                         "* EWALD * or "
                                         "* REACTION FIELD / REACTION_FIELD / RF * or "
                                         "* SHIFT / SHIFT ELECTROSTATICS / SHIFT_ELECTROSTATICS * or "
                                         "* NONE / NO CUTOFF / NO_CUTOFF * or "
                                         "* DEFAULT / CUTOFF")
                                         .arg(s),
                                     CODELOC);

        return ForceFieldInfo::CUTOFF;
    }
    else
        return it->second;
}

QString ForceFieldInfo::cutoff_type_to_string(ForceFieldInfo::CUTOFF_TYPE typ)
{
    switch (typ)
    {
    case ForceFieldInfo::NO_CUTOFF:
        return "NO_CUTOFF";
    case ForceFieldInfo::CUTOFF:
        return "CUTOFF";
    case ForceFieldInfo::EWALD:
        return "EWALD";
    case ForceFieldInfo::PME:
        return "PME";
    case ForceFieldInfo::REACTION_FIELD:
        return "REACTION_FIELD";
    case ForceFieldInfo::SHIFT_ELECTROSTATICS:
        return "SHIFT_ELECTROSTATICS";
    default:
        return "CUTOFF";
    }
}

QDataStream &operator<<(QDataStream &ds, const ForceFieldInfo &ffinfo)
{
    writeHeader(ds, r_ffinfo, 1);

    SharedDataStream sds(ds);

    sds << ffinfo.spc << ffinfo.params << ffinfo.dtl
        << ffinfo.ctff.to(SireUnits::angstrom)
        << ForceFieldInfo::cutoff_type_to_string(ffinfo.ctff_typ);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, ForceFieldInfo &ffinfo)
{
    VersionID v = readHeader(ds, r_ffinfo);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        double l;
        QString cutoff_type;

        sds >> ffinfo.spc >> ffinfo.params >> ffinfo.dtl >> l >> cutoff_type;

        ffinfo.ctff = l * SireUnits::angstrom;
        ffinfo.ctff_typ = ForceFieldInfo::string_to_cutoff_type(cutoff_type);
    }
    else
        throw version_error(v, "1", r_ffinfo, CODELOC);

    return ds;
}

// this (specific) value is the largest we can have on Windows
// to still calculate an energy. Anything higher gives and energy of 0...
const auto MAX_CUTOFF = 347557 * angstrom;

ForceFieldInfo::ForceFieldInfo()
    : ConcreteProperty<ForceFieldInfo, Property>(),
      spc(Cartesian()), dtl(MMDetail()), ctff(MAX_CUTOFF), ctff_typ(NO_CUTOFF)
{
}

ForceFieldInfo::ForceFieldInfo(const System &system,
                               const SireBase::PropertyMap &map)
    : ConcreteProperty<ForceFieldInfo, Property>()
{
    this->operator=(ForceFieldInfo());

    const auto mols = system.molecules();

    const auto space_property = map["space"];

    if (space_property.hasValue())
    {
        this->setSpace(space_property.value().asA<Space>());
    }
    else if (system.containsProperty(space_property.source()))
    {
        this->setSpace(system.property(space_property.source()).asA<Space>());
    }

    const auto cutoff_prop = map["cutoff"];

    if (cutoff_prop.hasValue())
    {
        this->setCutoff(cutoff_prop.value().asA<GeneralUnitProperty>().toUnit<Length>());
    }
    else if (cutoff_prop.source() != "cutoff")
    {
        throw SireError::invalid_arg(QObject::tr(
                                         "The cutoff property should have a value. It cannot be the string "
                                         "'%1'. If you want to specify the cutoff type, using "
                                         "the 'cutoff_type' property.")
                                         .arg(cutoff_prop.source()),
                                     CODELOC);
    }

    const auto cutoff_type = map["cutoff_type"];

    if (cutoff_type.hasSource() and cutoff_type.source() != "cutoff_type")
    {
        this->setCutoffType(cutoff_type.source(), map);
    }

    const auto ff_prop = map["forcefield"];

    if (ff_prop.hasValue())
    {
        dtl = ff_prop.value().asA<FFDetail>();
    }
    else
    {
        // find the first FFDetail object
        auto it = mols.constBegin();

        for (; it != mols.constEnd(); ++it)
        {
            const auto &data = it.value().data();

            if (data.hasProperty(ff_prop.source()))
            {
                dtl = data.property(ff_prop.source()).asA<FFDetail>();
                ++it;
                break;
            }
        }

        const auto &ffdetail = dtl.read().asA<FFDetail>();

        // now check that any other forcefields are compatible
        for (; it != mols.constEnd(); ++it)
        {
            const auto &data = it.value().data();

            if (data.hasProperty(ff_prop.source()))
            {
                data.property(ff_prop.source()).asA<FFDetail>().assertCompatibleWith(ffdetail);
            }
        }
    }
}

ForceFieldInfo::ForceFieldInfo(const SireMol::MoleculeView &mol,
                               const SireBase::PropertyMap &map)
    : ConcreteProperty<ForceFieldInfo, Property>()
{
    System system;
    MoleculeGroup molgroup("all");
    molgroup.add(mol);
    system.add(molgroup);

    const auto space_prop = map["space"];

    if (not space_prop.hasValue())
    {
        const auto &moldata = mol.data();

        if (moldata.hasProperty(space_prop.source()))
        {
            system.setProperty("space", moldata.property(space_prop.source()));
        }
    }

    this->operator=(ForceFieldInfo(system, map));
}

ForceFieldInfo::ForceFieldInfo(const SireMol::SelectorMol &mols,
                               const SireBase::PropertyMap &map)
    : ConcreteProperty<ForceFieldInfo, Property>()
{
    System system;
    system.add(mols.toMoleculeGroup());

    const auto space_prop = map["space"];

    if (not space_prop.hasValue())
    {
        // use the first space found from the molecule
        for (int i = 0; i < mols.count(); ++i)
        {
            const auto &moldata = mols[i].data();

            if (moldata.hasProperty(space_prop.source()))
            {
                system.setProperty("space", moldata.property(space_prop.source()));
                break;
            }
        }
    }

    this->operator=(ForceFieldInfo(system, map));
}

ForceFieldInfo::ForceFieldInfo(const ForceFieldInfo &other)
    : ConcreteProperty<ForceFieldInfo, Property>(other),
      spc(other.spc), params(other.params), dtl(other.dtl),
      ctff(other.ctff), ctff_typ(other.ctff_typ)
{
}

ForceFieldInfo::~ForceFieldInfo()
{
}

ForceFieldInfo &ForceFieldInfo::operator=(const ForceFieldInfo &other)
{
    if (this != &other)
    {
        spc = other.spc;
        params = other.params;
        dtl = other.dtl;
        ctff = other.ctff;
        ctff_typ = other.ctff_typ;
    }

    return *this;
}

bool ForceFieldInfo::operator==(const ForceFieldInfo &other) const
{
    return spc == other.spc and params == other.params and dtl == other.dtl and
           ctff_typ == other.ctff_typ and ctff == other.ctff;
}

bool ForceFieldInfo::operator!=(const ForceFieldInfo &other) const
{
    return not this->operator==(other);
}

QString ForceFieldInfo::toString() const
{
    QStringList parts;

    parts.append(QObject::tr("space=%1").arg(spc.read().toString()));
    parts.append(QObject::tr("cutoff_type=%1").arg(cutoff_type_to_string(ctff_typ)));

    if (this->hasCutoff())
        parts.append(QObject::tr("cutoff=%1").arg(ctff.toString()));

    if (not params.isEmpty())
    {
        parts.append(QObject::tr("params=%1").arg(params.toString()));
    }

    parts.append(QObject::tr("detail=%1").arg(dtl.read().toString()));

    return QObject::tr("ForceFieldInfo(\n  %1\n)").arg(parts.join(",\n  "));
}

const char *ForceFieldInfo::what() const
{
    return ForceFieldInfo::typeName();
}

const char *ForceFieldInfo::typeName()
{
    return QMetaType::typeName(qMetaTypeId<ForceFieldInfo>());
}

ForceFieldInfo *ForceFieldInfo::clone() const
{
    return new ForceFieldInfo(*this);
}

const SireVol::Space &ForceFieldInfo::space() const
{
    return spc.read();
}

void ForceFieldInfo::setSpace(const SireVol::Space &space)
{
    spc = space;

    if (this->hasNoCutoff())
    {
        this->setLargestCutoff();
    }
}

SireUnits::Dimension::Length ForceFieldInfo::cutoff() const
{
    return ctff;
}

void ForceFieldInfo::setCutoff(SireUnits::Dimension::Length length)
{
    if (length.value() <= 0)
    {
        ctff = MAX_CUTOFF;
        ctff_typ = NO_CUTOFF;
    }
    else
    {
        if (this->space().isPeriodic())
        {
            if (length > this->space().maximumCutoff())
            {
                throw SireError::incompatible_error(QObject::tr(
                                                        "You cannot set a cutoff (%1) that is larger than "
                                                        "the maximum possible cutoff (%2) for the periodic "
                                                        "space %3.")
                                                        .arg(length.toString())
                                                        .arg(this->space().maximumCutoff().toString())
                                                        .arg(this->space().toString()),
                                                    CODELOC);
            }
        }

        ctff = length;

        if (ctff_typ == NO_CUTOFF)
            ctff_typ = CUTOFF;
    }
}

void ForceFieldInfo::setLargestCutoff()
{
    if (spc.read().isPeriodic())
    {
        this->setCutoff(spc.read().maximumCutoff() - 1 * angstrom);
    }
    else
    {
        this->setNoCutoff();
    }
}

bool ForceFieldInfo::hasCutoff() const
{
    return not this->hasNoCutoff();
}

bool ForceFieldInfo::hasNoCutoff() const
{
    return ctff.value() == 0 or ctff == MAX_CUTOFF;
}

void ForceFieldInfo::setNoCutoff()
{
    ctff = MAX_CUTOFF;
    ctff_typ = NO_CUTOFF;
}

QString ForceFieldInfo::cutoffType() const
{
    return cutoff_type_to_string(ctff_typ);
}

void ForceFieldInfo::setCutoffType(QString cutoff_type)
{
    this->setCutoffType(cutoff_type, PropertyMap());
}

void ForceFieldInfo::setCutoffType(QString s_cutoff_type,
                                   const PropertyMap &map)
{
    auto cutoff_type = string_to_cutoff_type(s_cutoff_type);

    if (cutoff_type == NO_CUTOFF)
    {
        if (this->space().isPeriodic())
        {
            throw SireError::incompatible_error(QObject::tr(
                                                    "You cannot disable the cutoff when using a periodic space: %1")
                                                    .arg(this->space().toString()),
                                                CODELOC);
        }

        params.clear();
        ctff = MAX_CUTOFF;
    }
    else
    {
        params.clear();

        if (ctff.value() == 0 or ctff == MAX_CUTOFF)
        {
            this->setLargestCutoff();
        }

        if (cutoff_type == EWALD or cutoff_type == PME)
        {
            // set the EWALD parameters
            const auto tolerance_prop = map["tolerance"];

            if (tolerance_prop.hasValue())
            {
                this->setParameter("tolerance", GeneralUnit(tolerance_prop.value().asADouble()));
            }
            else
            {
                this->setParameter("tolerance", GeneralUnit(0.0001));
            }
        }
        else if (cutoff_type == REACTION_FIELD)
        {
            const auto dielectric_prop = map["dielectric"];

            if (dielectric_prop.hasValue())
            {
                this->setParameter("dielectric", GeneralUnit(dielectric_prop.value().asADouble()));
            }
            else
            {
                this->setParameter("dielectric", GeneralUnit(78.3));
            }
        }
    }

    ctff_typ = cutoff_type;
}

QStringList ForceFieldInfo::cutoffTypes()
{
    QStringList types;

    types.append(cutoff_type_to_string(NO_CUTOFF));
    types.append(cutoff_type_to_string(CUTOFF));
    types.append(cutoff_type_to_string(EWALD));
    types.append(cutoff_type_to_string(PME));
    types.append(cutoff_type_to_string(REACTION_FIELD));
    types.append(cutoff_type_to_string(SHIFT_ELECTROSTATICS));

    return types;
}

GeneralUnit ForceFieldInfo::getParameter(const QString &parameter) const
{
    if (params.hasProperty(parameter))
        return params.property(parameter).asA<GeneralUnit>();
    else
        return GeneralUnit(0);
}

void ForceFieldInfo::setParameter(const QString &parameter,
                                  const GeneralUnit &value)
{
    params.setProperty(parameter, SireBase::wrap(value));
}

SireBase::Properties ForceFieldInfo::parameters() const
{
    return params;
}

void ForceFieldInfo::setDetail(const FFDetail &d)
{
    dtl = d;
}

const FFDetail &ForceFieldInfo::detail() const
{
    return dtl.read().asA<FFDetail>();
}
