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
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits::Dimension;
using namespace SireUnits;

static const RegisterMetaType<ForceFieldInfo> r_ffinfo;

QDataStream &operator<<(QDataStream &ds, const ForceFieldInfo &ffinfo)
{
    writeHeader(ds, r_ffinfo, 1);

    SharedDataStream sds(ds);

    sds << ffinfo.spc << ffinfo.params
        << ffinfo.ctff_typ << ffinfo.ctff.to(SireUnits::angstrom);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, ForceFieldInfo &ffinfo)
{
    VersionID v = readHeader(ds, r_ffinfo);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        double l;

        sds >> ffinfo.spc >> ffinfo.params >> ffinfo.ctff_typ >> l;

        ffinfo.ctff = l * SireUnits::angstrom;
    }
    else
        throw version_error(v, "1", r_ffinfo, CODELOC);

    return ds;
}

const QString NO_CUTOFF = "NO_CUTOFF";

ForceFieldInfo::ForceFieldInfo()
    : ConcreteProperty<ForceFieldInfo, Property>(),
      spc(Cartesian()), ctff_typ(NO_CUTOFF), ctff(0)
{
}

ForceFieldInfo::ForceFieldInfo(const System &system,
                               const SireBase::PropertyMap &map)
    : ConcreteProperty<ForceFieldInfo, Property>()
{
    this->operator=(ForceFieldInfo());

    const auto space_property = map["space"];

    if (space_property.hasValue())
    {
        this->setSpace(space_property.value().asA<Space>());
    }
    else if (system.containsProperty(space_property.source()))
    {
        this->setSpace(system.property(space_property.source()).asA<Space>());
    }

    // construct from the system
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
      spc(other.spc), params(other.params), ctff_typ(other.ctff_typ),
      ctff(other.ctff)
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
        ctff_typ = other.ctff_typ;
        ctff = other.ctff;
    }

    return *this;
}

bool ForceFieldInfo::operator==(const ForceFieldInfo &other) const
{
    return spc == other.spc and params == other.params and
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
    parts.append(QObject::tr("cutoff_type=%1").arg(ctff_typ));

    if (this->hasCutoff())
        parts.append(QObject::tr("cutoff=%1").arg(ctff.toString()));

    if (not params.isEmpty())
    {
        parts.append(QObject::tr("params=%1").arg(params.toString()));
    }

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
        ctff = Length(0);
        ctff_typ = NO_CUTOFF;
    }
    else if (ctff_typ == NO_CUTOFF)
    {
        ctff_typ = "CUTOFF";
        ctff = length;
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
    return ctff.value() > 0;
}

bool ForceFieldInfo::hasNoCutoff() const
{
    return ctff.value() == 0;
}

void ForceFieldInfo::setNoCutoff()
{
    ctff = Length(0);
    ctff_typ = NO_CUTOFF;
}

QString ForceFieldInfo::cutoffType() const
{
    return ctff_typ;
}

void ForceFieldInfo::setCutoffType(const QString &cutoff_type)
{
    // will likely want to validate this...
    ctff_typ = cutoff_type;
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
