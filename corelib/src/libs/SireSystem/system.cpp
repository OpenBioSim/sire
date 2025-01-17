/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <QHash>
#include <QMutex>

#include "delta.h"
#include "monitorname.h"
#include "system.h"
#include "systemtrajectory.h"

#include "SireFF/energytable.h"
#include "SireFF/ff.h"
#include "SireFF/ffidx.h"
#include "SireFF/ffmolgroup.h"
#include "SireFF/fieldtable.h"
#include "SireFF/forcefield.h"
#include "SireFF/forcetable.h"
#include "SireFF/potentialtable.h"
#include "SireFF/probe.h"

#include "SireMM/ljparameterdb.h"

#include "SireMol/atomcoords.h"
#include "SireMol/molecule.h"
#include "SireMol/moleculegroup.h"
#include "SireMol/molecules.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/selectormol.h"

#include "SireBase/savestate.h"
#include "SireBase/timeproperty.h"
#include "SireBase/lazyevaluator.h"

#include "SireError/errors.h"
#include "SireMol/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>
#include <QTime>

using namespace SireSystem;
using namespace SireFF;
using namespace SireMol;
using namespace SireCAS;
using namespace SireVol;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

static VersionRegistry<QUuid> *system_registry(0);

static VersionRegistry<QUuid> &systemRegistry()
{
    if (system_registry == 0)
    {
        // need to make this thread-safe...
        system_registry = new VersionRegistry<QUuid>();
    }

    return *system_registry;
}

////////
//////// Implementation of SystemData
////////

/** Rebuild the index of molecule groups and molecules in the system */
void System::reindex()
{
    MolGroupsBase::clearIndex();

    QList<MGNum> mgnums0 = molgroups[0]->mgNums();
    QList<MGNum> mgnums1 = molgroups[1]->mgNums();

    foreach (MGNum mgnum, mgnums0)
    {
        MolGroupsBase::addToIndex(molgroups[0]->at(mgnum));
    }

    foreach (MGNum mgnum, mgnums1)
    {
        MolGroupsBase::addToIndex(molgroups[1]->at(mgnum));
    }
}

/** Rebuild the index of molecule groups in this System */
void System::rebuildIndex()
{
    // rebuild our index...
    mgroups_by_num.clear();

    QList<MGNum> mgnums0 = molgroups[0]->mgNums();
    QList<MGNum> mgnums1 = molgroups[1]->mgNums();

    mgroups_by_num.reserve(mgnums0.count() + mgnums1.count());

    foreach (MGNum mgnum, mgnums0)
    {
        mgroups_by_num.insert(mgnum, 0);
    }

    foreach (MGNum mgnum, mgnums1)
    {
        if (mgroups_by_num.contains(mgnum))
            throw SireError::program_bug(QObject::tr("It should not be possible for a System to contain two molecule "
                                                     "groups that both contain the same number... (%1)")
                                             .arg(mgnum),
                                         CODELOC);

        mgroups_by_num.insert(mgnum, 1);
    }

    // now rebuild the index of molecules
    MolGroupsBase::clearIndex();

    foreach (MGNum mgnum, mgnums0)
    {
        MolGroupsBase::addToIndex(molgroups[0]->at(mgnum));
    }

    foreach (MGNum mgnum, mgnums1)
    {
        MolGroupsBase::addToIndex(molgroups[1]->at(mgnum));
    }
}

static const RegisterMetaType<System> r_system;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const System &system)
{
    if (system.needsAccepting())
        qDebug() << "SYSTEM NEEDS ACCEPTING";

    writeHeader(ds, r_system, 4);

    if (system.subversion != 0)
        throw SireError::program_bug(QObject::tr("It is a mistake to try and save a system that is in a "
                                                 "temporarily invalid state (i.e. has non-zero subversion number). "
                                                 "The subversion number for %1 is %2.")
                                         .arg(system.toString())
                                         .arg(system.subversion),
                                     CODELOC);

    SharedDataStream sds(ds);

    // first try to save all of the loaded LJ parameter types. This will
    // help ensure that the LJID of parameters don't change too much between
    // save and loads, which will help with memory consumption...
    SireMM::LJDBIOLock dblock = SireMM::LJParameterDB::saveParameters(sds);

    sds << system.uid << system.sysname << system.molgroups[0] << system.molgroups[1] << system.sysmonitors
        << system.cons << system.shared_properties << static_cast<const MolGroupsBase &>(system);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, System &system)
{
    VersionID v = readHeader(ds, r_system);

    if (v == 4)
    {
        SharedDataStream sds(ds);

        SireMM::LJDBIOLock dblock = SireMM::LJParameterDB::loadParameters(sds);

        sds >> system.uid >> system.sysname >> system.molgroups[0] >> system.molgroups[1] >> system.sysmonitors >>
            system.cons >> system.shared_properties >> static_cast<MolGroupsBase &>(system);

        system.rebuildIndex();

        system.sysversion = systemRegistry().registerObject(system.uid);
        system.subversion = 0;
        system.shared_properties.clear();
    }
    else if (v == 3)
    {
        SharedDataStream sds(ds);

        SireMM::LJDBIOLock dblock = SireMM::LJParameterDB::loadParameters(sds);

        sds >> system.uid >> system.sysname >> system.molgroups[0] >> system.molgroups[1] >> system.sysmonitors >>
            system.cons >> static_cast<MolGroupsBase &>(system);

        system.rebuildIndex();

        system.sysversion = systemRegistry().registerObject(system.uid);
        system.subversion = 0;
        system.shared_properties.clear();
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> system.uid >> system.sysname >> system.molgroups[0] >> system.molgroups[1] >> system.sysmonitors >>
            system.cons >> static_cast<MolGroupsBase &>(system);

        system.rebuildIndex();

        system.sysversion = systemRegistry().registerObject(system.uid);
        system.subversion = 0;
        system.shared_properties.clear();
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> system.uid >> system.sysname >> system.molgroups[0] >> system.molgroups[1] >> system.sysmonitors >>
            static_cast<MolGroupsBase &>(system);

        system.rebuildIndex();

        system.sysversion = systemRegistry().registerObject(system.uid);
        system.subversion = 0;
        system.shared_properties.clear();
    }
    else
        throw version_error(v, "1-4", r_system, CODELOC);

    return ds;
}

/** Construct an unnamed System */
System::System()
    : ConcreteProperty<System, MolGroupsBase>(), uid(QUuid::createUuid()),
      sysversion(systemRegistry().registerObject(uid)),
      subversion(0)
{
    molgroups[0] = ForceFields();
    molgroups[1] = MoleculeGroups();
}

const System &System::null()
{
    return *(create_shared_null<System>());
}

/** Construct a named System */
System::System(const QString &name)
    : ConcreteProperty<System, MolGroupsBase>(), uid(QUuid::createUuid()), sysname(name),
      sysversion(systemRegistry().registerObject(uid)),
      subversion(0)
{
    molgroups[0] = ForceFields();
    molgroups[1] = MoleculeGroups();
}

/** Copy constructor */
System::System(const System &other)
    : ConcreteProperty<System, MolGroupsBase>(other), uid(other.uid), sysname(other.sysname),
      sysversion(other.sysversion), sysmonitors(other.sysmonitors), cons(other.cons),
      mgroups_by_num(other.mgroups_by_num),
      shared_properties(other.shared_properties),
      system_trajectory(other.system_trajectory),
      subversion(other.subversion)
{
    molgroups[0] = other.molgroups[0];
    molgroups[1] = other.molgroups[1];
}

/** Destructor */
System::~System()
{
}

/** Copy assignment operator */
System &System::operator=(const System &other)
{
    if (this != &other)
    {
        uid = other.uid;
        sysname = other.sysname;
        sysversion = other.sysversion;
        molgroups[0] = other.molgroups[0];
        molgroups[1] = other.molgroups[1];
        sysmonitors = other.sysmonitors;
        cons = other.cons;
        mgroups_by_num = other.mgroups_by_num;
        shared_properties = other.shared_properties;
        system_trajectory = other.system_trajectory;
        subversion = other.subversion;

        MolGroupsBase::operator=(other);
    }

    return *this;
}

/** Comparison operator - two systems are equal if they have the
    same UID and version */
bool System::operator==(const System &other) const
{
    return uid == other.uid and sysversion == other.sysversion and subversion == other.subversion;
}

/** Comparison operator - two systems are equal if they have the
    same UID and version */
bool System::operator!=(const System &other) const
{
    return not this->operator==(other);
}

/** Return a modifiable reference to all of the forcefields in
    this system */
ForceFields &System::_pvt_forceFields()
{
    BOOST_ASSERT(molgroups[0]->isA<ForceFields>());
    return molgroups[0].edit().asA<ForceFields>();
}

/** Return a reference to all of the forcefields in this system */
const ForceFields &System::_pvt_forceFields() const
{
    BOOST_ASSERT(molgroups[0]->isA<ForceFields>());
    return molgroups[0]->asA<ForceFields>();
}

/** Return a reference to all of the forcefields in this system */
const ForceFields &System::_pvt_constForceFields() const
{
    return this->_pvt_forceFields();
}

/** Return a modifiable reference to all of the non-forcefield
    molecule groups in this system */
MoleculeGroups &System::_pvt_moleculeGroups()
{
    BOOST_ASSERT(molgroups[1]->isA<MoleculeGroups>());
    return molgroups[1].edit().asA<MoleculeGroups>();
}

/** Return a reference to all of the non-forcefield
    molecule groups in this system */
const MoleculeGroups &System::_pvt_moleculeGroups() const
{
    BOOST_ASSERT(molgroups[1]->isA<MoleculeGroups>());
    return molgroups[1]->asA<MoleculeGroups>();
}

/** Return a reference to all of the non-forcefield
    molecule groups in this system */
const MoleculeGroups &System::_pvt_constMoleculeGroups() const
{
    return this->_pvt_moleculeGroups();
}

void System::_pvt_throwMissingGroup(MGNum mgnum) const
{
    throw SireMol::missing_group(QObject::tr("There is no molecule group with number %1 in this system. "
                                             "Available groups have numbers %2.")
                                     .arg(mgnum)
                                     .arg(Sire::toString(this->mgNums())),
                                 CODELOC);
}

/** Return the set of molecule groups that contains the molecule group
    with number 'mgnum'

    \throw SireMol::missing_group
*/
MolGroupsBase &System::_pvt_moleculeGroups(MGNum mgnum)
{
    int idx = mgroups_by_num.value(mgnum, -1);

    if (idx == -1)
        this->_pvt_throwMissingGroup(mgnum);

    return molgroups[idx].edit();
}

/** Return the set of molecule groups that contains the molecule group
    with number 'mgnum'

    \throw SireMol::missing_group
*/
const MolGroupsBase &System::_pvt_moleculeGroups(MGNum mgnum) const
{
    int idx = mgroups_by_num.value(mgnum, -1);

    if (idx == -1)
        this->_pvt_throwMissingGroup(mgnum);

    return molgroups[idx].read();
}

/** Return the set of molecule groups that contains the molecule group
    with number 'mgnum'

    \throw SireMol::missing_group
*/
const MolGroupsBase &System::_pvt_constMoleculeGroups(MGNum mgnum) const
{
    return this->_pvt_moleculeGroups(mgnum);
}

/** Return the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
MoleculeGroup &System::_pvt_moleculeGroup(MGNum mgnum)
{
    int idx = mgroups_by_num.value(mgnum, -1);

    if (idx == -1)
        this->_pvt_throwMissingGroup(mgnum);

    return SireMol::detail::get_editable_group(molgroups[idx].edit(), mgnum);
}

/** Return the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup &System::_pvt_moleculeGroup(MGNum mgnum) const
{
    int idx = mgroups_by_num.value(mgnum, -1);

    if (idx == -1)
        this->_pvt_throwMissingGroup(mgnum);

    return molgroups[idx]->at(mgnum);
}

/** Internal function used to get the group with number 'mgnum' */
MoleculeGroup &System::getGroup(MGNum mgnum)
{
    return this->_pvt_moleculeGroup(mgnum);
}

/** Internal function used to get the group with number 'mgnum' */
const MoleculeGroup &System::getGroup(MGNum mgnum) const
{
    return this->_pvt_moleculeGroup(mgnum);
}

/** Internal function used to get lots of groups */
void System::getGroups(const QList<MGNum> &mgnums, QVarLengthArray<const MoleculeGroup *, 10> &groups) const
{
    groups.clear();

    foreach (MGNum mgnum, mgnums)
    {
        groups.append(&(this->_pvt_moleculeGroup(mgnum)));
    }
}

/** Internal function used to get pointers to all of the groups in this system */
QHash<MGNum, const MoleculeGroup *> System::getGroups() const
{
    QHash<MGNum, const MoleculeGroup *> groups;

    QList<MGNum> mgnums0 = molgroups[0]->mgNums();
    QList<MGNum> mgnums1 = molgroups[1]->mgNums();

    groups.reserve(mgnums0.count() + mgnums1.count());

    foreach (MGNum mgnum, mgnums0)
    {
        groups.insert(mgnum, &(molgroups[0]->at(mgnum)));
    }

    foreach (MGNum mgnum, mgnums1)
    {
        groups.insert(mgnum, &(molgroups[1]->at(mgnum)));
    }

    return groups;
}

/** Return the forcefield that matches the ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FF &System::operator[](const FFID &ffid) const
{
    return this->_pvt_forceFields()[ffid];
}

/** Return the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
const SystemMonitor &System::operator[](const MonitorID &monid) const
{
    return sysmonitors[monid];
}

/** Return the molecule group at ID 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
const MoleculeGroup &System::operator[](const MGID &mgid) const
{
    return MolGroupsBase::operator[](mgid);
}

/** Return the molecule with number 'molnum'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
*/
ViewsOfMol System::operator[](MolNum molnum) const
{
    return MolGroupsBase::operator[](molnum);
}

/** Return the molecule with ID 'molid'

    \throw SireMol::missing_molecule
    \throw SireMol::duplicate_molecule
    \throw SireError::invalid_index
*/
ViewsOfMol System::operator[](const MolID &molid) const
{
    return MolGroupsBase::operator[](molid);
}

ViewsOfMol System::operator[](int i) const
{
    return MolGroupsBase::operator[](i);
}

ViewsOfMol System::operator[](const QString &name) const
{
    return MolGroupsBase::operator[](name);
}

QList<MolViewPtr> System::operator[](const SireBase::Slice &slice) const
{
    return MolGroupsBase::operator[](slice);
}

/** Convienient syntax for System::add */
System &System::operator+=(const FF &forcefield)
{
    this->add(forcefield);
    return *this;
}

/** Convienient syntax for System::add */
System &System::operator+=(const MoleculeGroup &molgroup)
{
    this->add(molgroup);
    return *this;
}

/** Convienient syntax for System::add */
System &System::operator+=(const Constraint &constraint)
{
    this->add(constraint);
    return *this;
}

/** Convienient syntax for System::add */
System &System::operator+=(const Constraints &constraints)
{
    this->add(constraints);
    return *this;
}

/** Convienient syntax for System::remove */
System &System::operator-=(const FF &forcefield)
{
    this->remove(forcefield);
    return *this;
}

/** Convienient syntax for System::remove */
System &System::operator-=(const MoleculeGroup &molgroup)
{
    this->remove(molgroup);
    return *this;
}

/** Convienient syntax for System::remove */
System &System::operator-=(const FFID &ffid)
{
    this->remove(ffid);
    return *this;
}

/** Convienient syntax for System::remove */
System &System::operator-=(const MGID &mgid)
{
    this->remove(mgid);
    return *this;
}

/** Convienient syntax for System::remove */
System &System::operator-=(const MolID &molid)
{
    this->remove(molid);
    return *this;
}

/** Convienient syntax for System::remove */
System &System::operator-=(const Constraint &constraint)
{
    this->remove(constraint);
    return *this;
}

/** Convienient syntax for System::remove */
System &System::operator-=(const Constraints &constraints)
{
    this->remove(constraints);
    return *this;
}

/** Set the name of this system */
void System::setName(const SysName &newname)
{
    if (sysname != newname)
    {
        sysname = newname;
        sysversion.incrementMajor();
    }
}

/** Set the name of this system */
void System::setName(const QString &newname)
{
    this->setName(SysName(newname));
}

/** Collect statistics about the current configuration */
void System::collectStats()
{
    sysmonitors.monitor(*this);
    // sysversion.incrementMajor();
    // cons.committed(*this);
}

/** Return the forcefield with ID 'ffid' in this system

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FF &System::at(const FFID &ffid) const
{
    return this->_pvt_forceFields().at(ffid);
}

/** Return the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
const SystemMonitor &System::at(const MonitorID &monid) const
{
    return sysmonitors.at(monid);
}

/** Return the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
const SystemMonitor &System::monitor(const MonitorID &monid) const
{
    return sysmonitors.monitor(monid);
}

/** Return the monitors with ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireError::invalid_index
*/
QList<SysMonPtr> System::monitors(const MonitorID &monid) const
{
    return sysmonitors.monitors(monid);
}

/** Return the forcefield with ID 'ffid' in this system

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FF &System::forceField(const FFID &ffid) const
{
    return this->_pvt_forceFields().forceField(ffid);
}

/** Return the forcefield that contains the molecule group
    identified by the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireMol::duplicate_group
    \throw SireError::invalid_index
*/
const FF &System::forceField(const MGID &mgid) const
{
    return this->_pvt_forceFields().forceField(this->getGroupNumber(mgid));
}

/** Return the number of monitors in this system */
int System::nMonitors() const
{
    return sysmonitors.count();
}

/** Return the number of forcefields in this system */
int System::nForceFields() const
{
    return this->_pvt_forceFields().nForceFields();
}

/** Return the number of constraints on the system */
int System::nConstraints() const
{
    return cons.count();
}

/** Return the index of the forcefield with ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
FFIdx System::ffIdx(const FFID &ffid) const
{
    return this->_pvt_forceFields().ffIdx(ffid);
}

/** Return the name of the monitor at ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireSystem::duplicate_monitor
    \throw SireError::invalid_index
*/
MonitorName System::monitorName(const MonitorID &monid) const
{
    return sysmonitors.monitorName(monid);
}

/** Return the name of the forcefield with ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireFF::duplicate_forcefield
    \throw SireError::invalid_index
*/
const FFName &System::ffName(const FFID &ffid) const
{
    return this->_pvt_forceFields().ffName(ffid);
}

/** Return a string representation of this system */
QString System::toString() const
{
    QStringList parts;

    parts.append(QString("name=%1").arg(this->name()));
    parts.append(QString("num_molecules=%1").arg(this->nMolecules()));

    if (this->nForceFields() > 0)
        parts.append(QString("num_force_fields=%1").arg(this->nForceFields()));

    if (this->nMonitors() > 0)
        parts.append(QString("num_monitors=%1").arg(this->nMonitors()));

    if (parts.count() < 4)
    {
        parts.insert(2, QString("num_atoms=%1").arg(this->nAtoms()));
    }

    if (parts.count() < 5)
    {
        parts.insert(2, QString("num_residues=%1").arg(this->nResidues()));
    }

    return QString("System( %1 )").arg(parts.join(" "));
}

/** Return the symbol that represents the total energy component
    of the system */
const Symbol &System::totalComponent() const
{
    return this->_pvt_forceFields().totalComponent();
}

/** Return the total energy of this system. */
MolarEnergy System::energy()
{
    return this->_pvt_forceFields().energy();
}

/** Return the total energy of the energy component in this system
    that is identified by the energy component 'component'

    \throw SireFF::missing_component
*/
MolarEnergy System::energy(const Symbol &component)
{
    return this->_pvt_forceFields().energy(component);
}

/** Return the total energytable in this system
 */
void System::energy(EnergyTable &energytable, double scale_energy)
{
    this->_pvt_forceFields().energy(energytable, scale_energy);
}
/** Return the total energytable of the energy component in this system
    that is identified by the energy component 'component'
 */
void System::energy(EnergyTable &energytable, const Symbol &component, double scale_energy)
{
    this->_pvt_forceFields().energy(energytable, component, scale_energy);
}

/** Return the energies of the energy components of this system whose
    symbols are in 'components'

    \throw SireFF::missing_component
*/
Values System::energies(const QSet<Symbol> &components)
{
    return this->_pvt_forceFields().energies(components);
}

/** Return the energies of all energy components in this system */
Values System::energies()
{
    return this->_pvt_forceFields().energies();
}

/** Return whether or not the component 'component' is an energy component

    \throw SireFF::missing_component
*/
bool System::isEnergyComponent(const Symbol &component) const
{
    return this->_pvt_forceFields().isEnergyComponent(component);
}

/** Return whether or not this system has an energy component 'component' */
bool System::hasEnergyComponent(const Symbol &component) const
{
    return this->_pvt_forceFields().hasEnergyComponent(component);
}

/** Set the energy component 'symbol' equal to the expression 'expression' */
void System::setEnergyComponent(const Symbol &symbol, const Expression &expression)
{
    if (this->hasComponent(symbol))
    {
        if (this->componentExpression(symbol) == expression)
            return;
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        this->_pvt_forceFields().setEnergyComponent(symbol, expression);
        sysversion.incrementMajor();
        this->applyConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Return the symbols that represent the energy expressions of this system */
QSet<Symbol> System::energySymbols() const
{
    return this->_pvt_forceFields().energySymbols();
}

/** Return all of the energy components of this system */
Values System::energyComponents()
{
    return this->_pvt_forceFields().energyComponents();
}

/** Return the energy expression for the energy component 'component'

    \throw SireFF::missing_component
*/
Expression System::energyExpression(const Symbol &component) const
{
    return this->_pvt_forceFields().energyExpression(component);
}

/** Return the energy expressions for the energy components whose
    symbols are in 'symbols'

    \throw SireFF::missing_component
*/
QHash<Symbol, Expression> System::energyExpressions(const QSet<Symbol> &symbols) const
{
    return this->_pvt_forceFields().energyExpressions(symbols);
}

/** Return all of the energy expressions in this system */
QHash<Symbol, Expression> System::energyExpressions() const
{
    return this->_pvt_forceFields().energyExpressions();
}

/** Return the constant value for the constant component 'component'

    \throw SireFF::missing_component
*/
double System::constant(const Symbol &component) const
{
    return this->_pvt_forceFields().constant(component);
}

/** Return the values of all constant components in this system */
Values System::constants() const
{
    return this->_pvt_forceFields().constants();
}

/** Return the values of the constant components whose symbols
    are in 'components'

    \throw SireFF::missing_component
*/
Values System::constants(const QSet<Symbol> &components) const
{
    return this->_pvt_forceFields().constants(components);
}

/** Return whether or not the system component 'component'
    is a constant component

    \throw SireFF::missing_component
*/
bool System::isConstantComponent(const Symbol &component) const
{
    return this->_pvt_forceFields().isConstantComponent(component);
}

/** Return whether or not this system has a constant
    component with symbol 'component' */
bool System::hasConstantComponent(const Symbol &component) const
{
    return this->_pvt_forceFields().hasConstantComponent(component);
}

/** Set the constant component 'symbol' to the value 'value'

    \throw SireError::incompatible_error
*/
void System::setConstantComponent(const Symbol &symbol, double value)
{
    Delta delta(*this);
    delta.update(symbol, value);
    this->operator=(delta.apply());
}

/** Set the constant component 'symbol' to the 'expression'

    \throw SireError::incompatible_error
*/
void System::setConstantComponent(const Symbol &symbol, const Expression &expression)
{
    if (this->hasConstantComponent(symbol))
    {
        Expression ex = this->constantExpression(symbol);

        if (ex == expression)
            return;
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        this->_pvt_forceFields().setConstantComponent(symbol, expression);
        sysversion.incrementMajor();
        this->applyConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

void System::setConstant(const Symbol &symbol, double value)
{
    this->setConstantComponent(symbol, value);
}

void System::setConstant(const Symbol &symbol, const Expression &expression)
{
    this->setConstantComponent(symbol, expression);
}

/** Return the symbols that represent constant components of this system */
QSet<Symbol> System::constantSymbols() const
{
    return this->_pvt_forceFields().constantSymbols();
}

/** Return the values of all constant components of this system */
Values System::constantComponents() const
{
    return this->_pvt_forceFields().constantComponents();
}

/** Return the expression that defines the constant component with
    symbol 'symbol'

    \throw SireFF::missing_component
*/
Expression System::constantExpression(const Symbol &symbol) const
{
    return this->_pvt_forceFields().constantExpression(symbol);
}

/** Return the expressions that define the constant components
    whose symbols are in 'symbols'

    \throw SireFF::missing_component
*/
QHash<Symbol, Expression> System::constantExpressions(const QSet<Symbol> &symbols) const
{
    return this->_pvt_forceFields().constantExpressions(symbols);
}

/** Return all of the expressions that define the constant components
    of this system */
QHash<Symbol, Expression> System::constantExpressions() const
{
    return this->_pvt_forceFields().constantExpressions();
}

/** Synonym for System::setConstantComponent(symbol, value) */
void System::setComponent(const Symbol &symbol, double value)
{
    this->setConstantComponent(symbol, value);
}

/** Synonym for System::setEnergyComponent(symbol, expression) */
void System::setComponent(const Symbol &symbol, const SireCAS::Expression &expression)
{
    if (expression.isConstant())
        this->setConstantComponent(symbol, expression);
    else
        this->setEnergyComponent(symbol, expression);
}

/** Return all of the symbols that represent the constant and
    energy components of this system */
QSet<Symbol> System::componentSymbols() const
{
    return this->_pvt_forceFields().componentSymbols();
}

/** Return the values of all components of this system
    (constant components and energies) */
Values System::componentValues()
{
    return this->_pvt_forceFields().componentValues();
}

/** Return the value of the energy or constant component
    with symbol 'symbol'

    \throw SireFF::missing_component
*/
double System::componentValue(const Symbol &symbol)
{
    return this->_pvt_forceFields().componentValue(symbol);
}

/** Retunr the value of the energy or constant component values
    whose symbols are in 'symbols'

    \throw SireFF::missing_component
*/
Values System::componentValues(const QSet<Symbol> &symbols)
{
    return this->_pvt_forceFields().componentValues(symbols);
}

/** Return whether or not this system has a constant or energy
    component represented by the symbol 'symbol' */
bool System::hasComponent(const Symbol &symbol) const
{
    return this->_pvt_forceFields().hasComponent(symbol);
}

/** Return the expression that defines the component represented
    by the symbol 'symbol'

    \throw SireFF::missing_component
*/
Expression System::componentExpression(const Symbol &symbol) const
{
    return this->_pvt_forceFields().componentExpression(symbol);
}

/** Return the expressions that define the components whose
    symbols are in 'symbols'

    \throw SireFF::missing_component
*/
QHash<Symbol, Expression> System::componentExpressions(const QSet<Symbol> &symbols) const
{
    return this->_pvt_forceFields().componentExpressions(symbols);
}

/** Return all of the expressions that define all of the components
    of this system */
QHash<Symbol, Expression> System::componentExpressions() const
{
    return this->_pvt_forceFields().componentExpressions();
}

/** Add the forces acting on the molecules in the forcetable 'forcetable'
    from this system onto this forcetable, scaled by the optionally
    supplied 'scale_force' */
void System::force(ForceTable &forcetable, double scale_force)
{
    this->_pvt_forceFields().force(forcetable, scale_force);
}

/** Add the forces acting on the molecules in the forcetable 'forcetable'
    from the component of this system identified by 'component' onto
    this forcetable, scaled by the optionally supplied 'scale_force' */
void System::force(ForceTable &forcetable, const Symbol &component, double scale_force)
{
    this->_pvt_forceFields().force(forcetable, component, scale_force);
}

/** Add the fields acting on the molecules in the fieldtable 'fieldtable'
    from this system onto this fieldtable, scaled by the optionally
    supplied 'scale_field' */
void System::field(FieldTable &fieldtable, double scale_field)
{
    this->_pvt_forceFields().field(fieldtable, scale_field);
}

/** Add the fields acting on the molecules in the fieldtable 'fieldtable'
    from the component of this system identified by 'component' onto
    this fieldtable, scaled by the optionally supplied 'scale_field' */
void System::field(FieldTable &fieldtable, const Symbol &component, double scale_field)
{
    this->_pvt_forceFields().field(fieldtable, component, scale_field);
}

/** Add the fields acting on the molecules in the fieldtable 'fieldtable'
    from this system onto this fieldtable, scaled by the optionally
    supplied 'scale_field' */
void System::field(FieldTable &fieldtable, const Probe &probe, double scale_field)
{
    this->_pvt_forceFields().field(fieldtable, probe, scale_field);
}

/** Add the fields acting on the molecules in the fieldtable 'fieldtable'
    from the component of this system identified by 'component' onto
    this fieldtable, scaled by the optionally supplied 'scale_field' */
void System::field(FieldTable &fieldtable, const Symbol &component, const Probe &probe, double scale_field)
{
    this->_pvt_forceFields().field(fieldtable, component, probe, scale_field);
}

/** Add the potentials acting on the molecules in the potential table 'pottable'
    from this system onto this potential table, scaled by the optionally
    supplied 'scale_potential' */
void System::potential(PotentialTable &pottable, double scale_potential)
{
    this->_pvt_forceFields().potential(pottable, scale_potential);
}

/** Add the potentials acting on the molecules in the potential table 'pottable'
    from the component of this system identified by 'component' onto
    this potential table, scaled by the optionally supplied 'scale_potential' */
void System::potential(PotentialTable &pottable, const Symbol &component, double scale_potential)
{
    this->_pvt_forceFields().potential(pottable, component, scale_potential);
}

/** Add the potentials acting on the molecules in the potential table 'pottable'
    from this system onto this potential table, scaled by the optionally
    supplied 'scale_potential' */
void System::potential(PotentialTable &pottable, const Probe &probe, double scale_potential)
{
    this->_pvt_forceFields().potential(pottable, probe, scale_potential);
}

/** Add the potentials acting on the molecules in the potential table 'pottable'
    from the component of this system identified by 'component' onto
    this potential table, scaled by the optionally supplied 'scale_potential' */
void System::potential(PotentialTable &pottable, const Symbol &component, const Probe &probe, double scale_potential)
{
    this->_pvt_forceFields().potential(pottable, component, probe, scale_potential);
}

/** Internal function used to copy a shared property to all of the
    contained molecules
*/
void System::_updateSharedProperty(const QString &name, const Property &value)
{
    if (this->needsAccepting())
    {
        this->accept();
    }

    Molecules mols_to_update;

    for (auto mol : this->molecules())
    {
        const auto &moldata = mol.constData();

        bool needs_updating = false;

        if (moldata.hasProperty(name))
        {
            needs_updating = not(value.equals(moldata.property(name)));
        }
        else
            needs_updating = true;

        if (needs_updating)
        {
            auto data = MoleculeData(moldata);
            data.setProperty(name, value);
            mol.update(data);

            mols_to_update.add(mol);
        }
    }

    if (not mols_to_update.isEmpty())
    {
        Delta delta(*this, true);

        // this ensures that only a single copy of System is used - prevents
        // unnecessary copying
        this->operator=(System());
        delta.update(mols_to_update);
        this->operator=(delta.apply());

        if (this->needsAccepting())
        {
            delta = Delta();
            this->accept();
        }
    }
}

/** Internal function used to set the value of the property, which has
 *  an extra "update_shared" flag to optionally update (or not) the
 *  value of any shared properties in the molecules.
 */
void System::_setProperty(const QString &name, const Property &value,
                          bool update_shared)
{
    this->_pvt_forceFields().setProperty(name, value);

    if (update_shared and this->shared_properties.contains(name))
    {
        // also update this property in all of the molecules
        this->_updateSharedProperty(name, value);
    }

    sysversion.incrementMajor();
    this->applyConstraints();
}

/** Set the value of the property called 'name' to the value 'value' in
    all forcefields that have this property

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setProperty(const QString &name, const Property &value)
{
    this->_setProperty(name, value, true);
}

/** Set the value of the property called 'name' in the forcefields identified
    by the ID 'ffid' to the value 'value'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setProperty(const FFID &ffid, const QString &name, const Property &value)
{
    this->_pvt_forceFields().setProperty(ffid, name, value);
    sysversion.incrementMajor();
    this->applyConstraints();
}

/** Remove the property with name 'name'. Note that this can only
    remove user-level properties - it cannot remove built-in properties
    of the system. This does nothing if there is no user-level
    property with this name */
void System::removeProperty(const QString &name)
{
    this->_pvt_forceFields().removeProperty(name);

    if (this->shared_properties.contains(name))
    {
        // remove this from all the contained molecules too!
        const Molecules mols = this->molecules();

        Molecules updated_mols = mols;

        for (const auto &mol : mols)
        {
            MoleculeData moldata(mol.data());
            moldata.removeProperty(name);

            updated_mols.update(moldata);
        }

        this->update(updated_mols);
    }

    sysversion.incrementMajor();
    this->applyConstraints();
}

/** Return whether or not the property 'name' exists and is a compound
    property (either a link or a combined property) */
bool System::isCompoundProperty(const QString &name) const
{
    return this->_pvt_forceFields().isCompoundProperty(name);
}

/** Return whether or not the property 'name' exists and is a user
    supplied property (either a compound property or an extra
    System property) */
bool System::isUserProperty(const QString &name) const
{
    return this->_pvt_forceFields().isUserProperty(name);
}

/** Return whether or not the property 'name' exists and is a builtin
    property of one of the forcefields in this System */
bool System::isBuiltinProperty(const QString &name) const
{
    return this->_pvt_forceFields().isBuiltinProperty(name);
}

/** Return the raw compound property with name 'name' - this returns
    the property representing the link, or the combined property,
    and raises an exception if a compound property with this name
    does not exist

    \throw SireBase::missing_property
*/
const Property &System::compoundProperty(const QString &name) const
{
    return this->_pvt_forceFields().compoundProperty(name);
}

/** Return the user-supplied property at 'name'. This raises an
    exception if there is no user-supplied property with this name

    \throw SireBase::missing_property
*/
const Property &System::userProperty(const QString &name) const
{
    return this->_pvt_forceFields().userProperty(name);
}

/** Return the built-in property at 'name'. This will by-pass any
    user-supplied property with this name, and will raise an
    exception if there is no built-in property with this name

    \throw SireBase::missing_property
*/
const Property &System::builtinProperty(const QString &name) const
{
    return this->_pvt_forceFields().builtinProperty(name);
}

/** Internal function called to apply all of the constraints to the system */
void System::applyAllConstraints()
{
    Delta delta(*this);
    this->operator=(delta.apply());
}

/** Apply the system (and molecule) constraints */
void System::applyConstraints()
{
    this->applyAllConstraints();
}

/** Return whether or not the constraints are satisfied */
bool System::constraintsSatisfied() const
{
    return cons.areSatisfied(*this);
}

/** Return the values of the property called 'name' in all of the
    forcefields that contain this property

    \throw SireBase::missing_property
    \throw SireBase::duplicate_property
*/
const Property &System::property(const PropertyName &name) const
{
    if (name.hasSource())
    {
        return this->_pvt_forceFields().property(name.source());
    }
    else
    {
        return name.value();
    }
}

/** Return the value of the property 'name' in the forcefield identified
    by the ID 'ffid'

    \throw SireBase::duplicate_property
    \throw SireFF::missing_forcefield
    \throw SireBase::missing_property
    \throw SireError::invalid_index
*/
const Property &System::property(const FFID &ffid, const PropertyName &name) const
{
    if (name.hasSource())
        return this->_pvt_forceFields().property(ffid, name.source());
    else
        return name.value();
}

/** Return the names of all of the properties of this system */
QStringList System::propertyKeys() const
{
    return this->_pvt_forceFields().propertyKeys();
}

/** Return the names of all of the properties of the forcefields in
    this system that match the ID 'ffid'

    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
QStringList System::propertyKeys(const FFID &ffid) const
{
    return this->_pvt_forceFields().propertyKeys(ffid);
}

/** Return whether or not any of the forcefields contain a property called 'name' */
bool System::containsProperty(const QString &name) const
{
    return this->_pvt_forceFields().containsProperty(name);
}

/** Return whether or not any of the forcefields identified by the ID 'ffid'
    contain a property called 'name' */
bool System::containsProperty(const FFID &ffid, const QString &name) const
{
    return this->_pvt_forceFields().containsProperty(ffid, name);
}

/** Return whether or not any of the forcefields contain a property called 'name' */
bool System::containsProperty(const PropertyName &name) const
{
    return this->_pvt_forceFields().containsProperty(name);
}

/** Return whether or not any of the forcefields identified by the ID 'ffid'
    contain a property called 'name' */
bool System::containsProperty(const FFID &ffid, const PropertyName &name) const
{
    return this->_pvt_forceFields().containsProperty(ffid, name);
}

/** Return the values of all of the properties of this system

    \throw SireBase::duplicate_property
*/
Properties System::properties() const
{
    return this->_pvt_forceFields().properties();
}

/** Return the values of all of the properties of this system that
    are in the forcefields that match the ID 'ffid'

    \throw SireBase::duplicate_property
    \throw SireFF::missing_forcefield
    \throw SireError::invalid_index
*/
Properties System::properties(const FFID &ffid) const
{
    return this->_pvt_forceFields().properties(ffid);
}

/** Return the values of all user-level properties of this
    system */
Properties System::userProperties() const
{
    return this->_pvt_forceFields().userProperties();
}

/** Return the values of all built-in properties of this system */
Properties System::builtinProperties() const
{
    return this->_pvt_forceFields().builtinProperties();
}

/** Specify that the property called 'name' should be a shared property */
void System::addSharedProperty(const QString &name)
{
    if (this->shared_properties.contains(name))
    {
        return;
    }

    this->shared_properties.append(name);

    if (this->containsProperty(name))
    {
        this->_updateSharedProperty(name, this->property(name));
    }
}

/** Specify that the passed property called 'name' with value 'value'
 *  should be a shared property
 */
void System::addSharedProperty(const QString &name, const Property &value)
{
    if (this->shared_properties.contains(name))
    {
        this->setProperty(name, value);
        return;
    }

    this->shared_properties.append(name);
    this->_updateSharedProperty(name, value);
    this->_setProperty(name, value, false);
}

/** Set the shared property called 'name' to the value 'value' */
void System::setSharedProperty(const QString &name, const Property &value)
{
    if (not this->shared_properties.contains(name))
    {
        throw SireError::invalid_key(QObject::tr(
                                         "There is no shared property called '%1'. The shared properties "
                                         "are: %2")
                                         .arg(name)
                                         .arg(Sire::toString(this->shared_properties)),
                                     CODELOC);
    }

    this->setProperty(name, value);
}

/** Remove the specified shared property from this system. This will
 *  remove the property, and will also remove its value from all of the
 *  contained molecules.
 */
void System::removeSharedProperty(const QString &name)
{
    if (not this->shared_properties.contains(name))
    {
        return;
    }

    this->removeProperty(name);
    this->shared_properties.removeAll(name);
}

void System::removeAllSharedProperties()
{
    auto to_remove = this->shared_properties;

    for (const auto &name : to_remove)
    {
        this->removeSharedProperty(name);
    }
}

Properties System::sharedProperties() const
{
    Properties props;

    for (const auto &name : this->shared_properties)
    {
        if (this->containsProperty(name))
        {
            props.setProperty(name, this->property(name));
        }
        else
        {
            props.setProperty(name, NullProperty());
        }
    }

    return props;
}

/** Return the list of all monitors of this system */
const SystemMonitors &System::monitors() const
{
    return sysmonitors;
}

/** Return an array of all of the forcefields in this system */
const ForceFields &System::forceFields() const
{
    return this->_pvt_forceFields();
}

/** Return all of the extra non-forcefield molecule groups in this system */
const MoleculeGroups &System::extraGroups() const
{
    return this->_pvt_moleculeGroups();
}

/** Return all of the contraints that are applied to the system */
const Constraints &System::constraints() const
{
    return cons;
}

/** Completely clear all statistics held in the monitors */
void System::clearStatistics()
{
    sysmonitors.clearStatistics();
}

/** Clear the statistics of the monitors that match the ID 'monid'.
    This does nothing if there are no matching monitors */
void System::clearStatistics(const MonitorID &monid)
{
    sysmonitors.clearStatistics(monid);
}

/** Tell all of the forcefields that they will need to recalculate
    their energies from scratch. This can speed up calculations where
    you know that the majority (or all) of the molecules will be
    changing */
void System::mustNowRecalculateFromScratch()
{
    this->_pvt_forceFields().mustNowRecalculateFromScratch();
}

/** Return whether or not any part of the forcefield is using temporary
    workspaces that need to be accepted */
bool System::needsAccepting() const
{
    return this->_pvt_forceFields().needsAccepting() or this->_pvt_moleculeGroups().needsAccepting();
}

/** Tell all of the forcefields that the last move was accepted. This allows
    any cacheing or use of temporary workspaces to be committed */
void System::accept()
{
    this->_pvt_forceFields().accept();
    this->_pvt_moleculeGroups().accept();
}

/** Return whether or not any of the forcefields are dirty */
bool System::isDirty() const
{
    return this->_pvt_forceFields().isDirty();
}

/** Return whether or not all of the forcefields are clean */
bool System::isClean() const
{
    return this->_pvt_forceFields().isClean();
}

/** Add a system monitor 'monitor', identified by the name 'name', which
    will be updated every 'frequency' steps.

    \throw SireSystem::duplicate_monitor
*/
void System::add(const QString &name, const SystemMonitor &monitor, int frequency)
{
    sysmonitors.add(name, monitor, frequency);
    sysversion.incrementMajor();
}

/** Add the monitors in 'monitors' to this system

    \throw SireSystem::duplicate_monitor
*/
void System::add(const SystemMonitors &monitors)
{
    sysmonitors.add(monitors);
    sysversion.incrementMajor();
}

/** Add the monitors in 'monitors', setting the frequency of the
    new monitors to 'frequency'

    \throw SireSystem::duplicate_monitor
*/
void System::add(const SystemMonitors &monitors, int frequency)
{
    sysmonitors.add(monitors, frequency);
    sysversion.incrementMajor();
}

/** Set the monitors of this system to 'monitors' */
void System::setMonitors(const SystemMonitors &monitors)
{
    if (sysmonitors != monitors)
    {
        sysmonitors = monitors;
        sysversion.incrementMajor();
    }
}

/** Set the monitors of the system to 'monitors', and reset the
    frequency of all of the monitors so that they are triggered
    every 'frequency' steps */
void System::setMonitors(const SystemMonitors &monitors, int frequency)
{
    SystemMonitors new_monitors(monitors);
    new_monitors.setAllFrequency(frequency);
    this->setMonitors(new_monitors);
}

/** Add the forcefield 'forcefield' to this system. This will raise
    an exception if this forcefield (or one with the same name)
    is already present in this set. Note that if the added
    forcefield will be updated to contain the versions of
    any molecules that are already present in any of the
    other forcefields.

    \throw SireFF::duplicate_forcefield
    \throw SireMol::duplicate_group
*/
void System::add(const FF &forcefield)
{
    FFPtr ff(forcefield);
    ff.edit().update(this->matchToExistingVersion(forcefield.molecules()));

    SaveState old_state = SaveState::save(*this);

    try
    {
        this->_pvt_forceFields().add(ff.read());
        this->rebuildIndex();
        sysversion.incrementMajor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the molecule group 'molgroup' to this system. If this is
    a molecule group that is part of a forcefield, then the entire
    forcefield will be added to this system. This will raise
    an exception if this molecule group is already present in
    this system. Note that the added molecule group will be
    updated to contain the version of the any molecules that
    are already present in this system

    \throw SireFF::duplicate_forcefield
    \throw SireMol::duplicate_group
*/
void System::add(const MoleculeGroup &molgroup)
{
    if (molgroup.isA<FFMolGroup>())
    {
        this->add(molgroup.asA<FFMolGroup>().forceField());
    }
    else if (molgroup.isA<SireFF::detail::FFMolGroupPvt>())
    {
        this->add(molgroup.asA<SireFF::detail::FFMolGroupPvt>().forceField());
    }
    else
    {
        MolGroupPtr mgroup(molgroup);
        mgroup.edit().update(this->matchToExistingVersion(molgroup.molecules()));

        SaveState old_state = SaveState::save(*this);

        try
        {
            this->_pvt_moleculeGroups().add(mgroup);
            this->rebuildIndex();
            sysversion.incrementMajor();

            this->applyAllConstraints();
        }
        catch (...)
        {
            old_state.restore(*this);
            throw;
        }
    }
}

/** Add the passed constraint to the system */
void System::add(const Constraints &constraints)
{
    SaveState old_state = SaveState::save(*this);

    try
    {
        int nconstraints = cons.count();

        cons.add(constraints);

        if (cons.count() != nconstraints)
        {
            sysversion.incrementMajor();
            this->applyAllConstraints();
        }
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the passed constraint to the system */
void System::add(const Constraint &constraint)
{
    this->add(Constraints(constraint));
}

/** Set the constraints for the system equal to 'constraints' */
void System::setConstraints(const Constraints &constraints)
{
    if (cons == constraints)
        return;

    SaveState old_state = SaveState::save(*this);

    try
    {
        cons = constraints;
        sysversion.incrementMajor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Remove all monitors that match the ID 'monid'

    \throw SireSystem::missing_monitor
    \throw SireError::invalid_index
*/
void System::remove(const MonitorID &monid)
{
    sysmonitors.remove(monid);
    sysversion.incrementMajor();
}

/** Remove the forcefield(s) that match the ID 'ffid'

    \throw SireError::invalid_index
*/
void System::remove(const FFID &ffid)
{
    SaveState old_state = SaveState::save(*this);

    try
    {
        this->_pvt_forceFields().remove(ffid);
        this->rebuildIndex();
        sysversion.incrementMajor();
        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Remove the forcefield 'ff'. Note that this removes the forcefield
    in this system that has the same name as 'ff'

    \throw SireFF::missing_forcefield
*/
void System::remove(const FF &ff)
{
    this->remove(ff.name());
}

/** Remove the molecule group(s) that match the ID 'mgid'.
    Note that you can't remove molecule groups that are part
    of a forcefield

    \throw SireMol::missing_group
    \throw SireError::invalid_index
    \throw SireError::invalid_arg
*/
bool System::remove(const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    try
    {
        bool removed = false;

        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
            {
                const MoleculeGroup &molgroup = this->_pvt_moleculeGroup(mgnum);
                BOOST_ASSERT(molgroup.isA<SireFF::detail::FFMolGroupPvt>());

                const FF &ff = molgroup.asA<SireFF::detail::FFMolGroupPvt>().forceField();

                throw SireError::invalid_arg(QObject::tr("You cannot remove the molecule group with number %1 "
                                                         "(that matches the ID %2) as it is part of the forcefield "
                                                         "%3 in this system.")
                                                 .arg(mgnum)
                                                 .arg(mgid.toString())
                                                 .arg(ff.name()),
                                             CODELOC);
            }

            this->_pvt_moleculeGroups().remove(mgnum);
            this->removeFromIndex(mgnum);
            mgroups_by_num.remove(mgnum);
            removed = true;
        }

        if (removed)
        {
            sysversion.incrementMajor();
            this->applyAllConstraints();

            return true;
        }
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return false;
}

/** Remove the molecules contained in the molecule group 'molgroup'.
    This doesn't remove the molecule group itself though. If you
    want to remove the molecule group, use System::remove(molgroup.number())

    \throw SireMol::missing_group
*/
bool System::remove(const MoleculeGroup &molgroup)
{
    return this->remove(molgroup.molecules());
}

/** Remove all molecules from this system that match the ID 'molid'

    \throw SireMol::missing_molecule
*/
bool System::remove(const MolID &molid)
{
    QList<MolNum> molnums = molid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MolNum molnum, molnums)
        {
            foreach (MGNum mgnum, this->groupsContaining(molnum))
            {
                molgroups[mgroups_by_num.value(mgnum)].edit().remove(molnum, mgnum);
                mols_removed = true;
            }

            this->removeFromIndex(molnum);
        }

        sysversion.incrementMinor();
        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove the constraints in 'constraints' from the list of constraints
    that are applied to this system */
void System::remove(const Constraints &constraints)
{
    int nconstraints = cons.count();
    cons.remove(constraints);

    if (cons.count() != nconstraints)
    {
        sysversion.incrementMajor();
        this->applyAllConstraints();
    }
}

/** Remove the passed constraint from this system */
void System::remove(const Constraint &constraint)
{
    this->remove(Constraints(constraint));
}

/** Completely remove all molecules from this system */
bool System::removeAllMolecules()
{
    return this->removeAll();
}

/** Completely remove all non-forcefield molecule groups
    from this system */
void System::removeAllMoleculeGroups()
{
    SaveState old_state = SaveState::save(*this);

    try
    {
        QList<MGNum> mgnums = this->_pvt_constMoleculeGroups().mgNums();

        if (not mgnums.isEmpty())
        {
            this->_pvt_moleculeGroups().remove(IDOrSet<MGID>(mgnums));

            foreach (MGNum mgnum, mgnums)
            {
                this->removeFromIndex(mgnum);
                mgroups_by_num.remove(mgnum);
            }
        }

        sysversion.incrementMajor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Completely remove all monitors from this system */
void System::removeAllMonitors()
{
    if (not sysmonitors.isEmpty())
    {
        sysmonitors.removeAll();
        sysversion.incrementMajor();
    }
}

/** Completely remove all of the forcefields (and their contained
    molecule groups) from this system */
void System::removeAllForceFields()
{
    if (this->nForceFields() == 0)
        return;

    SaveState old_state = SaveState::save(*this);

    try
    {
        QList<MGNum> mgnums = this->_pvt_forceFields().mgNums();
        this->_pvt_forceFields().removeAllForceFields();

        foreach (MGNum mgnum, mgnums)
        {
            this->removeFromIndex(mgnum);
            mgroups_by_num.remove(mgnum);
        }

        sysversion.incrementMajor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Remove all constraints from this system */
void System::removeAllConstraints()
{
    if (cons.count() == 0)
        return;

    cons = Constraints();

    sysversion.incrementMajor();
}

/** Return the molecule group with number 'mgnum'

    \throw SireMol::missing_group
*/
const MoleculeGroup &System::at(MGNum mgnum) const
{
    return this->_pvt_moleculeGroup(mgnum);
}

/** Update the passed molecule so that it contains the current
 *  versions of the shared properties.
 */
void System::_updateSharedProperties(MoleculeData &data) const
{
    if (shared_properties.isEmpty())
        return;

    for (const auto &name : this->shared_properties)
    {
        if (this->containsProperty(name))
        {
            const auto &prop = this->property(name);

            if (data.hasProperty(name))
            {
                if (not prop.equals(data.property(name)))
                {
                    data.setProperty(name, prop);
                }
            }
            else
            {
                data.setProperty(name, prop);
            }
        }
    }
}

/** Update the passed molecule so that it contains the current
 *  versions of the shared properties.
 */
void System::_updateSharedProperties(MoleculeView &view) const
{
    if (shared_properties.isEmpty())
        return;

    const auto &data = view.constData();

    bool needs_update = false;

    for (const auto &name : this->shared_properties)
    {
        if (this->containsProperty(name))
        {
            const auto &prop = this->property(name);

            if (data.hasProperty(name))
            {
                if (not prop.equals(data.property(name)))
                {
                    needs_update = true;
                    break;
                }
            }
            else
            {
                needs_update = true;
                break;
            }
        }
    }

    if (needs_update)
    {
        MoleculeData moldata(data);
        this->_updateSharedProperties(moldata);
        view.update(moldata);
    }
}

/** Add the molecule viewed in 'molview' to the molecule groups
    identified by the ID 'mgid'. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to. The version of the molecule
    already present in this system is used, if such a molecule exists.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const MoleculeView &molview, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    PartialMolecule view(molview);

    if (this->contains(molview.data().number()))
    {
        // Update so that the added moleucle has the same version and data
        // as the current copy in this system.
        // This will make sure that any shared properties are up to data
        view.update(this->matchToExistingVersion(molview.data()));
    }
    else
    {
        // Update the molecule so that it contains the shared properties
        this->_updateSharedProperties(view);
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().add(view, mgnum, map);
            else
                this->_pvt_moleculeGroups().add(view, mgnum);

            this->addToIndex(mgnum, view.data().number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the views of the molecule in 'molviews' to the molecule groups
    identified by the ID 'mgid'. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to. The version of the molecule
    already present in this system is used, if such a molecule exists.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const ViewsOfMol &molviews, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    ViewsOfMol views(molviews);

    if (this->contains(views.data().number()))
    {
        // Update so that the added moleucle has the same version and data
        // as the current copy in this system.
        // This will make sure that any shared properties are up to data
        views.update(this->matchToExistingVersion(molviews.data()));
    }
    else
    {
        // Update the molecule so that it contains the shared properties
        this->_updateSharedProperties(views);
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().add(views, mgnum, map);
            else
                this->_pvt_moleculeGroups().add(views, mgnum);

            this->addToIndex(mgnum, views.number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the molecules viewed in 'molecules' to the molecule groups
    identified by the ID 'mgid'. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to. The version of the molecule
    already present in this system is used, if such a molecule exists.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const Molecules &molecules, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    Molecules mols = molecules;

    // update the molecules to match what is in this system, including
    // the shared properties for new molecules
    for (auto mol : molecules)
    {
        if (this->contains(mol.data().number()))
        {
            mol.update(this->matchToExistingVersion(mol.data()));
        }
        else
        {
            this->_updateSharedProperties(mol);
        }

        mols.update(mol.data());
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().add(mols, mgnum, map);
            else
                this->_pvt_moleculeGroups().add(mols, mgnum);

            this->addToIndex(mgnum, mols.molNums());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the molecules in the molecule group 'molgroup' to the molecule groups
    identified by the ID 'mgid'. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to. The version of the molecule
    already present in this system is used, if such a molecule exists.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const MoleculeGroup &molgroup, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    MolGroupPtr group(molgroup);

    for (auto mol : molgroup.molecules())
    {
        if (this->contains(mol.data().number()))
        {
            mol.update(this->matchToExistingVersion(mol.data()));
        }
        else
        {
            this->_updateSharedProperties(mol);
        }

        group.edit().update(mol.data());
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().add(group, mgnum, map);
            else
                this->_pvt_moleculeGroups().add(group, mgnum);

            this->addToIndex(mgnum, molgroup.molNums().toList());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the view of the molecule in 'molview' to the groups
    identified by 'mgid'. This only adds the view to a group
    if it doesn't already exist in the group. The version
    of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const MoleculeView &molview, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    PartialMolecule view(molview);

    if (this->contains(molview.data().number()))
    {
        view.update(this->matchToExistingVersion(molview.data()));
    }
    else
    {
        this->_updateSharedProperties(view);
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().addIfUnique(view, mgnum, map);
            else
                this->_pvt_moleculeGroups().addIfUnique(view, mgnum);

            this->addToIndex(mgnum, view.data().number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the views of the molecule in 'molviews' to the groups
    identified by 'mgid'. This only adds the view to a group
    if it doesn't already exist in the group. The version
    of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    ViewsOfMol views(molviews);

    if (this->contains(molviews.data().number()))
    {
        views.update(this->matchToExistingVersion(molviews.data()));
    }
    else
    {
        this->_updateSharedProperties(views);
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().addIfUnique(views, mgnum, map);
            else
                this->_pvt_moleculeGroups().addIfUnique(views, mgnum);

            this->addToIndex(mgnum, views.number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the views of the molecules in 'molecules' to the groups
    identified by 'mgid'. This only adds the view to a group
    if it doesn't already exist in the group. The version
    of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const Molecules &molecules, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    Molecules mols = molecules;

    for (auto mol : molecules)
    {
        if (this->contains(mol.data().number()))
        {
            mol.update(this->matchToExistingVersion(mol.data()));
        }
        else
        {
            this->_updateSharedProperties(mol);
        }

        mols.update(mol);
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().addIfUnique(mols, mgnum, map);
            else
                this->_pvt_moleculeGroups().addIfUnique(mols, mgnum);

            this->addToIndex(mgnum, mols.molNums());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Add the view of the molecules in the group 'molgroup' to the groups
    identified by 'mgid'. This only adds the view to a group
    if it doesn't already exist in the group. The version
    of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    MolGroupPtr group(molgroup);

    for (auto mol : molgroup.molecules())
    {
        if (this->contains(mol.data().number()))
        {
            mol.update(this->matchToExistingVersion(mol.data()));
        }
        else
        {
            this->_updateSharedProperties(mol);
        }

        group.edit().update(mol.data());
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().addIfUnique(group, mgnum, map);
            else
                this->_pvt_moleculeGroups().addIfUnique(group, mgnum);

            this->addToIndex(mgnum, molgroup.molNums().toList());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Convenient overload of System::add that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const MoleculeView &molview, const MGID &mgid)
{
    this->add(molview, mgid, PropertyMap());
}

/** Convenient overload of System::add that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->add(molviews, mgid, PropertyMap());
}

/** Convenient overload of System::add that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const Molecules &molecules, const MGID &mgid)
{
    this->add(molecules, mgid, PropertyMap());
}

/** Convenient overload of System::add that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::add(const MoleculeGroup &molgroup, const MGID &mgid)
{
    this->add(molgroup, mgid, PropertyMap());
}

/** Convenient overload of System::addIfUnique that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const MoleculeView &molview, const MGID &mgid)
{
    this->addIfUnique(molview, mgid, PropertyMap());
}

/** Convenient overload of System::addIfUnique that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const ViewsOfMol &molviews, const MGID &mgid)
{
    this->addIfUnique(molviews, mgid, PropertyMap());
}

/** Convenient overload of System::addIfUnique that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const Molecules &molecules, const MGID &mgid)
{
    this->addIfUnique(molecules, mgid, PropertyMap());
}

/** Convenient overload of System::addIfUnique that uses the default locations
    to find any necessary properties.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::addIfUnique(const MoleculeGroup &molgroup, const MGID &mgid)
{
    this->addIfUnique(molgroup, mgid, PropertyMap());
}

/** Remove all molecules from the molecule groups identified by the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::removeAll(const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (this->_pvt_moleculeGroups(mgnum).removeAll(mgnum))
            {
                this->clearIndex(mgnum);
                mols_removed = true;
            }
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove the view 'molview' from the specified groups in this
    forcefield. Note that this only removes the specific view
    (and indeed only the first copy of this view if there
    are duplicates) - it does not remove the atoms in this
    view from all of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::remove(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            MolGroupsBase &molgroups = this->_pvt_moleculeGroups(mgnum);

            bool mol_removed = molgroups.remove(molview, mgnum);

            mols_removed = mols_removed or mol_removed;

            if (not molgroups.at(mgnum).contains(molview.data().number()))
                this->removeFromIndex(mgnum, molview.data().number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove the views in 'molviews' from the specified groups in this
    forcefield. Note that this only removes the specific views
    (and indeed only the first copy of this view if there
    are duplicates) - it does not remove the atoms in this
    view from all of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::remove(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            MolGroupsBase &molgroups = this->_pvt_moleculeGroups(mgnum);

            bool mol_removed = molgroups.remove(molviews, mgnum);
            mols_removed = mols_removed or mol_removed;

            if (not molgroups.at(mgnum).contains(molviews.number()))
                this->removeFromIndex(mgnum, molviews.number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove them molecules in 'molecules' from the specified groups in this
    forcefield. Note that this only removes the specific views
    (and indeed only the first copy of this view if there
    are duplicates) - it does not remove the atoms in this
    view from all of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::remove(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            MolGroupsBase &molgroups = this->_pvt_moleculeGroups(mgnum);

            bool mol_removed = molgroups.remove(molecules, mgnum);
            mols_removed = mols_removed or mol_removed;

            const MoleculeGroup &molgroup = molgroups.at(mgnum);

            for (Molecules::const_iterator it = molecules.constBegin(); it != molecules.constEnd(); ++it)
            {
                if (not molgroup.contains(it->number()))
                    this->removeFromIndex(mgnum, it->number());
            }
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove the views in the molecule group 'molgroup' from the specified
    groups in this forcefield. Note that this only removes the specific views
    (and indeed only the first copy of this view if there
    are duplicates) - it does not remove the atoms in this
    view from all of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::remove(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->remove(molgroup.molecules(), mgid);
}

/** Remove the all copies of the view in 'molview' from the specified
    groups in this forcefield. Note that this only removes the specific views
    - it does not remove the atoms in this view from all of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::removeAll(const MoleculeView &molview, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            MolGroupsBase &molgroups = this->_pvt_moleculeGroups(mgnum);

            bool mol_removed = molgroups.removeAll(molview, mgnum);
            mols_removed = mols_removed or mol_removed;

            if (not molgroups.at(mgnum).contains(molview.data().number()))
                this->removeFromIndex(mgnum, molview.data().number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove the all copies of the views in 'molviews' from the specified
    groups in this forcefield. Note that this only removes the specific views
    - it does not remove the atoms in this view from all of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::removeAll(const ViewsOfMol &molviews, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            MolGroupsBase &molgroups = this->_pvt_moleculeGroups(mgnum);

            bool mol_removed = molgroups.removeAll(molviews, mgnum);
            mols_removed = mols_removed or mol_removed;

            if (not molgroups.at(mgnum).contains(molviews.number()))
                this->removeFromIndex(mgnum, molviews.number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove the all copies of the molecules in 'molecules' from the specified
    groups in this forcefield. Note that this only removes the specific views
    - it does not remove the atoms in this view from all of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::removeAll(const Molecules &molecules, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            MolGroupsBase &molgroups = this->_pvt_moleculeGroups(mgnum);

            bool mol_removed = molgroups.removeAll(molecules, mgnum);
            mols_removed = mols_removed or mol_removed;

            const MoleculeGroup &molgroup = molgroups.at(mgnum);

            for (Molecules::const_iterator it = molecules.constBegin(); it != molecules.constEnd(); ++it)
            {
                if (not molgroup.contains(it->number()))
                    this->removeFromIndex(mgnum, it->number());
            }
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove the all copies of the molecules in the molecule group 'molgroup'
    from the specified groups in this forcefield. Note that this only removes
    the specific views - it does not remove the atoms in this view from all
    of the other views

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::removeAll(const MoleculeGroup &molgroup, const MGID &mgid)
{
    return this->removeAll(molgroup.molecules(), mgid);
}

/** Remove all views of the molecule with number 'molnum' from the molecule
    groups identified by 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::remove(MolNum molnum, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (this->_pvt_moleculeGroups(mgnum).remove(molnum, mgnum))
            {
                this->removeFromIndex(mgnum, molnum);
                mols_removed = true;
            }
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Remove all of the molecules whose numbers are in 'molnums' from
    all of the molecule groups identified by the ID 'mgid'

    \throw SireMol::missing_group
    \throw SireError::invalid_index
*/
bool System::remove(const QSet<MolNum> &molnums, const MGID &mgid)
{
    QList<MGNum> mgnums = mgid.map(*this);

    SaveState old_state = SaveState::save(*this);

    bool mols_removed = false;

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (this->_pvt_moleculeGroups(mgnum).remove(molnums, mgnum))
            {
                this->removeFromIndex(mgnum, molnums);
                mols_removed = true;
            }
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }

    return mols_removed;
}

/** Update this system so that it uses the version of the molecule
    available in 'moldata'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::update(const MoleculeData &moldata, bool auto_commit)
{
    MoleculeData data = moldata;
    this->_updateSharedProperties(data);

    Delta delta(*this, auto_commit);

    // this ensures that only a single copy of System is used - prevents
    // unnecessary copying
    this->operator=(System());
    delta.update(data);
    this->operator=(delta.apply());

    if (auto_commit and this->needsAccepting())
    {
        delta = Delta();
        this->accept();
    }
}

void System::update(const MoleculeView &molview, bool auto_commit)
{
    this->update(molview.data(), auto_commit);
}

/** Update this system so that it uses the same version of the molecules
    present in 'molecules'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::update(const Molecules &molecules, bool auto_commit)
{
    Molecules mols = molecules;

    if (not shared_properties.isEmpty())
    {
        for (auto mol : molecules)
        {
            this->_updateSharedProperties(mol);
            mols.update(mol.data());
        }
    }

    Delta delta(*this, auto_commit);

    // this ensures that only a single copy of System is used - prevents
    // unnecessary copying
    this->operator=(System());
    delta.update(mols);
    this->operator=(delta.apply());

    if (auto_commit and this->needsAccepting())
    {
        delta = Delta();
        this->accept();
    }
}

void System::update(const SelectorMol &molecules, bool auto_commit)
{
    this->update(molecules.toMolecules(), auto_commit);
}

/** Update this system so that it uses the same version of the molecules
    present in the molecule group 'molgroup'

    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::update(const MoleculeGroup &molgroup, bool auto_commit)
{
    this->update(molgroup.molecules(), auto_commit);
}

/** Set the contents of the molecule group(s) identified by the ID 'mgid'
    so that they contain just the view of the molecule in 'molview'.
    The version of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const MoleculeView &molview, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    PartialMolecule view(molview);

    if (this->contains(view.data().number()))
    {
        view.update(this->matchToExistingVersion(molview.data()));
    }
    else
    {
        this->_updateSharedProperties(view);
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().setContents(mgnum, view, map);
            else
                this->_pvt_moleculeGroups(mgnum).setContents(mgnum, view);

            this->clearIndex(mgnum);
            this->addToIndex(mgnum, view.data().number());
        }

        sysversion.incrementMinor();
        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Set the contents of the molecule group(s) identified by the ID 'mgid'
    so that they contain just the views of the molecule in 'molviews'.
    The version of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const ViewsOfMol &molviews, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    ViewsOfMol views(molviews);

    if (this->contains(views.data().number()))
    {
        views.update(this->matchToExistingVersion(molviews.data()));
    }
    else
    {
        this->_updateSharedProperties(views);
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().setContents(mgnum, views, map);
            else
                this->_pvt_moleculeGroups(mgnum).setContents(mgnum, views);

            this->clearIndex(mgnum);
            this->addToIndex(mgnum, views.number());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Set the contents of the molecule group(s) identified by the ID 'mgid'
    so that they contain just the views of the molecules in 'molecules'.
    The version of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const Molecules &molecules, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    Molecules mols = molecules;

    for (auto mol : molecules)
    {
        if (this->contains(mol.data().number()))
        {
            mol.update(this->matchToExistingVersion(mol.data()));
        }
        else
        {
            this->_updateSharedProperties(mol);
        }

        mols.update(mol.data());
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().setContents(mgnum, mols, map);
            else
                this->_pvt_moleculeGroups(mgnum).setContents(mgnum, mols);

            this->clearIndex(mgnum);
            this->addToIndex(mgnum, mols.molNums());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Set the contents of the molecule group(s) identified by the ID 'mgid'
    so that they contain just the molecules in the group 'molgroup'.
    The version of the molecule already present in this set is used if
    such a molecule already exists. The supplied property map
    is used to find the properties required by any forcefields
    that this molecule may be added to.

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const MoleculeGroup &molgroup, const PropertyMap &map)
{
    QList<MGNum> mgnums = mgid.map(*this);

    MolGroupPtr group(molgroup);

    for (auto mol : molgroup.molecules())
    {
        if (this->contains(mol.data().number()))
        {
            mol.update(this->matchToExistingVersion(mol.data()));
        }
        else
        {
            this->_updateSharedProperties(mol);
        }

        group.edit().update(mol.data());
    }

    SaveState old_state = SaveState::save(*this);

    try
    {
        foreach (MGNum mgnum, mgnums)
        {
            if (mgroups_by_num.value(mgnum) == 0)
                this->_pvt_forceFields().setContents(mgnum, group, map);
            else
                this->_pvt_moleculeGroups(mgnum).setContents(mgnum, group);

            this->clearIndex(mgnum);
            this->addToIndex(mgnum, group->molNums().toList());
        }

        sysversion.incrementMinor();

        this->applyAllConstraints();
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

/** Convenient overload of System::setContents that uses the default
    property locations to find the properties required by the forcefields

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const MoleculeView &molview)
{
    this->setContents(mgid, molview, PropertyMap());
}

/** Convenient overload of System::setContents that uses the default
    property locations to find the properties required by the forcefields

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const ViewsOfMol &molviews)
{
    this->setContents(mgid, molviews, PropertyMap());
}

/** Convenient overload of System::setContents that uses the default
    property locations to find the properties required by the forcefields

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const Molecules &molecules)
{
    this->setContents(mgid, molecules, PropertyMap());
}

/** Convenient overload of System::setContents that uses the default
    property locations to find the properties required by the forcefields

    \throw SireMol::missing_group
    \throw SireBase::missing_property
    \throw SireError::invalid_index
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void System::setContents(const MGID &mgid, const MoleculeGroup &molgroup)
{
    this->setContents(mgid, molgroup, PropertyMap());
}

bool System::deltaUpdate(const Symbol &component, double value)
{
    this->_pvt_forceFields().setConstantComponent(component, value);
    ++subversion;
    return true;
}

bool System::deltaUpdate(const QString &property, const Property &value)
{
    this->_pvt_forceFields().setProperty(property, value);
    ++subversion;
    return true;
}

bool System::deltaUpdate(const QString &property, const FFID &ffid, const Property &value)
{
    this->_pvt_forceFields().setProperty(ffid, property, value);
    ++subversion;
    return true;
}

bool System::deltaUpdate(const QString &property, const QList<FFIdx> &ffidxs, const Property &value)
{
    if (ffidxs.isEmpty())
        return false;

    else if (ffidxs.count() == 1)
    {
        this->_pvt_forceFields().setProperty(ffidxs.at(0), property, value);
        ++subversion;
        return true;
    }
    else
    {
        SaveState old_state = SaveState::save(*this);

        try
        {
            foreach (const FFIdx &ffidx, ffidxs)
            {
                this->_pvt_forceFields().setProperty(ffidx, property, value);
            }

            ++subversion;
        }
        catch (...)
        {
            old_state.restore(*this);
            throw;
        }

        return true;
    }
}

bool System::deltaUpdate(const MoleculeData &moldata, bool auto_commit)
{
    bool in_molgroup = this->_pvt_constMoleculeGroups().contains(moldata.number());
    bool in_ffields = this->_pvt_constForceFields().contains(moldata.number());

    if (in_molgroup or in_ffields)
    {
        if (in_molgroup)
            this->_pvt_moleculeGroups().update(moldata, auto_commit);

        if (in_ffields)
            this->_pvt_forceFields().update(moldata, auto_commit);

        ++subversion;

        return true;
    }
    else
        return false;
}

QList<MolNum> System::deltaUpdate(const Molecules &molecules, bool auto_commit)
{
    if (molecules.isEmpty())
        return QList<MolNum>();

    else if (molecules.count() == 1)
    {
        QList<MolNum> molnums;

        if (this->deltaUpdate(molecules.constBegin()->data(), auto_commit))
            molnums.append(molecules.constBegin()->data().number());

        return molnums;
    }

    QList<MolNum> changed_mols;

    for (Molecules::const_iterator it = molecules.constBegin(); it != molecules.constEnd(); ++it)
    {
        if (this->contains(it.key()) and this->getMoleculeVersion(it.key()) != it.value().version())
        {
            changed_mols.append(it.key());
        }
    }

    if (changed_mols.isEmpty())
        return changed_mols;

    bool in_molgroup = this->_pvt_constMoleculeGroups().contains(changed_mols);
    bool in_ffields = this->_pvt_constForceFields().contains(changed_mols);

    if (in_molgroup or in_ffields)
    {
        if (in_ffields)
            this->_pvt_forceFields().update(molecules, auto_commit);

        if (in_molgroup)
            this->_pvt_moleculeGroups().update(molecules, auto_commit);

        ++subversion;

        return changed_mols;
    }
    else
        return QList<MolNum>();
}

void System::commitDelta(const Constraints &constraints, bool is_minor_change, bool is_major_change)
{
    SaveState old_state = SaveState::save(*this);

    try
    {
        cons = constraints;

        if (is_major_change)
            sysversion.incrementMajor();
        else if (is_minor_change)
            sysversion.incrementMinor();

        subversion = 0;

        cons.committed(*this);
    }
    catch (...)
    {
        old_state.restore(*this);
        throw;
    }
}

int System::nFrames() const
{
    return this->nFrames(PropertyMap());
}

int System::nFrames(const SireBase::PropertyMap &map) const
{
    return MolGroupsBase::nFrames(map);
}

void System::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

void System::loadFrame(int frame, const LazyEvaluator &evaluator)
{
    this->loadFrame(frame, evaluator, PropertyMap());
}

void System::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

void System::saveFrame()
{
    this->saveFrame(PropertyMap());
}

void System::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

void System::deleteAllFrames()
{
    this->deleteAllFrames(PropertyMap());
}

static SireUnits::Dimension::Time get_time_from_property(const Property &prop)
{
    try
    {
        if (prop.isA<TimeProperty>())
            return prop.asA<TimeProperty>().value();
        else
            return prop.asA<GeneralUnitProperty>().toUnit<SireUnits::Dimension::Time>();
    }
    catch (...)
    {
        return SireUnits::Dimension::Time(0);
    }
}

void System::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    LazyEvaluator evaluator;
    this->loadFrame(frame, evaluator, map);
}

void System::loadFrame(int frame, const LazyEvaluator &evaluator,
                       const SireBase::PropertyMap &map)
{
    this->accept();
    this->mustNowRecalculateFromScratch();
    MolGroupsBase::loadFrame(frame, evaluator, map);

    const auto space_property = map["space"];
    const auto time_property = map["time"];

    QString time_property_source("time");
    QString space_property_source("space");

    if (time_property.hasSource())
        time_property_source = QString(time_property.source());

    if (space_property.hasSource())
        space_property_source = QString(space_property.source());

    // we must get the space property and a time property
    SpacePtr old_space;
    SireUnits::Dimension::Time old_time;

    try
    {
        old_space = this->property(space_property_source).asA<Space>();
    }
    catch (...)
    {
        // no default space
        old_space = Cartesian();
    }

    try
    {
        old_time = get_time_from_property(this->property(time_property_source));
    }
    catch (const std::exception &e)
    {
        old_time = SireUnits::Dimension::Time(0);
    }

    SpacePtr new_space;
    SireUnits::Dimension::Time new_time;
    bool found_space = false;
    bool found_time = false;

    const auto groups = this->getGroups();

    if (space_property.hasValue())
    {
        new_space = space_property.value().asA<Space>();
        found_space = true;
    }

    if (time_property.hasValue())
    {
        new_time = get_time_from_property(time_property.value());
        found_time = true;
    }

    for (auto it = groups.constBegin(); it != groups.constEnd(); ++it)
    {
        for (const auto &mol : (*it)->molecules())
        {
            if (not found_space)
            {
                try
                {
                    new_space = mol.data().property(space_property_source).asA<Space>();
                    found_space = true;
                }
                catch (...)
                {
                }
            }

            if (not found_time)
            {
                try
                {
                    new_time = get_time_from_property(mol.data().property(time_property_source));
                    found_time = true;
                }
                catch (...)
                {
                }
            }

            if (found_space and found_time)
                break;
        }

        if (found_space and found_time)
            break;
    }

    if (found_space and (new_time != old_time) and
        this->shared_properties.contains(time_property_source))
    {
        this->_setProperty(time_property_source,
                           GeneralUnitProperty(new_time), false);
    }

    if (found_space and (not old_space.read().equals(*new_space)) and
        this->shared_properties.contains(space_property_source))
    {
        this->_setProperty(space_property_source, new_space, false);
    }
}

void System::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    if (frame == this->nFrames(map))
    {
        this->saveFrame(map);
        return;
    }

    this->accept();
    this->mustNowRecalculateFromScratch();
    MolGroupsBase::saveFrame(frame, map);
    system_trajectory = 0;
}

void System::saveFrame(const SireBase::PropertyMap &map)
{
    auto traj_prop = map["trajectory"];

    if (not traj_prop.hasSource())
        return;

    this->accept();
    this->mustNowRecalculateFromScratch();

    // get all of the molecules in the system
    auto molnums = this->molNums();
    auto mols = this->molecules();

    // get the space and time values
    const auto space_property = map["space"];
    const auto time_property = map["time"];

    QString time_property_source("time");
    QString space_property_source("space");

    if (time_property.hasSource())
        time_property_source = QString(time_property.source());

    if (space_property.hasSource())
        space_property_source = QString(space_property.source());

    // we must get the space property and a time property
    SpacePtr space;
    SireUnits::Dimension::Time time;
    bool found_space = false;
    bool found_time = false;

    try
    {
        space = this->property(space_property_source).asA<Space>();
        found_space = true;
    }
    catch (...)
    {
    }

    try
    {
        time = get_time_from_property(this->property(time_property_source));
        found_time = true;
    }
    catch (...)
    {
    }

    if (not found_space)
    {
        if (space_property.hasValue())
        {
            try
            {
                space = space_property.value().asA<Space>();
                found_space = true;
            }
            catch (...)
            {
            }
        }

        if (not found_space)
        {
            for (const auto &mol : mols)
            {
                try
                {
                    space = mol.data().property(space_property_source).asA<Space>();
                    found_space = true;
                    break;
                }
                catch (...)
                {
                }
            }
        }
    }

    if (not found_time)
    {
        if (time_property.hasValue())
        {
            try
            {
                time = get_time_from_property(time_property.value());
                found_time = true;
            }
            catch (...)
            {
            }
        }

        if (not found_time)
        {
            for (const auto &mol : mols)
            {
                try
                {
                    time = get_time_from_property(mol.data().property(time_property_source));
                    found_time = true;
                    break;
                }
                catch (...)
                {
                }
            }
        }
    }

    if (not found_space)
    {
        space = Cartesian();
    }

    if (not found_time)
    {
        time = SireUnits::Dimension::Time(0);
    }

    // do we have an active SystemTrajectory?
    bool must_create = false;

    if (system_trajectory.constData() == 0)
    {
        must_create = true;
    }

    SystemTrajectory *traj = dynamic_cast<SystemTrajectory *>(system_trajectory.data());

    if (traj == 0)
    {
        must_create = true;
    }
    else
    {
        must_create = not traj->isCompatibleWith(molnums, mols, map);
    }

    if (must_create)
    {
        traj = new SystemTrajectory(molnums, mols, map);
        system_trajectory = traj;
    }

    // save the frame into the system_trajectory - this will automatically
    // update all molecules containing this trajectory
    traj->saveFrame(mols, space, time, Properties(), map);

    if (must_create)
    {
        // add this trajectory onto all of the molecules...
        auto mols2 = mols;

        for (auto it = mols.constBegin(); it != mols.constEnd(); ++it)
        {
            auto mol = it->molecule().data();

            Trajectory moltraj;

            if (mol.hasProperty(traj_prop))
            {
                moltraj = mol.property(traj_prop).asA<Trajectory>();
                moltraj.append(traj->getTrajectory(mol.number()));
            }
            else
            {
                moltraj = Trajectory(traj->getTrajectory(mol.number()));
            }

            mol.setProperty(traj_prop.source(), moltraj);
            mols2.update(mol);
        }

        this->update(mols2);
    }
}

void System::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    this->accept();
    this->mustNowRecalculateFromScratch();
    MolGroupsBase::deleteFrame(frame, map);
}

void System::deleteAllFrames(const SireBase::PropertyMap &map)
{
    this->accept();
    this->mustNowRecalculateFromScratch();
    this->system_trajectory = 0;
    MolGroupsBase::deleteAllFrames(map);
}

void System::makeWhole(const PropertyMap &map)
{
    if (this->needsAccepting())
    {
        this->accept();
    }

    if (not this->containsProperty(map["space"]))
        return;

    if (not this->property(map["space"]).isA<Space>())
        return;

    const auto &space = this->property(map["space"]).asA<Space>();

    if (not space.isPeriodic())
        return;

    PropertyMap m = map;
    m.set("space", space);

    // get a list of all molecules in the system
    const SelectorMol mols(*this);

    // find the center of the box
    AABox box;

    for (auto mol : mols)
    {
        box += mol.evaluate().aaBox(m);
    }

    const auto center = box.center();

    SelectorMol changed_mols;

    for (const auto &mol : mols)
    {
        auto new_mol = mol.move().makeWhole(center, m).commit();

        if (new_mol.data().version() != mol.data().version())
        {
            changed_mols.append(new_mol);
        }
    }

    if (not changed_mols.isEmpty())
    {
        Delta delta(*this, true);

        // this ensures that only a single copy of System is used - prevents
        // unnecessary copying
        this->operator=(System());
        delta.update(changed_mols.toMolecules());
        this->operator=(delta.apply());

        if (this->needsAccepting())
        {
            delta = Delta();
            this->accept();
        }
    }
}

void System::makeWhole()
{
    this->makeWhole(PropertyMap());
}

void System::makeWhole(const Vector &center, const PropertyMap &map)
{
    if (this->needsAccepting())
    {
        this->accept();
    }

    if (not this->containsProperty(map["space"]))
        return;

    if (not this->property(map["space"]).isA<Space>())
        return;

    const auto &space = this->property(map["space"]).asA<Space>();

    if (not space.isPeriodic())
        return;

    PropertyMap m = map;
    m.set("space", space);

    // get a list of all molecules in the system
    const SelectorMol mols(*this);

    SelectorMol changed_mols;

    for (const auto &mol : mols)
    {
        auto new_mol = mol.move().makeWhole(center, m).commit();

        if (new_mol.data().version() != mol.data().version())
        {
            changed_mols.append(new_mol);
        }
    }

    if (not changed_mols.isEmpty())
    {
        Delta delta(*this, true);

        // this ensures that only a single copy of System is used - prevents
        // unnecessary copying
        this->operator=(System());
        delta.update(changed_mols.toMolecules());
        this->operator=(delta.apply());

        if (this->needsAccepting())
        {
            delta = Delta();
            this->accept();
        }
    }
}

void System::makeWhole(const Vector &center)
{
    this->makeWhole(center, PropertyMap());
}

const char *System::typeName()
{
    return QMetaType::typeName(qMetaTypeId<System>());
}
