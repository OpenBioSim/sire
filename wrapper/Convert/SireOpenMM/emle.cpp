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

#include "emle.h"

using namespace SireOpenMM;

EMLECallback::EMLECallback()
{
}

EMLECallback::EMLECallback(bp::object py_object, QString callback) :
    py_object(py_object), callback(callback)
{
}

boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>
EMLECallback::call(
    QVector<int> numbers_qm,
    QVector<double> charges_mm,
    QVector<QVector<double>> xyz_qm,
    QVector<QVector<double>> xyz_mm)
{
    return bp::call_method<boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>>(
            this->py_object.ptr(),
            this->callback.toStdString().c_str(),
            numbers_qm,
            charges_mm,
            xyz_qm,
            xyz_mm
    );
}

EMLEEngine::EMLEEngine()
{
}

EMLEEngine::EMLEEngine(bp::object py_object, SireUnits::Dimension::Length cutoff, double lambda) :
    callback(py_object, "_sire_callback"),
    cutoff(cutoff),
    lambda(lambda)
{
}

void EMLEEngine::setLambda(double lambda)
{
    this->lambda = lambda;
}

double EMLEEngine::getLambda() const
{
    return this->lambda;
}

void EMLEEngine::setCutoff(SireUnits::Dimension::Length cutoff)
{
    this->cutoff = cutoff;
}

SireUnits::Dimension::Length EMLEEngine::getCutoff() const
{
    return this->cutoff;
}

QVector<int> EMLEEngine::getAtoms() const
{
    return this->atoms;
}

void EMLEEngine::setAtoms(QVector<int> atoms)
{
    this->atoms = atoms;
}

const char *EMLEEngine::typeName()
{
    return QMetaType::typeName(qMetaTypeId<EMLEEngine>());
}

const char *EMLEEngine::what() const
{
    return EMLEEngine::typeName();
}

OpenMM::ForceImpl *EMLEEngine::createImpl() const
{
#ifdef SIRE_USE_CUSTOMCPPFORCE
    return new EMLEEngineImpl(*this);
#else
    throw SireError::unsupported(QObject::tr(
                                     "Unable to create a EMLEEngine because OpenMM::CustomCPPForceImpl "
                                     "is not available. You need to use OpenMM 8.1 or later."),
                                 CODELOC);
    return 0;
#endif
}

EMLEEngineImpl::EMLEEngineImpl(const EMLEEngine &owner) :
    OpenMM::CustomCPPForceImpl(owner),
    owner(owner)
{
}

EMLEEngineImpl::~EMLEEngineImpl()
{
}

const EMLEEngine &EMLEEngineImpl::getOwner() const
{
    return this->owner;
}

double EMLEEngineImpl::computeForce(
    OpenMM::ContextImpl &context,
    const std::vector<OpenMM::Vec3> &positions,
    std::vector<OpenMM::Vec3> &forces)
{
    return 0;
}

boost::tuple<double, QVector<QVector<double>>, QVector<QVector<double>>>
EMLEEngine::call(
    QVector<int> numbers_qm,
    QVector<double> charges_mm,
    QVector<QVector<double>> xyz_qm,
    QVector<QVector<double>> xyz_mm)
{
    return this->callback.call(numbers_qm, charges_mm, xyz_qm, xyz_mm);
}