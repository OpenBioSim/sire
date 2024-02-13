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

#ifndef SIREOPENMM_QMMM_H
#define SIREOPENMM_QMMM_H

#include "OpenMM.h"
#include "openmm/Force.h"

#include <boost/tuple/tuple.hpp>

#include <QMap>
#include <QVector>

#include "sireglobal.h"

#include "SireBase/property.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireOpenMM
{
    class QMForce;
    class QMEngine;
    class NullQMEngine;

    class QMForce : public OpenMM::Force
    {
    public:
        QMForce();

        virtual ~QMForce();

        //! Set the lambda weighting factor.
        virtual void setLambda(double lambda) = 0;

    protected:
        virtual OpenMM::ForceImpl *createImpl() const = 0;
    };

    class QMEngine : public SireBase::Property
    {
    public:
        QMEngine();

        virtual ~QMEngine();

        //! Clone the QM engine.
        virtual QMEngine *clone() const = 0;

        //! Get the name of the QM engine.
        static const char *typeName();

        //! Get the name of the QM engine.
        const char *what() const;

        //! Create a QM force object.
        virtual QMForce* createForce() const=0;

        //! Get a null QM engine.
        static const NullQMEngine &null();
    };

    typedef SireBase::PropPtr<SireOpenMM::QMEngine> QMEnginePtr;

    class NullQMEngine : public SireBase::ConcreteProperty<NullQMEngine, QMEngine>
    {
    public:
        NullQMEngine();

        ~NullQMEngine();

        //! Get the name of the QM engine.
        static const char *typeName();

        //! Get the name of the QM engine.
        const char *what() const;

        //! Create a QM force object.
        QMForce* createForce() const;
    };
}

Q_DECLARE_METATYPE(SireOpenMM::NullQMEngine)

SIRE_EXPOSE_CLASS(SireOpenMM::QMForce)
SIRE_EXPOSE_CLASS(SireOpenMM::QMEngine)
SIRE_EXPOSE_CLASS(SireOpenMM::NullQMEngine)

SIRE_EXPOSE_PROPERTY(SireOpenMM::QMEnginePtr, SireOpenMM::QMEngine)

SIRE_END_HEADER

#endif
