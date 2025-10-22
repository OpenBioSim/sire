#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireIO_properties.h"

#include "SireBase/booleanproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/progressbar.h"
#include "SireBase/propertylist.h"
#include "SireBase/releasegil.h"
#include "SireBase/stringproperty.h"
#include "SireBase/timeproperty.h"
#include "SireError/errors.h"
#include "SireFF/ffdetail.h"
#include "SireIO/errors.h"
#include "SireMM/mmdetail.h"
#include "SireMol/core.h"
#include "SireMol/mgname.h"
#include "SireMol/mgnum.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/molidx.h"
#include "SireMol/trajectory.h"
#include "SireMol/trajectoryaligner.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireSystem/system.h"
#include "SireUnits/units.h"
#include "filetrajectory.h"
#include "filetrajectoryparser.h"
#include "moleculeparser.h"
#include "supplementary.h"
#include <QDebug>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QMutex>
#include <QTextStream>
#include "moleculeparser.h"
#include "SireError/errors.h"
#include "SireMol/cuttingfunction.h"
#include "SireMol/molecule.h"
#include "SireMol/molidx.h"
#include "SireMol/mover.hpp"
#include "SireStream/datastream.h"
#include "iobase.h"
#include <QDebug>
#include <QFile>
#include "iobase.h"
void register_SireIO_properties()
{
    register_property_container< SireIO::MoleculeParserPtr, SireIO::MoleculeParser >();
    register_property_container< SireIO::IOPtr, SireIO::IOBase >();
}
