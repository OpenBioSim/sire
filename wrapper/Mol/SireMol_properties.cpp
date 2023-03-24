#include <Python.h>
#include <boost/python.hpp>

#include "Base/convertproperty.hpp"
#include "SireMol_properties.h"

#include "SireError/errors.h"
#include "SireStream/datastream.h"
#include "SireUnits/units.h"
#include "atomidentifier.h"
#include "atomidx.h"
#include "atommatcher.h"
#include "atommatchers.h"
#include "atomname.h"
#include "atomselection.h"
#include "evaluator.h"
#include "moleculeinfodata.h"
#include "moleculeview.h"
#include "mover.hpp"
#include "tostring.h"
#include "atommatcher.h"
#include "SireBase/errors.h"
#include "SireError/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "atombeads.h"
#include "atomidx.h"
#include "atomselection.h"
#include "beadidx.h"
#include "beading.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"
#include <boost/noncopyable.hpp>
#include "beading.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireVol/coordgroup.h"
#include "atom.h"
#include "atomcoords.h"
#include "atomelements.h"
#include "atomselection.h"
#include "bondhunter.h"
#include "connectivity.h"
#include "molecule.h"
#include "moleculedata.h"
#include "moleculeinfodata.h"
#include "moleculeview.h"
#include "mover.hpp"
#include "selector.hpp"
#include <QDebug>
#include <QMutex>
#include "bondhunter.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "cuttingfunction.h"
#include "molecule.h"
#include "moleditor.h"
#include "mover.hpp"
#include "residuecutting.h"
#include "selector.hpp"
#include <QMutex>
#include "cuttingfunction.h"
#include "SireBase/incremint.h"
#include "SireBase/majorminorversion.h"
#include "SireBase/refcountdata.h"
#include "SireError/errors.h"
#include "SireID/index.h"
#include "SireMol/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "core.h"
#include "editor.hpp"
#include "mgname.h"
#include "mgnum.h"
#include "molecule.h"
#include "moleculegroup.h"
#include "molidentifier.h"
#include "molidx.h"
#include "molname.h"
#include "molnum.h"
#include "mover.hpp"
#include "partialmolecule.h"
#include "select.h"
#include "tostring.h"
#include <QDebug>
#include <QMutex>
#include <QVector>
#include <boost/tuple/tuple.hpp>
#include "moleculegroup.h"
#include "SireBase/slice.h"
#include "SireError/errors.h"
#include "SireMol/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "atom.h"
#include "atomid.h"
#include "cgid.h"
#include "chain.h"
#include "chainid.h"
#include "cutgroup.h"
#include "editor.hpp"
#include "mgidx.h"
#include "mgname.h"
#include "mgnum.h"
#include "molecule.h"
#include "moleculegroups.h"
#include "molecules.h"
#include "molidx.h"
#include "molname.h"
#include "molnum.h"
#include "mover.hpp"
#include "partialmolecule.h"
#include "resid.h"
#include "residue.h"
#include "segid.h"
#include "segment.h"
#include "select.h"
#include "selector.hpp"
#include "tostring.h"
#include "viewsofmol.h"
#include <QDebug>
#include <QMutex>
#include "moleculegroups.h"
#include "SireBase/errors.h"
#include "SireBase/generalunitproperty.h"
#include "SireBase/slice.h"
#include "SireError/errors.h"
#include "SireMol/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "SireVol/space.h"
#include "atom.h"
#include "atomselection.h"
#include "chain.h"
#include "cutgroup.h"
#include "molecule.h"
#include "molecules.h"
#include "moleculeview.h"
#include "residue.h"
#include "segment.h"
#include "select.h"
#include "selector.hpp"
#include "selectormol.h"
#include "trajectory.h"
#include <QDebug>
#include "moleculeview.h"
#include "SireCAS/identities.h"
#include "SireCAS/values.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "core.h"
#include "geometryperturbation.h"
#include "molecule.h"
#include "moleditor.h"
#include "mover.hpp"
#include "perturbation.h"
#include "perturbation.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"
#include "atomselection.h"
#include "core.h"
#include "editor.hpp"
#include "evaluator.h"
#include "moleculedata.h"
#include "moleculeview.h"
#include "mover.hpp"
#include "partialmolecule.h"
#include "weightfunction.h"
#include "weightfunction.h"
void register_SireMol_properties()
{
    register_property_container< SireMol::AtomMatcherPtr, SireMol::AtomMatcher >();
    register_property_container< SireMol::BeadingPtr, SireMol::Beading >();
    register_property_container< SireMol::BondHunterPtr, SireMol::BondHunter >();
    register_property_container< SireMol::CutFuncPtr, SireMol::CuttingFunction >();
    register_property_container< SireMol::MolGroupPtr, SireMol::MoleculeGroup >();
    register_property_container< SireMol::MolGroupsPtr, SireMol::MolGroupsBase >();
    register_property_container< SireMol::MolViewPtr, SireMol::MoleculeView >();
    register_property_container< SireMol::PerturbationPtr, SireMol::Perturbation >();
    register_property_container< SireMol::WeightFuncPtr, SireMol::WeightFunction >();
}
