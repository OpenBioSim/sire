/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007   Christopher Woods
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

#include <Python.h>
#include <boost/python.hpp>

#include <QVector>
#include <QSet>

#include <boost/tuple/tuple.hpp>

#include "Helpers/convertlist.hpp"
#include "Helpers/convertdict.hpp"
#include "Helpers/convertset.hpp"
#include "Helpers/tuples.hpp"
#include "Base/convertpackedarray.hpp"

#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"

#include "SireMol/core.h"

#include "SireMM/bond.h"
#include "SireMM/angle.h"
#include "SireMM/dihedral.h"
#include "SireMM/improper.h"

#include "SireMM/ljparameter.h"
#include "SireMM/twoatomfunctions.h"
#include "SireMM/threeatomfunctions.h"
#include "SireMM/fouratomfunctions.h"
#include "SireMM/cljatoms.h"
#include "SireMM/amberparams.h"
#include "SireMM/gromacsparams.h"

#include "SireMM/cljboxes.h"

#include "SireBase/packedarray2d.hpp"

using namespace SireMM;

using boost::python::register_tuple;

void register_SireMM_containers()
{
  register_PackedArray<SireBase::PackedArray2D<LJParameter>>();

  register_tuple<boost::tuple<int, int, LJ1264Parameter>>();

  register_list<QList<boost::tuple<int, int, LJ1264Parameter>>>();

  register_list<QVector<LJParameter>>();
  register_list<QList<LJParameter>>();

  register_list<QVector<TwoAtomFunction>>();
  register_list<QVector<ThreeAtomFunction>>();
  register_list<QVector<FourAtomFunction>>();

  register_list<QVector<CLJAtom>>();

  register_list<QVector<CLJBox>>();
  register_list<QVector<CLJBoxIndex>>();

  register_list<QList<GromacsBond>>();
  register_list<QList<GromacsAngle>>();
  register_list<QList<GromacsDihedral>>();

  register_list<QVector<AmberDihPart>>();

  register_list<QList<Bond>>();
  register_list<QList<Angle>>();
  register_list<QList<Dihedral>>();
  register_list<QList<Improper>>();

  register_dict<QHash<BondID, AmberBond>>();
  register_dict<QHash<AngleID, AmberAngle>>();
  register_dict<QHash<DihedralID, AmberDihedral>>();
  register_dict<QHash<ImproperID, AmberDihedral>>();
  register_dict<QHash<BondID, AmberNB14>>();

  register_dict<QHash<QString, GromacsAtomType>>();
  register_dict<QMultiHash<QString, GromacsAtomType>>();
  register_dict<QMultiHash<QString, GromacsBond>>();
  register_dict<QMultiHash<QString, GromacsAngle>>();
  register_dict<QMultiHash<QString, GromacsDihedral>>();

  register_dict<QMultiHash<SireMol::BondID, GromacsBond>>();
  register_dict<QMultiHash<SireMol::AngleID, GromacsAngle>>();
  register_dict<QMultiHash<SireMol::DihedralID, GromacsDihedral>>();
}
