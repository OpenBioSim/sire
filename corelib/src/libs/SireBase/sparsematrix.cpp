/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2024  Christopher Woods
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

#include "sparsematrix.hpp"

#include "sireglobal.h"

using namespace SireBase;
using namespace SireStream;

static const auto sparse_magic = Sire::getMagic("SireBase::SparseMatrix");

namespace SireBase
{
    namespace detail
    {
        void writeSparseMatrixMagic(QDataStream &ds, VersionID v)
        {
            ds << sparse_magic << v;
        }

        VersionID checkSparseMatrixMagic(QDataStream &ds)
        {
            Sire::MagicID id;
            VersionID v;

            // do this in a transaction, so that any exceptions thrown here
            // will not affect the stream (i.e. we can try to read this again)
            ds.startTransaction();

            ds >> id >> v;

            if (id != sparse_magic)
            {
                ds.rollbackTransaction();
                // we have to assume this is an older SparseMatrix
                return 1;
            }

            ds.commitTransaction();

            return v;
        }

        void throwSparseMatrixMagicError(VersionID v, const QString &supported)
        {
            throw SireStream::version_error(
                QObject::tr("Invalid version for SparseMatrix - got %1 but only version(s) %2 are supported")
                    .arg(v)
                    .arg(supported),
                CODELOC);
        }
    }
}