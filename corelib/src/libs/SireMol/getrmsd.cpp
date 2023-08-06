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

#include "getrmsd.h"

#include "SireBase/releasegil.h"
#include "SireBase/progressbar.h"
#include "SireBase/parallel.h"

#include "SireMaths/align.h"

#include "SireUnits/units.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireMaths;

namespace SireMol
{

    QVector<Length> get_rmsd(
        const SelectorM<Atom> &atoms,
        const SelectorM<Atom> &reference,
        const TrajectoryAligner &aligner,
        const QList<qint64> &frames,
        const PropertyMap &map)
    {
        QVector<Length> rmsds;

        if (frames.isEmpty())
            return rmsds;

        // release the GIL here so that progress bars can be displayed
        auto handle = SireBase::release_gil();

        rmsds = QVector<Length>(frames.count());
        auto rmsds_data = rmsds.data();

        const auto coords_property = map["coordinates"];

        SireBase::ProgressBar bar(frames.count(), "Calculate RMSD");

        const auto coords = reference.property<Vector>(coords_property).toVector();

        const bool needs_aligning = aligner.count() != 0;

        bar = bar.enter();

        PropertyMap my_map = map;
        my_map.set("coords_only", BooleanProperty(true));

        if (false) // should_run_in_parallel(frames.count(), map))
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, frames.count()),
                              [&](const tbb::blocked_range<int> &r)
                              {
                                  SelectorM<Atom> frame = atoms;
                                  PropertyMap frame_map = my_map;

                                  for (int i = r.begin(); i < r.end(); ++i)
                                  {
                                      if (needs_aligning)
                                      {
                                          auto transform = aligner[frames[i]];

                                          if (not transform.isNull())
                                              frame_map.set("transform", transform);
                                      }

                                      frame.loadFrame(frames[i], frame_map);

                                      const auto frame_coords = frame.property<Vector>(coords_property).toVector();

                                      rmsds_data[i] = SireMaths::getRMSD(coords, frame_coords) * angstrom;

                                      bar.tick();
                                  }
                              });
        }
        else
        {
            SelectorM<Atom> frame = atoms;
            PropertyMap frame_map = my_map;

            for (int i = 0; i < frames.count(); ++i)
            {
                if (needs_aligning)
                {
                    auto transform = aligner[frames[i]];

                    if (not transform.isNull())
                        frame_map.set("transform", transform);
                }

                frame.loadFrame(frames[i], frame_map);

                const auto frame_coords = frame.property<Vector>(coords_property).toVector();

                rmsds_data[i] = SireMaths::getRMSD(coords, frame_coords) * angstrom;

                bar.tick();
            }
        }

        bar.success();

        return rmsds;
    }
}
