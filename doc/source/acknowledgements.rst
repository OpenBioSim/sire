================
Acknowledgements
================

We give huge thanks to everyone who has contributed to :mod:`sire` development,
whether directly in the form of writing code, tests or documentation,
or indirectly via support, discussions or submitting issues or
bug reports.

We gratefully acknowledge funding from the
`EPSRC <https://epsrc.ukri.org>`__ and
`BBSRC <https://bbsrc.ukri.org>`__ who contributed
funding to the initial development of :mod:`sire`. We are also thankful
to `UCB <https://www.ucb.com>`__, `Cresset <https://www.cresset-group.com>`__,
`Exscientia <https://www.exscientia.ai>`__ and
`Evotec <https://www.evotec.com/en>`__ who have all either
directly funded development, or have funded researchers
who have contributed to :mod:`sire`.

We also thank the `Software Sustainability Institute <https://software.ac.uk>`__
for many useful discussions.

We thank the Universities of
`Bristol <https://bristol.ac.uk>`__ and
`Edinburgh <https://ed.ac.uk>`__ for providing the
time to the members of staff who have contributed to :mod:`sire`'s development.

We thank `CCP-BioSim <https://ccpbiosim.ac.uk>`__ who have also provided
guidance and encouragement during the development of this software.

Website
=======

This website was generated using `sphinx <https://www.sphinx-doc.org/en/master/index.html>`__,
using a modified version of the `furo theme <https://pradyunsg.me/furo/>`__.

Hosting
=======

:mod:`sire` is developed on `GitHub <https://github.com>`__, making extensive
use of its many excellent features. This include using
GitHub pages for hosting this website, and
`GitHub actions <https://github.com/OpenBioSim/sire/actions>`__ for
CI/CD.

:mod:`sire` binary packages are hosted on `conda-forge <https://conda-forge.org>`__.

:mod:`sire` containers are hosted on `docker hub <https://hub.docker.com>`__.

The :mod:`sire` `notebook service <https://try.openbiosim.org>`__ is hosted
in a `JupyterHub <https://jupyterhub.readthedocs.io/en/stable/>`__ cluster,
built following the instructions on
`Zero to JupyterHub <https://jupyterhub.readthedocs.io/en/stable/>`__.
This is hosted in a `kubernetes <https://kubernetes.io>`__ cluster
on `Microsoft Azure <https://azure.microsoft.com/en-gb/>`__.

Third Party Software
====================

:mod:`sire` depends on a lot of third party software, the details and licenses of
which can be found below. The software will be installed automatically
as part of the installing the :mod:`sire` conda package, so you shouldn't
have to do anything yourself.

:mod:`sire` is itself distributed under the terms of the GPL version 3
(or any later GPL license). The C++ source code is licensed
under the GPL 3 or later, and linking with GPL3 dependencies
(e.g. GSL) means that the entire package is licensed under GPL 3
or later.

C++ Dependencies
================

Qt 5
----

:mod:`sire` is built on top of Qt.Core from Qt 5. This is used under the terms
of the `LGPL 2 <http://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html>`__
or later. :mod:`sire` does not modify Qt, so this license allows both
commercial and non-commercial use without fee. You can find out more
about Qt and its license `from here <https://www.qt.io/terms-conditions/>`__.

Python
------

The :mod:`sire` C++ library is wrapped up and made available for use within Python 3.
This is used under the terms of the `PSF license <https://docs.python.org/3/license.html>`__,
which is compatible with
the GPLv3. The license allows both commercial and non-commercial use
without fee. You can find out more about Python and its license
`from here <https://www.python.org/>`__.

boost
-----

:mod:`sire` uses many of the components from the boost libraries, in particular
the boost::python module that is used to wrap up the C++ code.
This is used under the terms of the
`Boost Software License <http://www.boost.org/users/license.html>`__,
which allows both commercial and non-commercial
use without fee. This is compatible with the GPLc3. You can find out
more about boost and its license `from here <http://www.boost.org/>`__.

Py++
----

:mod:`sire` uses Py++ to auto-generate all of the C++ python wrappers. Py++ uses
either GCCXML or CastXML, and, as it is used as a tool, its license does
not affect :mod:`sire`. You can read more about Py++
`from here <http://pyplusplus.readthedocs.io/en/latest/>`__.

cmake
-----

:mod:`sire` uses cmake as its build system. As it is used as a tool, its license
does not affect :mod:`sire`. CMake is excellent. You can read more about it
`from here <https://cmake.org/>`__.

Anaconda
---------

:mod:`sire` uses Anaconda Python (specifically miniforge and conda-forge) to
simplify the management and installation of Python and the various
modules on which :mod:`sire` depends.

Anaconda (and miniforge) are distributed as
`open source projects <https://www.continuum.io/open-source-core-modern-software>`__.
As :mod:`sire` does not explicitly link with them, the license is not an issue.
You can find out more about Anaconda `from here <https://www.continuum.io/>`__.

Threading Building Blocks (tbb)
-------------------------------

:mod:`sire` uses the `Threading Building Blocks <https://www.threadingbuildingblocks.org/>`__
library for within-node
parallelisation. This is licensed under the open source
`Apache 2.0 license <https://www.threadingbuildingblocks.org/faq/10>`__.

Gnu Scientific Library (GSL)
----------------------------

:mod:`sire` uses some of the routines from the Gnu Scientific Library.
This is used under the terms of the `GPL v3 <http://www.gnu.org/copyleft/gpl.html>`__
license. More information
about GSL and its license can be `found here <http://www.gnu.org/software/gsl/>`__.

NetCDF
------

:mod:`sire` links to the `NetCDF <https://docs.unidata.ucar.edu/netcdf-c/current/copyright.html>`__
library so that it can read/write Amber binary files. NetCDF
is openly licensed under a BSD-style license, and is compatible
with the GPL.

OpenMM
------

:mod:`sire` links to `OpenMM <https://openmm.org>`__ to perform accelerated
dynamics (e.g. as part of the ``somd`` program), and also via an
interconversion layer to switch between sire and OpenMM system representations.
This is licensed under either the MIT or LGPL licenses,
so compatible with the GPL.

Regress
-------

:mod:`sire` uses the linear least squares regression library, `regress`, for
polynomial least squares fitting. This is used under the terms of
the GPLv3 license.

The source code for this module can be
`found here <https://github.com/openbiosim/sire/blob/devel/corelib/src/libs/:mod:`sire`Analysis/third_party/regress.cpp>`__.

eig3
----

:mod:`sire` uses the eig3 library for eigenvector/eigenmatrix calculations by
Connelly Barnes. This is in the public domain, and is derived itself
from the Java matrix library JAMA (also public domain).

Information about this can be
`found here <http://barnesc.blogspot.co.uk/2007/02/eigenvectors-of-3x3-symmetric-matrix.html>`__,
with the license within :mod:`sire` `found here <https://github.com/openbiosim/sire/blob/devel/corelib/src/libs/:mod:`sire`Maths/third_party/eig3/readme.txt>`__.

Mersenne Twister
----------------

:mod:`sire` uses the Mersenne Twister program by Richard Wagner for the generation
of random numbers. This is used under a BSD-style license, shown below.

::

 // Mersenne Twister random number generator -- a C++ class MTRand
 // Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
 // Richard J. Wagner  v1.0  15 May 2003  rjwagner@writeme.com

 // The Mersenne Twister is an algorithm for generating random numbers.  It
 // was designed with consideration of the flaws in various other generators.
 // The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
 // are far greater.  The generator is also fast; it avoids multiplication and
 // division, and it benefits from caches and pipelines.  For more information
 // see the inventors' web page at http://www.math.keio.ac.jp/~matumoto/emt.html

 // Reference
 // M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
 // Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
 // Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

 // Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 // Copyright (C) 2000 - 2003, Richard J. Wagner
 // All rights reserved.
 //
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions
 // are met:
 //
 //   1. Redistributions of source code must retain the above copyright
 //      notice, this list of conditions and the following disclaimer.
 //
 //   2. Redistributions in binary form must reproduce the above copyright
 //      notice, this list of conditions and the following disclaimer in the
 //      documentation and/or other materials provided with the distribution.
 //
 //   3. The names of its contributors may not be used to endorse or promote
 //      products derived from this software without specific prior written
 //      permission.
 //
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 // CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 // EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 // PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 // PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 // LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 // NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 // SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 // The original code included the following notice:
 //
 //     When you use this, send an email to: matumoto@math.keio.ac.jp
 //     with an appropriate reference to your work.
 //
 // It would be nice to CC: rjwagner@writeme.com and Cokus@math.washington.edu
 // when you write.

I must remember to send them an email…

More information about Mersenne Twister can be
`found here <http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/c-lang.html>`__.

sse_mathfun
-----------

:mod:`sire` uses sse_mathfun for vectorising intrinsic maths
functions on processors that support SSE. These libraries were written
by Julien Pommier, and released under the BSD-style zlib license,
which is given here.

::

 /* Copyright (C) 2007  Julien Pommier
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license)
 */


avx_mathfun
-----------

This is an AVX library inspired by sse_mathfun, that extends support to
processors with AVX instructions. It was written by Giovanni Garberoglio,
and is also under a BSD-style zlib license.

::

  AVX implementation of sin, cos, sincos, exp and log
   Based on "sse_mathfun.h", by Julien Pommier
   http://gruntthepeon.free.fr/ssemath/
   Copyright (C) 2012 Giovanni Garberoglio
   Interdisciplinary Laboratory for Computational Science (LISC)
   Fondazione Bruno Kessler and University of Trento
   via Sommarive, 18
   I-38123 Trento (Italy)
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license)

sse2neon
--------

:mod:`sire` uses `sse2neon <https://github.com/DLTcollab/sse2neon>`__ to
convert the hand-vectorised SSE2 code to ARM NEON instructions. This
worked amazingly well, and we sincerely thank
this project for producing such a well-documented and performant
header-only library. This is used under a MIT license.

LAP (Linear Assignment Problem Solver)
--------------------------------------

:mod:`sire` implements its own C++ version of the LAP library for solving the
linear assignment problem. This is `available here <https://github.com/openbiosim/sire/blob/devel/corelib/src/libs/SireMaths/linearap.cpp>`__.

The original code is Freeware, with more information about it available
`from here <http://www.assignmentproblems.com/linearAP.htm>`__.

MD5
---

:mod:`sire` uses the MD5 library written by L. Peter Deutsch.
It is used under a BSD-style license, given below.

::

  Copyright (C) 1999, 2002 Aladdin Enterprises.  All rights reserved.
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  L. Peter Deutsch
  ghost@aladdin.com

More information about MD5 libraries in general can be
`found here <http://userpages.umbc.edu/~mabzug1/cs/md5/md5.html>`__.

kabasch fitting
---------------

I have written a C++ implementation of the kabasch algorithm for alignment.
This was inspired by the calculate_rmsd python script written by
Jimmy Charnley Kromann and Lars Bratholm,
available https://github.com/charnley/rmsd, and under license;

::

        =====================
        Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

        1. Redistributions of source code must retain the above copyright notice, this
           list of conditions and the following disclaimer.
        2. Redistributions in binary form must reproduce the above copyright notice,
           this list of conditions and the following disclaimer in the documentation
           and/or other materials provided with the distribution.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        =======================

mdanalysis
----------

I have written a C++ implementation of the ``_infer_bo_and_charges``
algorithm from `mdanalysis <https://mdanalysis.org>`__.

The algorithm is `described here <https://blog.matteoferla.com/2020/02/guess-bond-order-in-rdkit-by-number-of.html>`__
and `here <https://docs.mdanalysis.org/2.0.0/_modules/MDAnalysis/converters/RDKit.html#_infer_bo_and_charges>`__.

The code was written by Matteo Ferla and was released under the GPL.

gemmi
-----

:mod:`sire` uses the `gemmi <https://gemmi.readthedocs.io/en/latest/>`__ library
for reading and writing PDBx/mmCIF files. This is used under the terms of
the LGPL3 license. Thanks to this project for providing such a useful
header-only C++ library that is easy to bring incorporate into other
programs, and that is so well documented and feature-rich. This is provided
via an integration that supports interconversion between sire and gemmi objects.

PEGTL
-----

`gemmi <https://gemmi.readthedocs.io/en/latest/>`__ uses the
`PEGTL <https://github.com/taocpp/PEGTL>`__ library
for parsing mmCIF files. This is used under the terms of the MIT
or Boost Software License (the version bundled with gemmi is MIT,
while the version on GitHub is BSL)

pyboost11
---------

Thanks to Yung-Yu's excellent
`blog post <https://yyc.solvcon.net/en/latest/writing/2021/pyboost11/pyboost11.html>`__
on how to interconvert between boost::python and pybind11, and for providing
`pyboost11.hpp <https://yyc.solvcon.net/en/latest/writing/2021/pyboost11/code.html#pyboost11-code-pyboost11-hpp>`__
as a header-only library to support interconversion. This was really useful
to let me mix pybind11-wrapped gemmi objects with boost::python-wrapped sire.
The code is used under the terms of the liberal open source license
given at the top of the header file.

RDKit
-----

:mod:`sire` provides an integration with the `RDKit <https://www.rdkit.org>`__
package, supporting interconversion between the sire and RDKit
molecular representations. This is used under the terms of the BSD3 license.

Python Dependencies
===================

ap (ascii plot)
---------------

:mod:`sire` bundles the Python “ap” library for drawing ascii graphs.
This is available as “Sire.Tools.ap”

Version 0.9 written by M. Fouesneau is included, available
freely from `GitHub here <https://github.com/mfouesneau/asciiplot>`__.
The only change I've made is running this through Python's 2to3
program to make this code work with Python 3.

The header documentation reads;

::

 Package that allows you to plot simple graphs in ASCII, a la matplotlib.
 This package is a inspired from Imri Goldberg's ASCII-Plotter 1.0
 (https://pypi.python.org/pypi/ASCII-Plotter/1.0)
 At a time I was enoyed by security not giving me direct access to my computer,
 and thus to quickly make figures from python, I looked at how I could make
 quick and dirty ASCII figures. But if I were to develop something, I wanted
 something that can be used with just python and possible standard-ish packages
 (numpy, scipy).
 So I came up with this package after many iterations based of ASCII-plotter.
 I added the feature to show multiple curves on one plot with different markers.
 And I also made the usage, close to matplotlib, such that there is a plot,
 hist, hist2d and imshow functions.

 TODO:
     imshow does not plot axis yet.
     make a correct documentation

lazy_import
-----------

:mod:`sire` uses `lazy_import <https://github.com/mnmelo/lazy_import>`__ to
lazy load the modules. This is licensed under the GPLv3.

rich
----

:mod:`sire` uses `rich <https://github.com/Textualize/rich>`__ to provide
rich console output when printing. This is licensed under the
GPL-compatible MIT license.

pandas
------

:mod:`sire` uses `pandas <https://pandas.pydata.org/docs/>`__ to output
data in DataFrames that can be more easily operated on and explored
by users. Pandas is BSD-licensed, so compatible with the GPL.
