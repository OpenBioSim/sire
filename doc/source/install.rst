============
Installation
============

Binary packages for :mod:`sire` are available on MacOS, Linux and Windows
running on Intel/AMD (X86-64) or ARM64 processors.

:mod:`sire` can be compiled on any UNIX or Windows-compatible operating system
running on X86-64, ARM64 or PowerPC processors.

You have a range of options for installing the software.

1. No-installation - Run in a Web Browser
=========================================

We run a completely free `JupyterHub <https://try.openbiosim.org>`__ on
which we have sire installed.

This is at `try.openbiosim.org <https://try.openbiosim.org>`__.
You only need a `GitHub account <https://github.com>`__, which is
used to log into the server.

Simply go to `try.openbiosim.org <https://try.openbiosim.org>`__ in your
web browser, and log in using your `GitHub account <https://github.com>`__.
This will start a Jupyter Lab instance. In here you can start a terminal,
and then use :mod:`sire` directly via ``ipython``. Or you can start a Jupyter
notebook and use :mod:`sire` there.

To import :mod:`sire`, at a Python prompt type

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use :mod:`sire` or the
:doc:`quickstart guide <quickstart/index>` if you want an overview.

.. note::

   This free JupyterHub server is limited. You only have up to 2 GB of
   memory and at most 1 processor core. Disk storage is temporary,
   and any data will be lost when you log out. Because it only
   supports a limited number of concurrent users, inactive sessions will be
   automatically stopped and logged out after 20 minutes. Please only
   use this service to explore and learn :mod:`sire`.
   Do not use it for production work.

2. Easy installation - Run in a conda environment
=================================================

The easiest way to install :mod:`sire` is in a new
`conda environment <https://anaconda.org>`__.

You can use any conda environment or installation. We recommend using
`miniforge3 <https://github.com/conda-forge/miniforge#miniforge3>`__,
as this is pre-configured to use `conda-forge <https://conda-forge.org>`__.

.. _Install_miniforge:
Either... Install a new copy of ``miniforge``
----------------------------------------------

To install a new copy of
`miniforge <https://github.com/conda-forge/miniforge#miniforge3>`__,
first download a ``miniforge3`` from
`this page <https://github.com/conda-forge/miniforge#miniforge3>`__ that
matches your operating system and processor.

Install ``miniforge`` following the
`instructions here <https://github.com/conda-forge/miniforge#install>`__.

Once installed, you should be able to run the ``conda`` command to
install other packages (e.g. ``conda -h`` will print out help on
how to use the ``conda`` command).

Or... Use an existing anaconda/miniconda install
------------------------------------------------

If you want to use an existing anaconda or miniconda installation,
then first open a terminal with that distribution activated.
For example, open a terminal via anaconda navigator, or
open a terminal and run
``source /path/to/conda/bin/activate``, where ``/path/to/conda`` is
the full path to your anaconda or miniconda installation.

You should now be able to run the ``conda`` command to install other
packages (e.g. ``conda -h`` will print out help on how to use the
``conda`` command).

And then... Install sire into a new environment
-----------------------------------------------

We recommend that :mod:`sire` is installed into a new (clean) environment.
This minimises the risk of failures caused by incompatible dependencies.

Sire is currently packaged for Python 3.10, 3.11 and Python 3.12. We will start
by creating a Python 3.12 environment that we will call ``openbiosim``.

.. code-block:: bash

   $ conda create -n openbiosim "python<3.13"

.. note::

   We use ``python<3.13`` as this will install the most recent 3.12
   release of python.

We can now install :mod:`sire` into that environment by typing

.. code-block:: bash

   $ conda install -n openbiosim -c conda-forge -c openbiosim sire

.. note::

   The option ``-n openbiosim`` tells ``conda`` to install :mod:`sire`
   into the ``openbiosim`` environment. The option ``-c conda-forge``
   tells ``conda`` to use the ``conda-forge`` channel for all
   dependencies. The option ``-c openbiosim``
   tells ``conda`` to install :mod:`sire` from the ``openbiosim``
   conda channel.

If you want the latest development release, then install by typing

.. code-block:: bash

   $ conda install -n openbiosim -c conda-forge -c "openbiosim/label/dev" sire

You can install a specific version of sire by specifying the version number
in the conda install command, e.g.

.. code-block:: bash

    conda install -n openbiosim -c conda-forge -c openbiosim sire==2024.1.0

Note that limited space means that we can only keep a small number of
versions of sire on the official openbiosim conda channel. Generally
these are all point releases of the latest major version, plus the latest
point release of the last major version.

We do provide an
`archive channel <https://openbiosim.blob.core.windows.net/archive/index.html>`__
of all previous releases. You can search this archive channel for the
release you are interested in using the following command:

.. code-block:: bash

    conda search -c https://openbiosim.blob.core.windows.net/archive sire

This will return a list of all versions of sire available in the archive.

You can install a specific version from the archive using a command like:

.. code-block:: bash

    conda install -n openbiosim -c https://openbiosim.blob.core.windows.net/archive sire==2023.2.3

You may (optionally) want to install additional tools such as
``ipython`` and ``jupyterlab``. To do this, type

.. code-block:: bash

   $ conda install -n openbiosim ipython jupyterlab

To run :mod:`sire`, you must now activate the ``openbiosim`` environment.
You can do this by typing

.. code-block:: bash

   $ conda activate openbiosim

You can now start a Python session (e.g. running ``python``, or
``ipython`` or ``jupyter lab`` if you installed those). At the
Python prompt you can import :mod:`sire` by typing

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use :mod:`sire` or the
:doc:`quickstart guide <quickstart/index>` if you want an overview.

3. Also easy installation - Run in a container
==============================================

.. warning::

   Because of low demand, pre-built containers are created only for
   major releases, and may be out of date compared to the newest release.
   Please `get in touch <https://github.com/OpenBioSim/sire/issues>`__
   if you want to use a container and would like us to build the latest
   release.

Another route to install :mod:`sire` is to download and run our
pre-built containers. These can be run via
`docker <https://www.docker.com>`__ (on Linux, MacOS and Windows)
or via `podman <https://podman.io>`__ (on Linux) on Intel (X86-64)
or ARM64 processors.

To run via `docker <https://www.docker.com>`__, simply type;

.. code-block:: bash

   $ docker run -p 8888:8888 -it openbiosim/sire:latest

or, via `podman <https://podman.io>`__, type;

.. code-block:: bash

   $ podman run -p 8888:8888 -it openbiosim/sire:latest

This will download the container from
`hub.docker.com <https://anaconda.org/openbiosim/sire>`__ and
will start a command prompt in that container.

You can now type ``python``, ``ipython`` or ``jupyter lab``
to start a python, ipython or jupyter lab session.

.. note::

   The option ``-p 8888:8888`` tells docker/podman to redirect
   port ``8888`` on your computer to port ``8888`` in the
   container. This will let you open a browser and navigate to
   the URL printed by ``jupyter lab`` if you are using jupyter.
   You can drop this option if you don't want to use
   ``jupyter lab``.

.. note::

   You can map directories from your computer into the container
   by using the ``-v`` option. For example,
   ``-v $HOME/input:/home/openbiosim/input`` would map your
   ``input`` folder in your home directory to the ``input`` folder
   in the home directory of the container. This will let :mod:`sire`
   read and write files on your computer.

You can now start a Python session (e.g. running ``python``, or
``ipython`` or ``jupyter lab`` if you installed those). At the
Python prompt you can import :mod:`sire` by typing

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use :mod:`sire` or the
:doc:`quickstart guide <quickstart/index>` if you want an overview.

4. Harder installation - Compile from source
============================================

Sometimes you will want to compile and run :mod:`sire` from source.
This could be because we don't distribute a binary package for your
operating system, or because you want to use a newer version
(e.g. code from the ``devel`` branch, or from your own feature
branch if you are a developer).

Prerequisites
-------------

You need `pixi <https://pixi.sh>`__ installed to manage the build
environment and dependencies. Follow the
`pixi installation instructions <https://pixi.sh/latest/#installation>`__
for your platform.

.. note::

   You need to have Visual Studio C++ (2017 or newer) installed to compile on Windows.
   The easiest way to do this is to install the free
   `Visual Studio 2022 Community Edition <https://visualstudio.microsoft.com/vs/community>`__.
   Make sure to install "Desktop development with C++",
   including the options "MSVC v143 - VS 2022 C++ x64/x86 build tools (v14.30)",
   "C++ CMake tools for Windows", and at least one of "Windows 11 SDK" and/or
   "Windows 10 SDK" (any version will do). Currently only the X64 compilers
   have been tested - we are interested to try Windows/ARM64 once more of
   the dependencies are available.

Download the source code
------------------------

Download the latest development version of :mod:`sire` by typing;

.. code-block:: bash

   $ git clone https://github.com/openbiosim/sire

This will download into a directory called :mod:`sire`. Navigate into
this directory (e.g. ``cd sire``).

You can change to a different branch using the ``git checkout BRANCH``
command, e.g.

.. code-block:: bash

   $ git checkout main

will check out the ``main`` branch of :mod:`sire`. This always corresponds
to the last released version of :mod:`sire`. Or, you can check out a
feature branch using

.. code-block:: bash

   $ git checkout feat_name

where ``feat_name`` should be replaced by the name of the feature
branch you want to compile.

Create the build environment
-----------------------------

Use ``pixi`` to create and activate the development environment. This will
install all required dependencies, including compilers:

.. code-block:: bash

   $ pixi install -e dev
   $ pixi shell -e dev

Several environments are available depending on your needs:

* ``default`` - core sire dependencies only
* ``bss`` - include `BioSimSpace <https://biosimspace.org>`__ compatibility dependencies
* ``emle`` - include EMLE engine dependencies
* ``full`` - include both BSS and EMLE dependencies
* ``dev`` - all of the above plus test dependencies

If you plan to install `BioSimSpace <https://biosimspace.org>`__ on
top of :mod:`sire`, use at least the ``bss`` or ``dev`` environment.
This ensures that incompatible versions of shared dependencies are not
accidentally installed.

Compile and install
-------------------

Compilation and installation of :mod:`sire` is managed via the
`setup.py <https://github.com/openbiosim/sire/blob/devel/setup.py>`__
script.

Run

.. code-block:: bash

   $ python setup.py --help

to get help on all of the options.

Typically, you just want to compile and install :mod:`sire`. To do this,
type

.. code-block:: bash

   $ python setup.py install

This will compile the :mod:`sire` C++ libraries and then the Python
wrappers. Be patient, as compilation can take quite a while!

Once :mod:`sire` has installed, you can import it in a ``python``,
``ipython`` or ``jupyter lab`` session by typing

>>> import sire as sr

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorial <tutorial/index>`
to learn how to use :mod:`sire` or the
:doc:`quickstart guide <quickstart/index>` if you want an overview.

Please take a look at our :doc:`developer guide <contributing/development>`
for more information on how to develop and contribute new code
to :mod:`sire`.

5. Hardest install - build your own conda packages
===================================================

You can build your own :mod:`sire` conda package using
`rattler-build <https://rattler-build.readthedocs.io>`__. This is useful if
you want a package with different dependencies to the one we distribute.

A. Check out the sire source code
----------------------------------

.. code-block:: bash

   $ git clone https://github.com/openbiosim/sire -b main
   $ cd sire

You can build a package for any branch of the code. The ``main`` branch
always corresponds to the last release.

B. Install rattler-build
-------------------------

Install ``rattler-build`` via ``pixi``:

.. code-block:: bash

   $ pixi global install rattler-build

Or follow the
`rattler-build installation instructions <https://rattler-build.readthedocs.io>`__.

C. Generate the recipe
-----------------------

Generate the rattler-build recipe from the ``pixi.toml`` dependency
definitions:

.. code-block:: bash

   $ python actions/generate_recipe.py --features bss emle

This creates ``recipes/sire/recipe.yaml``. The ``--features`` flag controls
which optional dependency groups are included. Omit features to build a
lighter package:

.. code-block:: bash

   $ python actions/generate_recipe.py              # core only
   $ python actions/generate_recipe.py --features bss   # core + BioSimSpace

You can edit the generated ``recipe.yaml`` to further customise the
dependency pins if needed.

D. Build the package
---------------------

.. code-block:: bash

   $ rattler-build build --recipe recipes/sire -c conda-forge -c openbiosim/label/dev

This will take a while. The built package will be placed in the ``output/``
directory.

You can then install it directly or upload it to a conda channel.

.. note::

   A full set of tests will be run on the package after it has been built.
   Some tests may fail if you have removed dependencies. If this happens,
   you can edit the generated ``recipe.yaml`` to remove or adjust the
   test section.

