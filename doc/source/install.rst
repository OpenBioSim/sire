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

You compile :mod:`sire` into an existing anaconda / miniconda environment.
Please create and activate an environment, e.g. by following
`the instructions <_Install_miniforge>` to install a fresh ``miniforge`` and
then creating and activating Python 3.11 environment called
``openbiosim``.

Next, download the source code. You could download the latest development
version of :mod:`sire` by typing;

.. code-block:: bash

   $ git clone https://github.com/openbiosim/sire

This will download into a directory called :mod:`sire`. Navigate into
this directory (e.g. ``cd sire``).

.. note::

   This will fail if ``git`` is not installed on your computer.
   You can easily install ``git`` using ``conda``, e.g.
   run ``conda install git``.

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

Compilation and installation of :mod:`sire` is managed via the
`setup.py <https://github.com/openbiosim/sire/blob/devel/setup.py>`__
script.

Run

.. code-block:: bash

   $ python setup.py --help

to get a help on all of the options.

Typically, you just want to compile and install :mod:`sire`. To do this,
type

.. code-block:: bash

   $ python setup.py install

This will download and install all of the dependencies via ``conda``. It will then compile
the :mod:`sire` C++ libraries, and then the Python wrappers. Be patient,
as compilation can take quite a while!

.. note::

   You need to have Visual Studio C++ (2017 or newer) installed to compile on Windows.
   The easiest way to do this is to install the free
   `Visual Studio 2022 Community Edition <https://visualstudio.microsoft.com/vs/community>`__.
   Make sure to install "Desktop development with C++",
   including the options "MSVC v143 - VS 2022 C++ x64/x86 build tools (v14.30)",
   "C++ CMake tools for Windows", and at least one of "Windows 11 SDK" and/or
   "Windows 10 SDK" (any version will do). You can, optionally, install the
   older C++ compilers too, e.g. "MSVC v142 - VS 2019 C++ x64/x86 build tools (v14.29)",
   and/or "MSVC v141 - VS 2017 C++ x64/x86 build tools (v14.16)". Currently
   only the X64 compilers have been tested - we are interested to try
   Windows/ARM64 once more of the dependencies are available.

If you plan to install `BioSimSpace <https://biosimspace.org>`__ on
top of :mod:`sire`, then you should install using;

.. code-block:: bash

   $ python setup.py --install-bss-deps install

This will use ``conda`` to download and install all of
BioSimSpace's dependencies as well. This ensures that incompatible versions
of shared dependencies are not accidentally installed.

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

5. Hardest install - build your own custom conda packages
=========================================================

The :mod:`sire` conda packages that we build have a lot of dependencies that
may conflict with your own environment. This is because we build :mod:`sire`
to be compatible with the latest version of `BioSimSpace <https://biosimspace.openbiosim.org>`__,
which itself optionally depends on a large number of simulation packages.

You can build your own :mod:`sire` conda package that has fewer dependencies,
or which is compatible with the packages already installed in your conda
environment. There are a few steps you need to complete.

A. Define your runtime environment
----------------------------------

The first step is to describe the desired runtime environment for the package.
The easiest way to do this is to create that environment, e.g. by installing
the packages you want, and to then create an ``environment.yml`` file
that describes that environment. You can do this by running

.. code-block:: bash

   $ conda env export -f environment.yml

This will create an environment file called ``environment.yml``
that creates pins for the exact version of all of the packages installed
in your environment.

If you want, you can edit this file to add or remove pins. Simply delete
lines describing the version of packages that you don't need pinned,
add new lines if there are additional packages that you do want pinned,
or even update the version number of the pins if you can allow more
flexibility for the installation.

B. Check out the sire source code
---------------------------------

The next step is to check out the :mod:`sire` source code (if you haven't
already).

.. code-block:: bash

   $ git clone https://github.com/openbiosim/sire -b main

This checks the ``main`` branch of the code out into a directory called
``sire``. You can build a package for any branch of the code. Typically,
you will want to choose the ``main`` branch, as this always corresponds to the
last release. You can checkout the ``main`` branch by changing into the
``sire`` directory and running;

.. code-block:: bash

   $ git checkout main

C. Create the conda build environment
-------------------------------------

While you could build in your existing environment, it is cleaner to
build in a dedicated build environment. Here, we will create a build
environment called ``build_sire``. You can use any name you want.

.. code-block:: bash

   $ conda env create -n build_sire -f environment.yml

Activate that environment

.. code-block:: bash

   $ conda activate build_sire

And then install the tools needed to run conda-build

.. code-block:: bash

   $ conda install -y -c conda-forge boa anaconda-client packaging=21 pip-requirements-parser

D. Create the conda recipe
--------------------------

Next, we need to create the conda recipe to build the package. We do this by
running the script ``actions/update_recipe.py``. You can add the path to
your ``environment.yml`` file as an argument. This tells the script to
create a recipe that includes all of the pins in the ``environment.yml``.
For example;

.. code-block:: bash

   $ python actions/update_recipe.py environment.yml

would create the recipe using the pins in ``environment.yml`` (assuming this
file was in the current directory).

The recipe is written to ``recipes/sire/meta.yaml``. You can (optionally)
edit the pins in this file too, if you want to do some fine-tuning.

.. note::

   You may need to edit the recipe to fix version inconsistencies.
   This is especially the case for ``rdkit`` - you need to to make
   sure that if you specify a version for ``rdkit`` in your
   ``environment.yml`` that you also use the same version
   for the ``rdkit-devel`` package.

E. Building the package
-----------------------

You can now run ``conda-build`` to create the package.

.. code-block:: bash

   $ conda build -c conda-forge -c openbiosim/label/dev recipes/sire

This will take a while. At the end, it will print out the location of the
sire conda package, e.g.

.. note::

   The above command assumes that you don't need any other channels included
   to install all of the packages included in your ``environment.yml``.
   The ``actions/update_recipe.py`` script will print out the correct
   ``conda build`` command at the end, which includes any extra
   channels that are needed.

::

   # To have conda build upload to anaconda.org automatically, use
   # conda config --set anaconda_upload yes
   anaconda upload \
       /path/to/miniforge/envs/build_sire/conda-bld/osx-64/sire-2023.3.0-py310hf95ea87_25.tar.bz2
   anaconda_upload is not set.  Not uploading wheels: []

   INFO :: The inputs making up the hashes for the built packages are as follows:
   {
     "sire-2023.3.0-py310hf95ea87_25": {
       "recipe": {
         "c_compiler": "clang",
         "cxx_compiler": "clangxx",
         "numpy": "1.22",
         "target_platform": "osx-64"
       }
     }
   }

In this case, you can see that the package is the file
``/path/to/miniforge/envs/build_sire/conda-bld/osx-64/sire-2023.3.0-py310hf95ea87_25.tar.bz2``.

Copy this conda package to wherever you need (e.g. into a channel, upload
to conda, etc.).

.. note::

   A full set of tests will be run on the package after it has been built.
   Some of these tests may fail if you have edited the recipe to remove
   some of the dependencies. If this happens, you can decide to ignore
   the tests, e.g. by removing them from the conda recipe (``meta.yml``)
   or by just copying the file that is produced and has been placed into
   the ``conda-bld/broken`` directory.

You can then install it, either via the channel you've uploaded to, or by
directly running ``conda install`` on the package file itself.

