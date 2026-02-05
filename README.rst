======================================
`Sire <https://sire.openbiosim.org>`__
======================================

Citation
========

If you use sire in your work, please cite the
`following paper <https://doi.org/10.1063/5.0200458>`__:

.. code-block:: bibtex

    @article{10.1063/5.0200458,
        author = {Woods, Christopher J. and Hedges, Lester O. and Mulholland, Adrian J. and Malaisree, Maturos and Tosco, Paolo and Loeffler, Hannes H. and Suruzhon, Miroslav and Burman, Matthew and Bariami, Sofia and Bosisio, Stefano and Calabro, Gaetano and Clark, Finlay and Mey, Antonia S. J. S. and Michel, Julien},
        title = "{Sire: An interoperability engine for prototyping algorithms and exchanging information between molecular simulation programs}",
        journal = {The Journal of Chemical Physics},
        volume = {160},
        number = {20},
        pages = {202503},
        year = {2024},
        month = {05},
        abstract = "{Sire is a Python/C++ library that is used both to prototype new algorithms and as an interoperability engine for exchanging information between molecular simulation programs. It provides a collection of file parsers and information converters that together make it easier to combine and leverage the functionality of many other programs and libraries. This empowers researchers to use sire to write a single script that can, for example, load a molecule from a PDBx/mmCIF file via Gemmi, perform SMARTS searches via RDKit, parameterize molecules using BioSimSpace, run GPU-accelerated molecular dynamics via OpenMM, and then display the resulting dynamics trajectory in a NGLView Jupyter notebook 3D molecular viewer. This functionality is built on by BioSimSpace, which uses sire’s molecular information engine to interconvert with programs such as GROMACS, NAMD, Amber, and AmberTools for automated molecular parameterization and the running of molecular dynamics, metadynamics, and alchemical free energy workflows. Sire comes complete with a powerful molecular information search engine, plus trajectory loading and editing, analysis, and energy evaluation engines. This, when combined with an in-built computer algebra system, gives substantial flexibility to researchers to load, search for, edit, and combine molecular information from multiple sources and use that to drive novel algorithms by combining functionality from other programs. Sire is open source (GPL3) and is available via conda and at a free Jupyter notebook server at https://try.openbiosim.org. Sire is supported by the not-for-profit OpenBioSim community interest company.}",
        issn = {0021-9606},
        doi = {10.1063/5.0200458},
        url = {https://doi.org/10.1063/5.0200458},
        eprint = {https://pubs.aip.org/aip/jcp/article-pdf/doi/10.1063/5.0200458/19969848/202503\_1\_5.0200458.pdf},
    }

About
=====

Sire is a molecular modelling framework that provides extensive
functionality to manipulate representations of biomolecular systems.

It is used as a key component of `BioSimSpace <https://biosimspace.org>`__,
and is distributed and supported as an open source community project by
`OpenBioSim <https://openbiosim.org>`__.

For more information about how to use Sire, and about application
built with Sire, please `visit the Sire website <https://sire.openbiosim.org>`__.

* `Features <https://sire.openbiosim.org/features.html>`__
* `Quick start guide <https://sire.openbiosim.org/quickstart/index.html>`__
* `Tutorial <https://sire.openbiosim.org/tutorial/index.html>`__

Installation
============

The easiest way to install Sire is using our `conda channel <https://anaconda.org/openbiosim/repo>`__.
Sire is built using dependencies from `conda-forge <https://conda-forge.org/>`__,
so please ensure that the channel takes strict priority. We recommend using
`miniforge3 <https://github.com/conda-forge/miniforge#miniforge3>`__.

To create a new environment:

.. code-block:: bash

    conda create -n openbiosim "python<3.13"
    conda activate openbiosim
    conda install -c conda-forge -c openbiosim sire

To install the latest development version you can use:

.. code-block:: bash

    conda create -n openbiosim-dev "python<3.13"
    conda activate openbiosim-dev
    conda install -c conda-forge -c openbiosim/label/dev sire

Installing older versions
-------------------------

You can install a specific version of sire by specifying the version number
in the conda install command, e.g.

.. code-block:: bash

    conda install -c conda-forge -c openbiosim sire==2024.1.0

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

    conda install -c https://openbiosim.blob.core.windows.net/archive sire==2023.2.3

Installation from source
------------------------

However, as you are here, it is likely you want to download the latest,
greatest version of the code, which you will need to compile. To compile
sire, you need a git client to download the source and
`pixi <https://pixi.sh>`__ to manage the build environment.

First, clone the sire source code and change into the directory:

.. code-block:: bash

    git clone https://github.com/OpenBioSim/sire
    cd sire

Next, use pixi to create and activate the development environment. This
will install all required dependencies, including compilers:

.. code-block:: bash

    pixi install -e dev
    pixi shell -e dev

Now compile and install sire:

.. code-block:: bash

    python setup.py install

A small word of warning, the compilation can easily take over an hour!

Other pixi environments are available depending on your needs:

* ``pixi install -e default`` - core sire dependencies only
* ``pixi install -e obs`` - include downstream OpenBioSim package dependencies
* ``pixi install -e emle`` - include EMLE engine dependencies
* ``pixi install -e dev`` - all of the above plus test dependencies

Support and Development
=======================

Bugs, Comments, Questions
-------------------------
For bug reports/suggestions/complaints please file an issue on
`GitHub <http://github.com/OpenBioSim/sire/issues>`__.

Developers guide
----------------
Please `visit the website <https://sire.openbiosim.org>`__ for information on how to
develop applications using sire.

GitHub actions
--------------
Since sire is quite large, a build can take quite long and might not be neccessary
if a commit is only fixing a couple of typos. Simply add ``ci skip``
to your commit message and GitHub actions will not invoke an autobuild.

Note that every time you commit to devel, it will trigger a build of sire,
full testing, construction of a Conda package and upload to our Anaconda
channel. Please think twice before committing directly to devel. You should
ideally be working in a _feature_ branch, and only commit to devel once you are
happy the code works on your branch. Use ``ci skip`` until you are happy that
you want to trigger a full build, test and deployment. This full pipeline will
take several hours to complete.

Have fun :-)
