=================
OpenBioSim | Sire
=================

Sire is a molecular modelling framework that provides extensive
functionality to manipulate representations of biomolecular systems.

It is used as a key component of `BioSimSpace <https://biosimspace.org>`__,
and is distributed and supported as an open source community project by
`OpenBioSim <https://openbiosim.org>`__.

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

Quick Start
===========

.. toctree::
   :maxdepth: 1

   features
   install
   quickstart/index

Tutorial
========

.. toctree::
   :maxdepth: 2

   tutorial/index

Blog Posts
==========

.. toctree::
   :maxdepth: 2

   blogposts

Detailed Guides
===============

.. toctree::
   :maxdepth: 2

   cheatsheet/index

API
===

.. toctree::
   :maxdepth: 2

   api/index

Support
============

.. toctree::
   :maxdepth: 1

   support

Contributing
============

.. toctree::
   :maxdepth: 2

   contributing/index
   contributors

.. toctree::
   :maxdepth: 1

   code_of_conduct

Changelog
=========

.. toctree::
   :maxdepth: 1

   changelog

Acknowledgements
================
.. toctree::
   :maxdepth: 1

   acknowledgements

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
