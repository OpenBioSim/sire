==============
Part 8 - QM/MM
==============

QM/MM is a method that combines the accuracy of quantum mechanics with the
speed of molecular mechanics. In QM/MM, a small region of the system is treated
at the quantum mechanical level, while the rest of the system is treated at the
molecular mechanical level. This allows us to perform accurate calculations on
the region of interest, while still being able to simulate the rest of the system
at a much lower computational cost. Due to the recent development of cheap and
accurate machine learning based QM models, there has been a resurgence of
interest in QM/MM methods. In this tutorial we will show how to perform
QM/MM simulations using ``sire``.

.. toctree::
   :maxdepth: 1

   part08/01_intro
   part08/02_emle
   part08/03_adp_pmf
   part08/04_diels_alder
