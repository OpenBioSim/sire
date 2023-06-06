=====
Units
=====

:mod:`sire` has a full units system, meaning that every value can be
properly associated with its accompanying physical unit.

The easiest way to create a value with a unit is to use the
:func:`sire.u` function.

>>> import sire as sr
>>> print(sr.u("5 nanometers"))
50 Å
>>> print(sr.u("6.5 kcal mol-1"))
6.5 kcal mol-1
>>> print(sr.u("35°"))
35°

There is a robust grammar and parser that is used to convert pretty much
all string representations of units. This can by via the long form of
the unit

>>> print(sr.u("10 kilocalories per mole"))
10 kcal mol-1

or the short form

>>> print(sr.u("10 kcal.mol-1"))
10 kcal mol-1

It recognises all of the SI prefixes, both in long and short form

>>> print(sr.u("10 microseconds"))
1e+07 ps
>>> print(sr.u("10µs"))
1e+07 ps

and has a specially-built grammar to combine units together

>>> print(sr.u("300 nm ps-2"))
3000 Å ps-2

Converting to other units
-------------------------

By default, the values are printed in internal units (see later). You
can convert them to a value of any compatible units using the
:func:`~sire.units.GeneralUnit.to` function e.g.

>>> print(sr.u("100m").to("kilometer"))
0.1
>>> print(sr.u("25 celsius").to("K"))
298.15
>>> print(sr.u("15 kcal.mol-1").to("J mol-1"))
62760.0

Supported units
---------------

:mod:`sire` supports a large number of units. You can use both the short
and long name of the unit, and can use the singular or plural form for
the long name. For example;

>>> print(sr.u("10 kcal"))
10 kcal
>>> print(sr.u("10 kilocalorie"))
10 kcal
>>> print(sr.u("10 kilocalories"))
10 kcal

Note that :mod:`sire` uses ``A`` as the short form for ``angstrom``. It is
not interpreted as ``ampere`` (you need to use the full name for ``ampere``).

This is because :mod:`sire` is focussed on molecular simulation, where
``A`` is commonly used to mean ``angstrom``. You can also use ``Å``.

Here is the full set of supported units.

*Long units*
``calorie, joule, hartree, mole, dozen, radian, degree, angstrom, meter, bohr, inch
  inches, foot, feet, yard, mile, second, minute, hour, day, week, fortnight, akma,
  dalton, gram, tonne, newton, ounce, pound, stone, hundredweight, pascal, bar, atm,
  atmosphere, psi, mmHg, kelvin, celsius, fahrenheit, amp, ampere, volt, farad,
  watt, electron, e_charge, mod_electron, faraday, coulomb, kcal_per_mol,
  kJ_per_mol``

*Short units*
``cal, J, Ha, mol, rad, °, Å, A, m, ", ', in, ft, mph, kph, s, g, N, Pa, K, °K, °C, °F, V, F, W, e, |e|, C``

.. note::

   Parsing is case-sensitive. So ``1 joule`` would parse correctly, while
   ``1 Joule`` would raise an error.

SI prefixes
-----------

All of the `SI prefixes <https://en.wikipedia.org/wiki/Metric_prefix>`__
from ``quecto`` to ``quetta`` are supported, in both long and short forms.

>>> print(sr.u("5 picometers"))
0.05 Å
>>> print(sr.u("5 pm"))
0.05 Å
>>> print(sr.u("10 megajoules"))
2390.06 kcal
>>> print(sr.u("10 MJ"))
2390.06 kcal

This includes using ``u``, ``µ`` or ``μ`` as the short version of ``micro``.

>>> print(sr.u("5 μs"))
5e+06 ps

Raising units to a power
------------------------

You can raise a unit to a power using the following symbols.

* `` `` e.g. ``m 3``, ``m3``, ``mol-1``
* ``**`` e.g. ``m**3``, ``m ** 3``, ``mol**-1``
* ``^`` e.g. ``m^3``, ``m ^ 3``, ``mol^-1``

Powers can be positive or negative, but must always be integers.
:mod:`sire` doesn't support raising units to fractional powers.

Combining units
---------------

The following symbols can be used to multiply units together.

* ``a space`` e.g. ``m s-1``, ``kcal mol-1``
* ``*`` e.g. ``m*s-1``, ``m * s-1``, ``kcal*mol-1``, ``kcal * mol-1``
* ``.`` e.g. ``m.s-1``, ``kcal.mol-1``

The following symbols can be used to divide units.

* ``/`` e.g. ``m/s``, ``m / s``, ``kcal/mol``, ``kcal / mol``
* ``per`` e.g. ``m per s``, ``kcal per mol``

Units are combined from right to left, meaning that ``kcal / mol / A**2``
is evaluated as ``kcal / (mol / A**2)``.

You can use round brackets to control the order of evaluation, e.g.
``(kcal / mol) / A**2`` would give the molar energy per square angstrom.

.. note::

   Note that ``per`` can only be used to combine individual units, e.g.
   ``kcal per mol``, not ``kcal per (mol / A**2)``. Also note that
   ``per`` is evaluated first, and only between the two units it
   is placed between. So ``kcal per mol / A**2`` will be evaluated
   as ``(kcal per mol) / A**2``.

Good rules of thumb are to use ``per`` when you want to create a derived
unit such as ``miles per hour`` or ``kcal per mol``, and to use ``/`` only
with round brackets to make sure that you get the order of evaluation
that you intend. Alternatively, do not use division at all, but instead
raise units to negative powers, e.g. ``miles hour-1`` or ``kcal.mol-1``.

Changing default units
----------------------

:mod:`sire` prints values out using default output units. You can change
these using the functions in :mod:`sire.units`, e.g.
:func:`sire.units.set_si_units` will change the output to SI units.

Changing the output units just changes how they are printed. It doens't
change their internal representation. For more info, see the section
below on ``Under the hood - GeneralUnit``.

Conversion from ``pint``
------------------------

The :func:`sire.u` function can auto-convert from other units systems.
For example, you can pass in units created via
`pint <https://pint.readthedocs.io/en/stable/>`__.

>>> import pint
>>> ureg = pint.UnitRegistry()
>>> distance = 24.0 * ureg.meter
>>> print(sr.u(distance))
2.4e+11 Å
>>> print(sr.u(distance).to(sr.u(ureg.centimeter)))
2400

Conversion from ``BioSimSpace``
-------------------------------

The :func:`sire.u` function can auto-convert from
`BioSimSpace <https://biosimspace.openbiosim.org>`__ too!

>>> import BioSimSpace as BSS
>>> import sire as sr
>>> distance = 3.5 * BSS.Units.Length.angstrom
>>> print(sr.u(distance))
3.5 Å

Conversion from other packages
------------------------------

Indeed, :func:`sire.u` can autoc-convert from any units package
that can convert to a standard units string. By default, if
:func:`sire.u` does not recognise the type, then it converts
the unit to a string, and then tries to parse it using the
in-built grammar. This should work for most cases, especially
if the other package can print units in a standard, human-readable way.

Under the hood - ``GeneralUnit``
--------------------------------

:func:`sire.u` works by parsing the string using a grammar that is built
on top of the :class:`sire.units.GeneralUnit` class. This class holds
the unit as a combination of a value and the physical dimension of the unit.

For example, ``5 m`` is ``5`` times a physical length (``L``). There
are seven physical dimensions:

1. Mass (``M``)
2. Length (``L``)
3. Time (``T``)
4. Charge (``C``)
5. temperature (``t``)
6. Quantity (``Q``)
7. Angle (``A``)

Every physical unit is a combination of these. For example, ``kcal``
is ``energy``, which is ``M2 L2 S-2`` (remember, ``E = mc2``).
Similarly, ``kcal mol-1`` is ``energy / Quantity``, so ``M2 L2 S-2 Q-1``.

The value of each physical dimension of each unit can be queried via the
functions of :class:`~sire.units.GeneralUnit`, e.g.
:func:`~sire.units.GeneralUnit.MASS` returns the power of the `M` dimension.

Internally, each dimension has a base unit which is used for scaling all
values along that dimension. The base units represent ``1.0`` for that
dimension. :mod:`sire` has base units chosen that lead to the highest precision
and best performance for the dimensional scale on which it operates
(namely the atomic scale). It uses the `AKMA <https://parmed.github.io/ParmEd/html/dimensional_analysis.html>`__
system, which is very common for molecular simulation codes.

1. Mass : ``dalton`` (chosen so ``1 g mol-1`` equals ``1.0``)
2. Length : ``angstrom``
3. Time : ``akma`` (chosen so that a time of ``1.0`` is compatible with the other units
   with no need for any scaling factors. It is approximately ``20.455 ps``)
4. Charge : ``absolute electron charge`` (chosen so a proton has charge ``1.0`` and
   an electron has charge ``-1.0``)
5. temperature : ``kelvin``
6. Quantity : ``1``
7. Angle : ``radian``

A value of ``5 meters`` is thus stored internally as ``5e10 * Length``,
while ``100 ps`` is stored internally as ``2045.48 Time``. You can
get the internal value of any unit by calling the
:func:`~sire.units.GeneralUnit.value` function, e.g.

>>> print(sr.u("100 ps").value())
2045.4828280872953

The choice of internal base units is almost invisible though, as
:mod:`sire` performs conversion from and to default output units whenever
a value is created or printed. The default output unit for time is
``picoseconds``, so ``100 ps`` when printed, will be converted from
``2045.48 Time`` to ``100 ps`` on output.

>>> print(sr.u("100 ps"))
100 ps

You can control the default output units for different functions using
the functions in :mod:`sire.units`. For example, calling
:func:`sire.units.set_si_units()` will change the default output units
to SI values.

>>> print(sr.u("10 kJ mol-1"))
2.39006 kcal mol-1
>>> sr.units.set_si_units()
>>> print(sr.u("10 kJ mol-1"))
10 kJ mol-1

.. note::

   Changing the output units does not change how the units are stored
   in :mod:`sire`. It just changes the scaling factors used to convert
   the units to/from input and output.

You can restore the default units using :func:`sire.units.set_internal_units()`

>>> sr.units.set_internal_units()
>>> print(sr.u("10 kJ mol-1"))
2.39006 kcal mol-1

You can also set individual units, e.g.
:func:`sire.units.set_mass_unit`, :func:`sire.units.set_energy_unit` etc.

Under the hood - Python to C++
------------------------------

In the Python layer, :mod:`sire` stores the value in a
:class:`~sire.units.GeneralUnit` object. This is a wrapper around the
C++ class of the same name. This C++ class is used as a temporary
intermediary to convert to templated ``PhysUnit<M,L,T,C,t,Q,A>`` objects.
These are template metaobjects, which store the physical dimension as
parameters held in the type of the C++ object (the ``M,L,T,C,t,Q,A``
parameters to the template). The object itself is just a standard
``double``, which holds the magnitude for the unit. This means that
a vector of units is just a vector of doubles. All of the unit checking
and unit code is handled via template metafunctions which are evaluated
at compile time. This means that unit types do not take up any more
space or any more compute time than plain double precision numbers.

These templated ``PhysUnit<M,L,T,C,t,Q,A>`` types are automatically
created from the C++ ``GeneralUnit`` class on function calls, and
are automatically converted back to a C++ ``GeneralUnit`` class
if a unit is returned (with this being wrapped up and exposed via
the :class:`~sire.units.GeneralUnit` Python wrapper).

In addition, the C++ ``Vector`` class, which represents a 3D point in space,
is automatically converted to hold ``Length`` types when it is queried
from the Python layer. Internally, it just holds three double precision
numbers. These are automatically converted to (or converted from) ``Length``
types when queried from C++ or Python. This minimises memory usage
and maximises compute speed.
