__all__ = []

from ..legacy import Move as _Move

from .. import use_new_api as _use_new_api

_use_new_api()


class Ensemble(_Move.Ensemble):
    """
    This class holds all information about the ensemble in which
    a simulation will be run.
    """

    def __init__(self, ensemble=None, temperature=None, pressure=None, map=None):
        from ..base import create_map

        map = create_map(map)

        from ..units import kelvin, atm

        default_temperature = 298.15 * kelvin
        default_pressure = 1.0 * atm

        if temperature is None:
            if map.specified("temperature"):
                temperature = map["temperature"].value()

        if temperature is not None:
            if not temperature.has_same_units(kelvin):
                raise TypeError(
                    "The temperature should be in units of temperature. You "
                    f"cannot use '{temperature}'"
                )

        if pressure is None:
            if map.specified("pressure"):
                pressure = map["pressure"].value()

        if pressure is not None:
            if not pressure.has_same_units(atm):
                raise TypeError(
                    "The pressure should be in units of pressure. You "
                    f"cannot use '{pressure}'"
                )

        if ensemble is None:
            if map.specified("ensemble"):
                ensemble = map["ensemble"]

                if ensemble.has_value():
                    ensemble = ensemble.value()
                else:
                    ensemble = ensemble.source()

        if ensemble is None:
            if pressure is None:
                if temperature is None:
                    ensemble = Ensemble.microcanonical()
                else:
                    ensemble = Ensemble.canonical(temperature)
            else:
                if temperature is None:
                    temperature = default_temperature

                ensemble = Ensemble.isothermal_isobaric(temperature, pressure)
        else:
            try:
                ensemble = ensemble.short_hand().lower()
            except Exception:
                ensemble = str(ensemble).lower()

            if ensemble == "nve":
                ensemble = Ensemble.microcanonical()
            elif ensemble == "nvt":
                if temperature is None:
                    temperature = default_temperature

                ensemble = Ensemble.canonical(temperature)

            elif ensemble == "npt":
                if temperature is None:
                    temperature = default_temperature

                if pressure is None:
                    pressure = default_pressure

                ensemble = Ensemble.isothermal_isobaric(temperature, pressure)

        super().__init__(ensemble)
