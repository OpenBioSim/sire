__all__ = ["to_alchemlyb"]


def to_alchemlyb(energy_trajectories, temperature=None, energy_unit="kcal/mol"):
    """
    Convert the passed list of energy trajectories into
    a single alchemlyb-compatible DataFrame, ready for
    free energy analysis via alchemlyb.

    Parameters
    ----------

    energy_trajectories : list of sire.maths.EnergyTrajectory objects,
                          list of filenames of s3 files,
                          globbed expression for list of filenames etc.
        A list of EnergyTrajectory objects, each containing the
        energy trajectory for a single simulation at a single
        lambda value.

    temperature : temperature, optional
        The temperature of the simulation. If not provided,
        the temperature will be taken from the values in
        each EnergyTrajectory

    energy_unit: str
        Whichever of the alchemlyb energy units you want the output
        DataFrame to use. This is in alchemlyb format, e.g.
        `kcal/mol`, `kJ/mol`, or `kT`

    Returns
    -------

    pandas.DataFrame
        A single DataFrame containing the energy trajectories
        from all simulations, ready for free energy analysis
        via alchemlyb.
    """
    if type(energy_trajectories) is str:
        # this could be a globbed file path
        import glob

        energy_trajectories = glob.glob(energy_trajectories)
    elif type(energy_trajectories) is not list:
        energy_trajectories = [energy_trajectories]

    # convert all of the energy trajectories to pandas DataFrames
    dataframes = {}

    for energy_trajectory in energy_trajectories:
        if type(energy_trajectory) is str:
            # load the energy trajectory from the s3 file
            from ..stream import load

            energy_trajectory = load(energy_trajectory)

        df = energy_trajectory.to_alchemlyb(
            temperature=temperature, energy_unit=energy_unit
        )

        # get the lambda value of the first row
        try:
            lambda_value = df.index[0][1]
        except IndexError:
            # this is a corrupt energy trajectory...
            from ..utils import Console

            Console.warning("Corrupt energy trajectory detected. Skipping...")
            continue

        if lambda_value in dataframes:
            dataframes[lambda_value].append(df)
        else:
            dataframes[lambda_value] = [df]

    # get the list of lambda values and sort...
    lambda_values = list(dataframes.keys())
    lambda_values.sort()

    # now create a list of dataframes in lambda order
    dfs = []

    for lambda_value in lambda_values:
        dfs.extend(dataframes[lambda_value])

    if len(dfs) == 0:
        return None
    elif len(dfs) == 1:
        return dfs[0]
    else:
        attrs = dfs[0].attrs

        # now concatenate the dataframes
        import pandas as pd

        combined = pd.concat(dfs)
        combined.attrs = attrs

        return combined
