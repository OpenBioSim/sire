####################################################################################################
#                                                                                                  #
#   Script for the analysis of Free energy simulations using MBAR                                  #
#                                                                                                  #
#   author: Antonia Mey <antonia.mey@ed.ac.uk>                                                     #
#                                                                                                  #
####################################################################################################
from Sire.Analysis import *
import Sire.Stream
import warnings
from Sire.Units import *

try:
    numpy = Sire.try_import("numpy")
except ImportError:
    raise ImportError(
        "Numpy is not installed. Please install numpy in order to use MBAR for your free energy analysis."
    )
try:
    scipy = Sire.try_import("scipy")
except ImportError:
    raise ImportError(
        "Scipy is not installed. Please install numpy in order to use MBAR for your free energy analysis."
    )

try:
    # This should force the installation of pymbar if it isn't
    # installed already
    _pymbar = Sire.try_import("pymbar")

    from pymbar import MBAR
    from pymbar import timeseries
except ImportError:
    import platform

    if platform.machine() in ["aarch64", "arm64"]:
        print("'pymbar' is not available on the 'aarch64' or 'arm64' platforms")

    raise ImportError(
        "'pymbar' is not installed. Please install pymbar in order to use MBAR for your free energy analysis.`"
    )

# Get the pymbar version number
_pymbar_version_no = None
try:
    _pymbar_version_no = int(_pymbar.__version__.split(".")[0])
except AttributeError:
    _pymbar_version_no = int(_pymbar.version.version.split(".")[0])

if _pymbar_version_no is None:
    raise ImportError("The pymbar version could not be determined.")

if _pymbar_version_no < 4:
    # Update the MBAR and timeseries modules to be compatible with the
    # pymbar 4 API. This requires updating the compute free energy
    # differences method to return a dictionary of results, and renaming
    # several methods.
    def compute_free_energy_differences(self, **kwargs):
        results = self.getFreeEnergyDifferences(**kwargs)
        return {
            "Delta_f": results[0],
            "dDelta_f": results[1],
        }

    MBAR.compute_free_energy_differences = compute_free_energy_differences
    MBAR.compute_overlap = MBAR.computeOverlap
    timeseries.statistical_inefficiency = timeseries.statisticalInefficiency
    timeseries.subsample_correlated_data = timeseries.subsampleCorrelatedData


import warnings


class FreeEnergies(object):
    r"""This class contains all the different pmf information
    The constructor expects subsampled MBAR and TI compatible data.
    Parameters
    ----------

    u_kln : ndarray(shape=(therm_states, therm_states, nsamples), dtype=float)
        reduced perturbed energies used for MBAR estimates
    N_K : ndarray(shape=(therm_states), dtype=int)
        number of samples per thermodynamic state
    lambda_array : ndarray(shape=(therm_states), dtype=float)
        lambda thermodynamic values
    gradients_kn : ndarray(shape=(therm_state, nsamples), dtype=float)
        reduced gradients
    """

    def __init__(self, u_kln=None, N_k=None, lambda_array=None, gradients_kn=None):
        r"""The data passed here is already subsampled"""

        self._u_kln = u_kln
        self._N_k = N_k
        self._lambda_array = lambda_array
        self._gradients_kn = gradients_kn

        # initialise results containers
        self._deltaF_mbar = None
        self._deltaF_ti = None
        self._dDeltaF_mbar = None
        self._f_k = None
        self._pmf_ti = None
        self._overlap_matrix = None
        self._pairwise_F = None

    def run_ti(self, cubic_spline=False):
        r"""Runs Thermodynamic integration free energy estimate
        Parameters
        ----------

        cubic_spline : bool
            Use cubic spline estimation instead of trapezium rule.
        """
        means = numpy.nanmean(self._gradients_kn, axis=1)
        if cubic_spline:
            NotImplementedError("Cubic Spline TI has not been implemented yet")
        else:
            self._pmf_ti = numpy.full(
                shape=(self._lambda_array.shape[0], 2), fill_value=numpy.nan
            )
            self._pmf_ti[:, 0] = self._lambda_array
            for i in range(1, self._lambda_array.shape[0]):
                self._pmf_ti[i - 1][1] = numpy.trapz(
                    means[0:i], self._lambda_array[0:i]
                )
            self._pmf_ti[-1][1] = numpy.trapz(means, self._lambda_array)
            self._deltaF_ti = numpy.trapz(means, self._lambda_array)

    def run_mbar(self, test_overlap=True):
        r"""Runs MBAR free energy estimate"""

        try:
            MBAR_obj = MBAR(self._u_kln, self._N_k, verbose=True)
            self._f_k = MBAR_obj.f_k
            results = MBAR_obj.compute_free_energy_differences()
            deltaF_ij = results["Delta_f"]
            dDeltaF_ij = results["dDelta_f"]
        except:
            solver_options = {"maximum_iterations": 10000, "verbose": True}
            solver_protocol = {"method": "BFGS", "options": solver_options}
            MBAR_obj = MBAR(self._u_kln, self._N_k, solver_protocol=(solver_protocol,))
            results = MBAR_obj.compute_free_energy_differences()
            deltaF_ij = results["Delta_f"]
            dDeltaF_ij = results["dDelta_f"]
        self._deltaF_mbar = deltaF_ij[0, self._lambda_array.shape[0] - 1]
        self._dDeltaF_mbar = dDeltaF_ij[0, self._lambda_array.shape[0] - 1]
        self._pmf_mbar = numpy.full(
            shape=(self._lambda_array.shape[0], 3), fill_value=numpy.nan
        )
        self._pmf_mbar[:, 0] = self._lambda_array
        self._pmf_mbar[:, 1] = self._f_k
        self._pmf_mbar[:, 2] = dDeltaF_ij[0]
        self._pairwise_F = numpy.full(
            shape=(self._lambda_array.shape[0] - 1, 4), fill_value=numpy.nan
        )
        self._pairwise_F[:, 0] = self._lambda_array[:-1]
        self._pairwise_F[:, 1] = self._lambda_array[1:]
        self._pairwise_F[:, 2] = numpy.diag(deltaF_ij, 1)
        self._pairwise_F[:, 3] = numpy.diag(dDeltaF_ij, 1)

        ##testing data overlap:
        if test_overlap:
            overlap_matrix = MBAR_obj.compute_overlap()
            self._overlap_matrix = overlap_matrix["matrix"]

    @property
    def pmf_ti(self):
        return self._pmf_ti

    @property
    def pmf_mbar(self):
        return self._pmf_mbar

    @property
    def deltaF_ti(self):
        return self._deltaF_ti

    @property
    def deltaF_mbar(self):
        return self._deltaF_mbar

    @property
    def errorF_mbar(self):
        return self._dDeltaF_mbar

    @property
    def overlap_matrix(self):
        return self._overlap_matrix

    @property
    def pairwise_F(self):
        return self._pairwise_F


class SubSample(object):
    r"""This class subsamples data based on the timeseries analysis or percentage of data ready for pmf use
    Parameters
    ----------
    gradients_kn : ndarray(shape=(therm_state, nsamples), dtype=float)
        reduced gradients
    energies : ndarray(shape=(therm_state, nsamples), trype=float)
        potential energies used to find statisitical inefficiency
    u_kln : ndarray(shape=(therm_states, therm_states, nsamples), dtype=float)
        reduced perturbed energies used for MBAR estimates
    N_K : ndarray(shape=(therm_states), dtype=int)
        number of samples per thermodynamic state
    lambda_array : ndarray(shape=(therm_states), dtype=float)
        lambda thermodynamic values
    percentage : int [0,100]
        percentage of the data that should be discarded from the beginning of the simulation
    subsample : string
        string idenfier for subsampling method, default='timeseries' from timeseries module in MBAR
    """

    def __init__(
        self,
        gradients_kn,
        energies,
        u_kln,
        N_k,
        percentage=100,
        subsample=True,
    ):
        self._gradients_kn = gradients_kn
        self._N_k = N_k
        self._energies_kn = energies
        self._u_kln = u_kln
        self._subsampled_u_kln = None
        self._subsampled_N_k_energies = None
        self._subsampled_N_k_gradients = None
        self._subsampled_grad_kn = None
        self._subsampled_energies_kn = None

        if N_k.shape[0] != u_kln.shape[0]:
            RuntimeError(
                "The number of thermodynamic states must be the same in u_kln and N_k!"
                "u_kln has size %d and N_k has size %d" % (u_kln.shape[0], N_k.shape[0])
            )
        self.subsample = subsample
        self.percentage = percentage
        assert (
            percentage > 0.0 and percentage <= 100.0
        ), "You must provide a percentage between 0 and 100"

    def subsample_gradients(self):
        r"""method to subsample gradients and get a better estiamte."""
        if self.percentage == 100 and not self.subsample:
            warnings.warn(
                "You are not subsampling your data according to the statistical inefficiency nor are "
                "you discarding initial data. Please set percentage to another value than 100!"
            )
        percentage_removal = (self._N_k * (1 - self.percentage / 100.0)).astype("int32")
        self._subsampled_N_k_gradients = self._N_k - percentage_removal
        N_max = int(numpy.max(self._subsampled_N_k_gradients))
        self._subsampled_grad_kn = numpy.full(
            shape=(self._N_k.shape[0], N_max), fill_value=numpy.nan
        )
        for p in range(percentage_removal.shape[0]):
            start = percentage_removal[p]
            finish = percentage_removal[p] + N_max
            self._subsampled_grad_kn[p, :] = self._gradients_kn[p, start:finish]
        if N_max <= 50:
            warnings.warn(
                "You have reduced your data to less than 50 samples, the results from these might not "
                "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option."
            )
        # if subsampling is percentage, then we are done here, otherwise we will now subsample according to timeseries

        if self.subsample:
            print("#Subsampling gradients according to statistical inefficiency")
            # first we compute statistical inefficiency
            self._gradients_kn = self._subsampled_grad_kn.copy()
            self._N_k = self._subsampled_N_k_gradients.copy()

            g_k = numpy.full(shape=(self._gradients_kn.shape[0]), fill_value=numpy.nan)
            self._subsampled_N_k_gradients = numpy.full(
                shape=(self._gradients_kn.shape[0]), fill_value=numpy.nan
            )
            for i in range(g_k.shape[0]):
                g_k[i] = timeseries.statistical_inefficiency(self._gradients_kn[i, :])
            g = int(numpy.max(g_k))
            # now we need to figure out what the indices in the data are for subsampling
            indices_k = []
            for i in range(g_k.shape[0]):
                indices_k.append(
                    timeseries.subsample_correlated_data(self._gradients_kn[i, :], g=g)
                )
                self._subsampled_N_k_gradients[i] = len(indices_k[i])
            N_max = int(numpy.max(self._subsampled_N_k_gradients))
            if N_max <= 50:
                warnings.warn(
                    "You have reduced your data to less than 50 samples, the results from these might not "
                    "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option."
                )
            self._subsampled_grad_kn = numpy.full(
                [self._gradients_kn.shape[0], N_max], fill_value=numpy.nan
            )
            for k in range(self._gradients_kn.shape[0]):
                self._subsampled_grad_kn[k, :] = self._gradients_kn[k, indices_k[k]]

    def subsample_energies(self):
        r"""This subsamples u_kln according to percentage, i.e. remove initial equilibration data and then can additionally subsample according to timeseries"""
        # removing percent
        if self.percentage == 100 and not self.subsample:
            warnings.warn(
                "You are not subsampling your data according to the statistical inefficiency nor are "
                "you discarding initial data. Please set percentage to another value than 100!"
            )

        percentage_removal = (self._N_k * (1 - self.percentage / 100.0)).astype("int32")
        self._subsampled_N_k_energies = self._N_k - percentage_removal
        N_max = int(numpy.max(self._subsampled_N_k_energies))
        self._subsampled_u_kln = numpy.full(
            shape=(self._N_k.shape[0], self._N_k.shape[0], N_max), fill_value=numpy.nan
        )
        self._subsampled_energies_kn = numpy.full(
            shape=(self._N_k.shape[0], N_max), fill_value=numpy.nan
        )
        for k in range(0, self._N_k.shape[0]):
            self._subsampled_u_kln[k] = self._u_kln[
                k, :, percentage_removal[k] : percentage_removal[k] + N_max
            ]
            self._subsampled_energies_kn[k] = self._energies_kn[
                k, percentage_removal[k] : percentage_removal[k] + N_max
            ]
        if N_max <= 50:
            warnings.warn(
                "You have reduced your data to less than 50 samples, the results from these might not "
                "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option."
            )

        # Now we are doing some additional subsampling according to timeseries analysis
        if self.subsample:
            print(
                "#Subsampling energies according to statistical inefficiency for pymbar"
            )

            self._u_kln = self._subsampled_u_kln.copy()
            self._N_k = self._subsampled_N_k_energies.copy()
            self._energies_kn = self._subsampled_energies_kn.copy()
            # first we compute statistical inefficiency
            g_k = numpy.full(shape=(self._energies_kn.shape[0]), fill_value=numpy.nan)
            for i in range(g_k.shape[0]):
                g_k[i] = timeseries.statistical_inefficiency(
                    self._energies_kn[i, percentage_removal[i] :]
                )
            g = numpy.max(g_k)
            # now we need to figure out what the indices in the data are for subsampling
            indices_k = []
            self._subsampled_N_k_energies = numpy.full(
                shape=(self._energies_kn.shape[0]), fill_value=numpy.nan
            )
            for i in range(g_k.shape[0]):
                indices_k.append(
                    timeseries.subsample_correlated_data(self._energies_kn[i, :], g=g)
                )
                self._subsampled_N_k_energies[i] = len(indices_k[i])
            # self._subsampled_N_k_energies = (numpy.ceil(self._N_k / g)).astype(int)
            N_max = int(numpy.max(self._subsampled_N_k_energies))
            if N_max <= 50:
                warnings.warn(
                    "You have reduced your data to less than 50 samples, the results from these might not "
                    "be trustworthy. If you don't want to add more samples consider rerunning the analysis using the percentage option."
                )
            self._subsampled_u_kln = numpy.full(
                [
                    self._gradients_kn.shape[0],
                    self._gradients_kn.shape[0],
                    N_max,
                ],
                fill_value=numpy.nan,
            )
            for k in range(self._gradients_kn.shape[0]):
                self._subsampled_u_kln[k, :, :] = self._u_kln[
                    k, :, indices_k[k]
                ].transpose()

    @property
    def u_kln(self):
        return self._subsampled_u_kln

    @property
    def gradients_kn(self):
        return self._subsampled_grad_kn

    @property
    def N_k_energies(self):
        return self._subsampled_N_k_energies

    @property
    def N_k_gradients(self):
        return self._subsampled_N_k_gradients
