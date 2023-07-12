from typing import List
from lmfit import Parameters, minimize
from lmfit.minimizer import MinimizerResult
from typing import Dict, Callable, Tuple
from numpy import isin, log, exp
from scipy.integrate import odeint
import numpy as np

from EnzymePynetics.core.modelresult import ModelResult
from EnzymePynetics.core.parameter import Parameter
from EnzymePynetics.core.correlation import Correlation


class KineticModel:
    def __init__(
        self,
        name: str,
        model: Callable,
        params: list,
        kcat_initial: float,
        Km_initial: float,
        y0: List[tuple],
        enzyme_inactivation: bool,
    ) -> None:
        self.name = name
        self.model = model
        self.params = params
        self.enzyme_inactivation = enzyme_inactivation
        self.y0 = y0
        self.kcat_initial = kcat_initial
        self.Km_initial = Km_initial
        self.parameters = self._set_parameters(params)
        self._fit_result: MinimizerResult = None
        self.result: ModelResult = None

    def _set_parameters(self, params: list) -> Parameters:
        """Initializes lmfit parameters, based on provided initial parameter guesses.

        Args:
            params (list): Parameter keys.

        Returns:
            Parameters: lmfit parameters with initil values and bounds.
        """

        parameters = Parameters()

        parameters.add(
            "k_cat",
            value=self.kcat_initial,
            min=self.kcat_initial / 100,
            max=self.kcat_initial * 100,
        )
        parameters.add(
            "Km",
            value=self.Km_initial,
            min=self.Km_initial / 100,
            max=self.Km_initial * 10000,
        )

        if "K_iu" in params:
            parameters.add("K_iu", value=0.1, min=0.0001, max=1000)
        if "K_ic" in params:
            parameters.add("K_ic", value=0.1, min=0.0001, max=1000)
        if self.enzyme_inactivation:
            parameters.add("K_ie", value=0.01, min=0.0001, max=0.9999)
        if "k_inact" in params:
            parameters.add("k_inact", value=0.01, min=0.0001, max=0.9999)

        return parameters

    def integrate(
        self, parameters: Parameters, time: list, y0: tuple
    ) -> np.ndarray:
        """Integrates model based on parameters for a given time array and initial conditions.

        Args:
            parameters (Parameters): lmfit parameters
            time (ndarray): time array for integration
            y0s (List[tuple]): initial conditions for the species in the model.

        Returns:
            result (np.ndarray): integrated model over given time.
        """

        result = np.array([odeint(func=self.model, y0=y, t=t, args=(
            parameters, self.enzyme_inactivation)) for y, t in zip(y0, time)])
        return result

    def residuals(
        self,
        parameters: Parameters,
        time: np.ndarray,
        y0s: List[tuple],
        ydata: np.ndarray,
    ) -> np.ndarray:
        """Calculates residuals between integrated model and measured data (substrate).

        Args:
            parameters (Parameters): LmFit parameters
            time (ndarray): time array, corresponding to ydata.
            y0s (List[tuple]): initial conditions of modeled species
            ydata (ndarray): measured substrate data, corresponding to time data.

        Returns:
            residuals (np.ndarray): List of substrate residuals.
        """

        y0s = np.array(y0s)
        model = self.integrate(parameters, time, y0s)
        residuals = model[:, :, 0] - ydata  # fitting to substrate data
        return residuals.flatten()

    def fit(self, ydata: np.ndarray, time: np.ndarray) -> MinimizerResult:
        """Fit model to substrate data.

        Args:
            ydata (ndarray): Experimental substrate data
            time (ndarray): Time array corresponding to measurement data

        Returns:
            MinimizerResult: least-squares minimization result.
        """

        fit_result: MinimizerResult = minimize(
            self.residuals, self.parameters, args=(time, self.y0, ydata)
        )
        self._fit_result = fit_result
        self.result = self._get_model_results(fit_result)

    def _calcualte_RMSD(self, lmfit_result: MinimizerResult) -> float:
        """Calculates root mean square deviation (RMSD) between model and experimental data.

        Args:
            lmfit_result (MinimizerResult): lmfit fitting result

        Returns:
            float: RMSD
        """
        residuals = lmfit_result.residual
        return np.sqrt(1 / residuals.size * np.sum(residuals**2))

    def _get_model_results(self, lmfit_result: MinimizerResult) -> ModelResult:
        """Extracts fitting parameters and statistics from the lmfit result.

        Args:
            lmfit_result (MinimizerResult): Result from lmfit minimization.

        Returns:
            ModelResult: Result parameters and statistics.
        """

        # Write lmfit results to ModelResult
        model_result = ModelResult()
        model_result.name = self.name
        model_result.equation = "#TODO"
        model_result.fit_success = lmfit_result.success

        if model_result.fit_success:
            model_result.AIC = lmfit_result.aic
            model_result.BIC = lmfit_result.bic
            model_result.RMSD = self._calcualte_RMSD(lmfit_result)

            # Get parameters and correlations between parameters
            parameters = []
            for key, value in lmfit_result.params.items():
                correlations = []
                for corr_key, corr_value in value.correl.items():
                    correlations.append(
                        Correlation(parameter=corr_key, value=corr_value)
                    )

                parameters.append(
                    Parameter(
                        name=key,
                        value=value.value,
                        standard_deviation=value.stderr,
                        upper_limit=value.max,
                        lower_limit=value.min,
                        correlations=correlations,
                    )
                )

            model_result.parameters = parameters

        return model_result
