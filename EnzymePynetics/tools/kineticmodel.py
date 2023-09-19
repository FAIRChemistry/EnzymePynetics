from ast import Call
from typing import List, Callable, Tuple
from lmfit import Parameters, minimize
from lmfit.minimizer import MinimizerResult
from scipy.integrate import odeint
import numpy as np
import sympy as sp

from EnzymePynetics.core.modelresult import ModelResult
from EnzymePynetics.core.parameter import Parameter
from EnzymePynetics.core.correlation import Correlation
from EnzymePynetics.enums.params import Params


class KineticModel:
    def __init__(
        self,
        name: str,
        substrate_rate_law: str,
        params: list,
        kcat_initial: float,
        Km_initial: float,
        enzyme_rate_law: str = None,
    ) -> None:
        self.name = name
        self.substrate_rate_law = substrate_rate_law
        self.enzyme_rate_law = enzyme_rate_law
        self.params = params
        self.kcat_initial = kcat_initial
        self.Km_initial = Km_initial
        self.substrate_callable: Callable = self._get_callable(substrate_rate_law)
        self.enzyme_callable: Callable = self._get_callable(enzyme_rate_law)
        self._fit_result: MinimizerResult = None
        self.result: ModelResult = None
        self.parameters = self._set_parameters(params)

    def _get_callable(self, equation: str) -> Callable:
        if isinstance(equation, str):
            # local_sympy_dict ensures that 'product' in the string expression is treated as
            # a symbol instead of a function
            local_sympy_dict = {"product": sp.Symbol("product")}
            expr_substrate = sp.parse_expr(equation, local_sympy_dict)
            free_symbols = list(expr_substrate.free_symbols)

            return sp.lambdify(free_symbols, expr_substrate)

        elif equation is None:
            return None

        else:
            raise TypeError(
                f"Equation of type {type(equation)} cannot be interpreted as string."
            )

    def _set_parameters(self, params: list) -> Parameters:
        """Initializes lmfit parameters, based on provided initial parameter guesses.

        Args:
            params (list): Parameter keys.

        Returns:
            Parameters: lmfit parameters with initial values and bounds.
        """

        parameters = Parameters()

        parameters.add(
            name=Params.k_cat.value,
            value=self.kcat_initial,
            min=self.kcat_initial / 50,
            max=self.kcat_initial * 50,
        )
        parameters.add(
            name=Params.K_m.value,
            value=self.Km_initial,
            min=self.Km_initial * 0.5,
            max=self.Km_initial * 5000,
        )

        if Params.K_iu.value in params:
            parameters.add(Params.K_iu.value, value=0.1, min=0.0001, max=1000)
        if Params.K_ic.value in params:
            parameters.add(Params.K_ic.value, value=0.1, min=0.0001, max=1000)
        if self.enzyme_rate_law:
            parameters.add(Params.k_ie.value, value=0.01, min=0.00001, max=0.9999)

        return parameters

    @staticmethod
    def model(w0, t, params, substrate_eq: Callable, enzyme_eq: Callable = None):
        species_keys = ["substrate", "enzyme", "product", "inhibitor"]
        observables_dict = dict(zip(species_keys, w0))

        params_dict = params.valuesdict()

        combined_dict = observables_dict | params_dict
        subtrate_eq_vars = substrate_eq.__code__.co_varnames
        substrate_dict = {k: combined_dict[k] for k in subtrate_eq_vars}

        d_substrate = substrate_eq(**substrate_dict)

        if enzyme_eq:
            enzyme_eq_vars = enzyme_eq.__code__.co_varnames
            enzyme_dict = {k: combined_dict[k] for k in enzyme_eq_vars}
            d_enzyme = enzyme_eq(**enzyme_dict)
        else:
            d_enzyme = 0

        d_product = -d_substrate
        d_inhibitor = 0

        return (d_substrate, d_enzyme, d_product, d_inhibitor)

    def integrate(self, parameters: Parameters, time: list, y0: tuple) -> np.ndarray:
        """Integrates model based on parameters for a given time array and initial conditions.

        Args:
            parameters (Parameters): lmfit parameters
            time (ndarray): time array for integration
            y0s (List[tuple]): initial conditions for the species in the model.

        Returns:
            result (np.ndarray): integrated model over given time.
        """

        result = np.array(
            [
                odeint(
                    func=self.model,
                    y0=y,
                    t=t,
                    args=(parameters, self.substrate_callable, self.enzyme_callable),
                )
                for y, t in zip(y0, time)
            ]
        )
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

        model = self.integrate(parameters, time, y0s)
        residuals = model[:, :, 0] - ydata  # fitting to substrate data
        return residuals.flatten()

    def fit(
        self, ydata: np.ndarray, time: np.ndarray, y0s: np.ndarray
    ) -> MinimizerResult:
        """Fit model to substrate data"""

        # y0s = np.array(y0s) y0s needs to be checked if it is ndarray?
        fit_result = minimize(
            self.residuals, self.parameters, args=(time, y0s, ydata), nan_policy="omit"
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
        model_result.fit_success = lmfit_result.success
        model_result.equations.append(self.substrate_rate_law)
        if self.enzyme_rate_law:
            model_result.equations.append(self.enzyme_rate_law)

        if model_result.fit_success:
            model_result.AIC = lmfit_result.aic
            model_result.BIC = lmfit_result.bic
            model_result.RMSD = self._calcualte_RMSD(lmfit_result)

            # Get parameters and correlations between parameters
            parameters = []
            for key, value in lmfit_result.params.items():
                correlations = []
                try:
                    for corr_key, corr_value in value.correl.items():
                        correlations.append(
                            Correlation(parameter=corr_key, value=corr_value)
                        )
                except AttributeError:
                    pass

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
