from typing import List
from lmfit import Parameters, minimize
from lmfit.minimizer import MinimizerResult
from typing import Dict, Callable, Tuple
from numpy import log, exp
from scipy.integrate import odeint

import numpy as np


class KineticModel():

    def __init__(self,
                 name: str,
                 model: Callable,
                 params: list,
                 kcat_initial: float,
                 Km_initial: float,
                 y0: List[tuple],
                 enzyme_inactivation: bool = False
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

    def _set_parameters(self, params: list) -> Parameters:
        """Initializes lmfit parameters, based on provided initial parameter guesses.

        Args:
            params (list): Parameter keys.

        Returns:
            Parameters: lmfit parameters with initil values and bounds.
        """

        parameters = Parameters()

        parameters.add('k_cat', value=self.kcat_initial,
                        min=self.kcat_initial/100, max=self.kcat_initial*100)
        parameters.add('Km', value=self.Km_initial, min=self.Km_initial/100,
                        max=self.Km_initial*10000)

        if "K_iu" in params:
            parameters.add("K_iu", value=0.1, min=0.0001, max=1000)
        if "K_ic" in params:
            parameters.add("K_ic", value=0.1, min=0.0001, max=1000)
        if self.enzyme_inactivation:
            parameters.add("K_ie", value=0.01, min=0.0001, max=0.9999)
        if "k_inact" in params:
            parameters.add("k_inact", value=0.01, min=0.0001, max=0.9999)

        return parameters
    
    def integrate(self, parameters: Parameters, time: np.ndarray, y0s: List[tuple]) -> np.ndarray:
        """Integrates model based on parameters for a given time array and initial conditions.

        Args:
            parameters (Parameters): lmfit parameters
            time (ndarray): time array for integration
            y0s (List[tuple]): initial conditions for the species in the model.

        Returns:
            result (np.ndarray): integrated model over given time.
        """
        result = [odeint(func=self.model, y0=y0, t=t, args=(parameters, self.enzyme_inactivation)) for y0, t in zip(y0s, time)]
        return np.array(result)
    
    def residuals(self, parameters: Parameters, time: np.ndarray, y0s: List[tuple], ydata: np.ndarray) -> np.ndarray:
        """Calculates residuals between integrated model and measured data (substrate).

        Args:
            parameters (Parameters): lmfit parameters
            time (ndarray): time array, corresponding to ydata.
            y0s (List[tuple]): initial conditions of modeled species
            ydata (ndarray): measured substrate data, corresponding to time data.

        Returns:
            residuals (np.ndarray): List of substrate residuals.
        """

        y0s = np.array(y0s)
        model = self.integrate(parameters, time, y0s)
        residuals = model[:,:,0] - ydata
        return residuals.flatten()
    
    def fit(self, ydata: np.ndarray, time: np.ndarray, y0s: List[tuple]) -> MinimizerResult:
        """Fit model to substrate data.

        Args:
            ydata (ndarray): Experimental substrate data
            time (ndarray): Time array corresponding to measurement data
            y0s (List[tuple]): intial conditions of modeled species

        Returns:
            MinimizerResult: Lest-squares minimization result.
        """
        result = minimize(self.residuals, self.parameters, args=(time, y0s, ydata))
        self._fit_result = result

        return result

