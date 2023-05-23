from typing import List
from lmfit import Parameters, minimize
from lmfit.minimizer import MinimizerResult
from typing import Dict, Callable, Tuple
from numpy import ndarray, log, exp
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
        self.result: MinimizerResult = None        

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
        if "t_0" in params:
            parameters.add("t_0", value=-1, min=-1000, max=1000)
        if "k_inact" in params:
            parameters.add("k_inact", value=0.01, min=0.0001, max=0.9999)

        return parameters
    
    def integrate(self, parameters: Parameters, time: ndarray, y0s: List[tuple]) -> np.ndarray:
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
    
    def residuals(self, parameters: Parameters, time: ndarray, y0s: List[tuple], ydata: np.ndarray) -> np.ndarray:
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
    
    def fit(self, ydata: ndarray, time: ndarray, y0s: List[tuple]) -> MinimizerResult:
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



def irreversible_model(w0: tuple, t, params: Parameters, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value

    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * cS / (Km+cS)
    dc_P = -dc_S
    dc_I = 0

    return (dc_S, dc_E, dc_P, dc_I)

### Product inhibition ###

def competitive_product_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_ic = params["K_ic"].value

    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * cS / (Km*(1+(cI / K_ic))+cS)
    dc_P = -dc_S
    dc_I = dc_P

    return (dc_S, dc_E, dc_P, dc_I)

def uncompetitive_product_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_iu = params["K_iu"].value

    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * cS / (cS*(1+(cI / K_iu))+Km)
    dc_P = -dc_S
    dc_I = dc_P

    return (dc_S, dc_E, dc_P, dc_I)

def noncompetitive_product_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_iu = params["K_iu"].value
    K_ic = params["K_ic"].value

    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * cS / (Km * (1+(cI/K_ic)) + (1+(cI/K_iu)) * cS)
    dc_P = -dc_S
    dc_I = dc_P

    return (dc_S, dc_E, dc_P, dc_I)

### Substrate inhibition ###

def substrate_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_iu = params["K_iu"].value

    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * cS / (Km + ((1+(cS/K_iu))*cS))
    dc_P = -dc_S
    dc_I = dc_S

    return (dc_S, dc_E, dc_P, dc_I)

### External inhibitor models ###

def competitive_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_ic = params["K_ic"].value

    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * cS / (Km*(1+(cI / K_ic))+cS)
    dc_P = -dc_S
    dc_I = 0

    return (dc_S, dc_E, dc_P, dc_I)


def uncompetitive_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_iu = params["K_iu"].value
    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * (cS) / (cS*(1+(cI / K_iu))+Km)
    dc_P = -dc_S
    dc_I = 0

    return (dc_S, dc_E, dc_P, dc_I)

def noncompetitive_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_iu = params["K_iu"].value
    K_ic = params["K_ic"].value
    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * (cS) / (Km * (1+(cI/K_ic)) + (1+(cI/K_iu)) * cS)
    dc_P = -dc_S
    dc_I = 0

    return (dc_S, dc_E, dc_P, dc_I)

def partially_competitive_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_iu = params["K_iu"].value
    K_ic = params["K_ic"].value
    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * (cS) / (Km * ((1+(cI/K_ic)) / (1+(cI/K_iu))) + cS)
    dc_P = -dc_S
    dc_I = 0

    return (dc_S, dc_E, dc_P, dc_I)

## Integrated Models ##

def integrated_MM_model(cS,
                  cE,
                  cS0,
                  params: Parameters,
                  k_inactivation: float = None,
                  enzyme_inactivation: bool = False,
                  ) -> list:
    
    params = params.valuesdict()
    K_m = params["K_m"]
    k_cat = params["k_cat"]
    t_0 = params["t_0"]


    if enzyme_inactivation:
        enzyme = exp(-k_inactivation * enzyme)
    
    return -1/(k_cat*cE)*(K_m* log(cS/cS0) + (cS-cS0)) + t_0