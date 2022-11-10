from lmfit import Parameters
from typing import Dict, Callable, Tuple
from numpy import ndarray

class KineticModel():
    def __init__(self,
                 name: str,
                 model: Callable,
                 params: list,
                 kcat_initial: float,
                 Km_initial: float,
                 w0: Dict[str, ndarray],
                 time: ndarray 
                 ) -> None:

        self.name = name
        self.model = model
        self.params = params
        self.time = time

        self.w0 = w0
        self.kcat_initial = kcat_initial
        self.Km_initial = Km_initial
        self.parameters = self.set_params(params)
        self.w0 = w0
        self.result = None

    def _set_params(self, params:list) -> Parameters:

        parameters = Parameters()

        parameters.add('k_cat', value=self.kcat_initial,
                        min=self.kcat_initial/100, max=self.kcat_initial*100)
        parameters.add('Km', value=self.Km_initial, min=self.Km_initial/100,
                        max=max(self.Km_initial)*1000)
        if "K_ie" in params:
            parameters.add("K_ie", value=0.01, min=0.0001, max=0.9999)
        if "K_iu" in params:
            parameters.add("K_iu", value=0.1, min=0.0001, max=1000)
        if "K_ic" in params:
            parameters.add("K_ic", value=0.1, min=0.0001, max=1000)


def irreversible_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
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


