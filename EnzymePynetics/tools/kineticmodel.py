from lmfit import Parameters
from lmfit.minimizer import MinimizerResult
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
                 enzyme_inactivation: bool = False
                 ) -> None:

        self.name = name
        self.model = model
        self.params = params
        self.enzyme_inactivation = enzyme_inactivation
        self.w0 = w0
        self.kcat_initial = kcat_initial
        self.Km_initial = Km_initial
        self.parameters = self._set_parameters(params)
        self.result: MinimizerResult = None

    def _set_parameters(self, params:list) -> Parameters:

        parameters = Parameters()

        parameters.add('k_cat', value=self.kcat_initial,
                        min=self.kcat_initial/100, max=self.kcat_initial*100)
        parameters.add('Km', value=self.Km_initial*100, min=self.Km_initial/100,
                        max=self.Km_initial*10000)
        if "K_iu" in params:
            parameters.add("K_iu", value=0.1, min=0.0001, max=1000)
        if "K_ic" in params:
            parameters.add("K_ic", value=0.1, min=0.0001, max=1000)
        if self.enzyme_inactivation:
            parameters.add("K_ie", value=0.01, min=0.0001, max=0.9999)

        return parameters

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
