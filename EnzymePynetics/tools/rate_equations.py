from lmfit import Parameters
from numpy import log, exp


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

    return -1/(k_cat*cE)*(K_m * log(cS/cS0) + (cS-cS0)) + t_0


# External inhibitor + substrate inhibition


def competitive_inhibition_with_substrate_inhibition_model(w0: tuple, t, params, flag_enzyme_inactivation: bool) -> tuple:
    cS, cE, cP, cI = w0

    k_cat = params['k_cat'].value
    Km = params['Km'].value
    K_ic = params["K_ic"].value
    K_iu = params["K_iu"].value

    if flag_enzyme_inactivation:
        K_ie = params["K_ie"].value
        dc_E = -K_ie * cE
    else:
        dc_E = 0

    dc_S = -k_cat * cE * cS / (Km*(1+(cI / K_ic)) + ((1+(cS/K_iu))*cS))
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
