from EnzymePynetics.modified.kineticparameter import KineticParameter
from EnzymePynetics.modified.sboterm import SBOTerm

K_CAT = KineticParameter(
    name="k_cat",
    value=0,
    unit="",
    ontology=SBOTerm.K_CAT,
)

K_M = KineticParameter(
    name="K_m",
    value=0,
    unit="",
    ontology=SBOTerm.K_M,
)

K_IC = KineticParameter(name="K_ic", value=0, unit="")

K_IU = KineticParameter(name="K_iu", value=0, unit="")

K_IE = KineticParameter(name="k_ie", value=0, unit="")

PARAMS = [K_CAT, K_M, K_IC, K_IU, K_IE]
