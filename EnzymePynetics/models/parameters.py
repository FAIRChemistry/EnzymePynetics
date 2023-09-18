from ..modified import KineticParameter, SBOTerm

K_CAT = KineticParameter(
    name="kcat",
    value=0,
    ontology=SBOTerm.K_CAT,
)

K_M = KineticParameter(
    name="Km",
    value=0,
    ontology=SBOTerm.K_M,
)

K_IC = KineticParameter(
    name="Kic",
    value=0
)

K_IU = KineticParameter(
    name="Kiu",
    value=0
)
