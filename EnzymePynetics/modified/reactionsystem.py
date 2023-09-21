from math import nan
import sdRDM
from lmfit import Parameters, minimize
from lmfit.minimizer import MinimizerResult
import numpy as np
from scipy.integrate import odeint

from typing import Callable, Optional, List
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator
from sdRDM.base.listplus import ListPlus


from .modelresult import ModelResult
from .reaction import Reaction
from .sboterm import SBOTerm


@forge_signature
class ReactionSystem(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("reactionsystemINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the reaction system",
    )

    reactions: List[Reaction] = Field(
        description="Reactions of the reaction system",
        default_factory=ListPlus,
        multiple=True,
    )

    result: Optional[ModelResult] = Field(
        default=ModelResult(),
        description="Result of the kinetic model fitting.",
    )

    @property
    def substrate(self) -> Reaction:
        for reaction in self.reactions:
            if reaction.educts[0].ontology == SBOTerm.SUBSTRATE.value:
                return reaction

        raise ValueError("No substrate found in reaction system")

    @property
    def enzyme(self) -> Reaction:
        for reaction in self.reactions:
            if reaction.educts[0].ontology == SBOTerm.CATALYST.value:
                return reaction

        return None

    def _create_lmfit_params(self) -> Parameters:
        parameters = Parameters()

        for reaction in self.reactions:
            for param in reaction.model.parameters:
                parameters.add(
                    name=param.name,
                    value=param.initial_value,
                    min=param.lower,
                    max=param.upper,
                )

        return parameters

    def fit_ode(substrate_data: np.ndarray, params: Parameters, time: np.ndarray):
        pass

    def _setup_ode_model(self) -> Callable:
        substrate_eq = self.substrate.model.function

        if self.enzyme:
            enzyme_eq = self.enzyme.model.function
        else:

            def enzyme_eq(enzyme):
                return enzyme

        def ode_model(species, time: np.ndarray, params: Parameters):
            species_dict = dict(zip(["substrate", "enzyme", "product"], species))
            params_dict = params.valuesdict()
            combined_dict = species_dict | params_dict

            subtrate_args = substrate_eq.__code__.co_varnames
            substrate_dict = {k: combined_dict[k] for k in subtrate_args}

            enzyme_args = enzyme_eq.__code__.co_varnames
            enzyme_dict = {k: combined_dict[k] for k in enzyme_args}

            d_substrate = substrate_eq(**substrate_dict)
            d_enzyme = enzyme_eq(**enzyme_dict)

            return np.array([d_substrate, d_enzyme, -d_substrate])

        return ode_model

    def simulate(
        self, times: np.ndarray, init_conditions: np.ndarray, params: Parameters
    ):
        return np.array(
            [
                odeint(
                    func=self._setup_ode_model(),
                    y0=init_condition,
                    t=time,
                    args=(params,),
                )
                for init_condition, time in zip(init_conditions, times)
            ]
        )

    def residuals(
        self,
        params: Parameters,
        times: np.ndarray,
        init_conditions: np.ndarray,
        subtrate_data: np.ndarray,
    ) -> np.ndarray:
        model_data = self.simulate(times, init_conditions, params)
        residuals = model_data[:, :, 0] - subtrate_data  # fitting to substrate data

        return residuals.flatten()

    def _get_init_conditions(
        self,
        substrate_data: np.ndarray,
        product_data: np.ndarray,
        enzyme_data: np.ndarray,
    ) -> np.ndarray:
        return np.array([substrate_data[:, 0], product_data[:, 0], enzyme_data[:, 0]]).T

    def fit(
        self,
        substrate_data: np.ndarray,
        product_data: np.ndarray,
        enzyme_data: np.ndarray,
        times: np.ndarray,
    ):
        params = self._create_lmfit_params()

        init_conditions = self._get_init_conditions(
            substrate_data, product_data, enzyme_data
        )

        lmfit_result = minimize(
            self.residuals,
            params,
            args=(times, init_conditions, substrate_data),
            method="leastsq",
            nan_policy="omit",
        )

        self._safe_results(lmfit_result)

        return lmfit_result

    def _safe_results(self, result: MinimizerResult):
        for reaction in self.reactions:
            for param in reaction.model.parameters:
                print(result.params[param.name].value)
                param.value = result.params[param.name].value
                param.stdev = result.params[param.name].stderr
