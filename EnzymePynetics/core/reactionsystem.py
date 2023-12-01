import sdRDM

import numpy as np
from typing import Callable, List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from lmfit import Parameters, minimize
from lmfit.minimizer import MinimizerResult
from scipy.integrate import odeint
from .reaction import Reaction
from .sboterm import SBOTerm
from .reactionelement import ReactionElement
from .modelresult import ModelResult
from .correlation import Correlation
from .kineticmodel import KineticModel
from .kineticparameter import KineticParameter
from .parameter import Parameter
from .paramtype import ParamType


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
        description="Result of the kinetic model fitting.",
        default_factory=ModelResult,
    )
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="e6ee5d208d59c44e46c402bcceafa385ca48b435"
    )

    def add_to_reactions(
        self,
        name: str,
        reversible: bool = False,
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
        ph: Optional[float] = None,
        ontology: SBOTerm = SBOTerm.BIOCHEMICAL_REACTION,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        model: Optional[KineticModel] = None,
        educts: List[ReactionElement] = ListPlus(),
        products: List[ReactionElement] = ListPlus(),
        modifiers: List[ReactionElement] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Reaction' to attribute reactions

        Args:
            id (str): Unique identifier of the 'Reaction' object. Defaults to 'None'.
            name (): Name of the reaction..
            reversible (): Whether the reaction is reversible or irreversible. Defaults to False
            temperature (): Numeric value of the temperature of the reaction.. Defaults to None
            temperature_unit (): Unit of the temperature of the reaction.. Defaults to None
            ph (): PH value of the reaction.. Defaults to None
            ontology (): Ontology defining the role of the given species.. Defaults to SBOTerm.BIOCHEMICAL_REACTION
            uri (): URI of the reaction.. Defaults to None
            creator_id (): Unique identifier of the author.. Defaults to None
            model (): Kinetic model decribing the reaction.. Defaults to None
            educts (): List of educts containing ReactionElement objects.. Defaults to ListPlus()
            products (): List of products containing ReactionElement objects.. Defaults to ListPlus()
            modifiers (): List of modifiers (Proteins, snhibitors, stimulators) containing ReactionElement objects.. Defaults to ListPlus()
        """
        params = {
            "name": name,
            "reversible": reversible,
            "temperature": temperature,
            "temperature_unit": temperature_unit,
            "ph": ph,
            "ontology": ontology,
            "uri": uri,
            "creator_id": creator_id,
            "model": model,
            "educts": educts,
            "products": products,
            "modifiers": modifiers,
        }
        if id is not None:
            params["id"] = id
        self.reactions.append(Reaction(**params))
        return self.reactions[-1]

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

    def _create_lmfit_params(self, fixed_params: List[str] = []) -> Parameters:
        parameters = Parameters()

        for reaction in self.reactions:
            for param in reaction.model.parameters:
                if param.name in fixed_params:
                    vary = False
                    value = param.value
                else:
                    vary = True
                    value = param.initial_value
                parameters.add(
                    name=param.name,
                    value=value,
                    min=param.lower,
                    max=param.upper,
                    vary=vary,
                )

        return parameters

    def _setup_ode_model(self) -> Callable:
        substrate_eq = self.substrate.model.function

        if self.enzyme:
            enzyme_eq = self.enzyme.model.function
        else:

            def enzyme_eq(catalyst):  # change to enzyme
                return 0

        def ode_model(species, time: np.ndarray, params: Parameters):
            species_dict = dict(zip(["substrate", "catalyst", "product"], species))

            try:
                params_dict = params.valuesdict()
            except AttributeError:
                params_dict = params

            combined_dict = species_dict | params_dict

            subtrate_args = substrate_eq.__code__.co_varnames
            substrate_dict = {k: combined_dict[k] for k in subtrate_args}

            enzyme_args = enzyme_eq.__code__.co_varnames
            # enzyme_args = ("catalyst",)
            enzyme_dict = {k: combined_dict[k] for k in enzyme_args}

            d_substrate = substrate_eq(**substrate_dict)
            d_enzyme = enzyme_eq(**enzyme_dict)

            return np.array([d_substrate, d_enzyme, -d_substrate])

        return ode_model

    def simulate(
        self, times: np.ndarray, init_conditions: np.ndarray, params: Parameters
    ):
        return np.array([
            odeint(
                func=self._setup_ode_model(),
                y0=init_condition,
                t=time,
                args=(params,),
            )
            for init_condition, time in zip(init_conditions, times)
        ])

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
        return np.array([substrate_data[:, 0], enzyme_data[:, 0], product_data[:, 0]]).T

    def fit(
        self,
        substrate_data: np.ndarray,
        enzyme_data: np.ndarray,
        product_data: np.ndarray,
        times: np.ndarray,
        fixed_params: List[str] = [],
    ):
        params = self._create_lmfit_params(fixed_params=fixed_params)

        init_conditions = self._get_init_conditions(
            substrate_data=substrate_data,
            enzyme_data=enzyme_data,
            product_data=product_data,
        )

        lmfit_result = minimize(
            self.residuals,
            params,
            args=(times, init_conditions, substrate_data),
            method="leastsq",
            max_nfev=300,
        )

        self._update_param_values(lmfit_result)
        self._update_fit_statistics(lmfit_result)

        return lmfit_result

    def get_parameter(self, param_name: str) -> KineticParameter:
        for reaction in self.reactions:
            return reaction.model.get_parameter(param_name)

        raise ValueError(f"Parameter '{param_name}' not found in kinetic model.")

    def _update_param_values(self, result: MinimizerResult):
        if not result.success:
            return

        for reaction in self.reactions:
            for param in reaction.model.parameters:
                param.value = result.params[param.name].value
                param.stdev = result.params[param.name].stderr

    def _update_fit_statistics(self, result: MinimizerResult):
        self.result.fit_success = result.success
        if not result.success:
            return

        self.result.AIC = result.aic
        self.result.BIC = result.bic
        self.result.RMSD = None

        for param in result.params.values():
            if param.correl:
                self.result.parameters.append(
                    Parameter(
                        name=param.name,
                        correlations=[
                            Correlation(parameter_name=key, value=value)
                            for key, value in param.correl.items()
                        ],
                    )
                )

    @property
    def fitted_params_dict(self):
        params = {}
        for reaction in self.reactions:
            for param in reaction.model.parameters:
                params[param.name] = param.value

        return params

    def _style_parameters(self):
        param_name_map = {
            ParamType.K_CAT.value: "<b><i>k</i><sub>cat</sub>:</b>",
            ParamType.K_M.value: "<b><i>K</i><sub>M</sub>:</b>",
            ParamType.K_IE.value: "<b><i>k</i><sub>ie</sub>:</b>",
            ParamType.K_IC.value: "<b><i>K</i><sub>ic</sub>:</b>",
            ParamType.K_IU.value: "<b><i>K</i><sub>iu</sub>:</b>",
        }
        params = ""
        for reaction in self.reactions:
            for parameter in reaction.model.parameters:
                params = (
                    params
                    + f"{param_name_map[parameter.name]} "
                    + f"{parameter.value:.3f} "
                    + f"{self._format_unit(parameter.unit)} \n\n"
                )
        return params

    def get_correlation(self, param_1: str, param_2: str):
        for parameter in self.result.parameters:
            if parameter.name == param_1:
                for correlation in parameter.correlations:
                    if correlation.parameter_name == param_2:
                        return correlation.value

        raise ValueError(f"No correlation found between {param_1} and {param_2}")

    @staticmethod
    def _format_unit(unit: str) -> str:
        unit = unit.replace(" / l", " L<sup>-1</sup>")
        unit = unit.replace("1 / s", "s<sup>-1</sup>")
        unit = unit.replace("1 / min", "min<sup>-1</sup>")
        unit = unit.replace("umol", "µmol")
        unit = unit.replace("ug", "µg")
        return unit
