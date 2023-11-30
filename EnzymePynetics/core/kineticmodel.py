import sdRDM

import sympy as sp
from typing import FrozenSet, List, Optional
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from sympy import pprint
from .kineticparameter import KineticParameter
from .sboterm import SBOTerm, ParamType


@forge_signature
class KineticModel(sdRDM.DataModel):
    """This object describes a kinetic model that was derived from the experiment."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("kineticmodelINDEX"),
        xml="@id",
    )

    name: str = Field(
        ...,
        description="Name of the kinetic law.",
    )

    equation: str = Field(
        ...,
        description="Equation for the kinetic law.",
    )

    parameters: List[KineticParameter] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="List of estimated parameters.",
    )

    ontology: Optional[SBOTerm] = Field(
        default=None,
        description="Type of the estimated parameter.",
    )
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="5dcc898a16a04c37e7fd62bb4b0d81bfd9103184"
    )

    def add_to_parameters(
        self,
        name: str,
        value: float,
        unit: str,
        initial_value: Optional[float] = None,
        upper: Optional[float] = None,
        lower: Optional[float] = None,
        is_global: bool = False,
        stdev: Optional[float] = None,
        constant: bool = False,
        ontology: Optional[SBOTerm] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'KineticParameter' to attribute parameters

        Args:
            id (str): Unique identifier of the 'KineticParameter' object. Defaults to 'None'.
            name (): Name of the estimated parameter..
            value (): Numerical value of the estimated parameter..
            unit (): Unit of the estimated parameter..
            initial_value (): Initial value that was used for the parameter estimation.. Defaults to None
            upper (): Upper bound of the estimated parameter.. Defaults to None
            lower (): Lower bound of the estimated parameter.. Defaults to None
            is_global (): Specifies if this parameter is a global parameter.. Defaults to False
            stdev (): Standard deviation of the estimated parameter.. Defaults to None
            constant (): Specifies if this parameter is constant. Defaults to False
            ontology (): Type of the estimated parameter.. Defaults to None
        """
        params = {
            "name": name,
            "value": value,
            "unit": unit,
            "initial_value": initial_value,
            "upper": upper,
            "lower": lower,
            "is_global": is_global,
            "stdev": stdev,
            "constant": constant,
            "ontology": ontology,
        }
        if id is not None:
            params["id"] = id
        self.parameters.append(KineticParameter(**params))
        return self.parameters[-1]

    @validator("equation")
    def check_equation_symbols(cls, v):
        sp_dict = {"product": sp.Symbol("product")}

        symbol_str, rate_law_str = v.split("=")

        symbol = sp.parse_expr(symbol_str, sp_dict)
        if not len(symbol.free_symbols) == 1:
            raise ValueError(f"Species equation must one of {SPECIES} left of '='.")

        rate_law = sp.parse_expr(rate_law_str, sp_dict)

        params = []
        observables = []
        unknowns = []
        for symbol in rate_law.free_symbols:
            if symbol.name in SPECIES:
                observables.append(symbol.name)
            elif symbol.name in [param.value for param in ParamType]:
                params.append(symbol.name)
            else:
                unknowns.append(symbol.name)

        if len(unknowns) > 0:
            raise ValueError(
                f"Equation '{v}' contains unknown symbols: {unknowns}",
                f"Allowed species are: {SPECIES}",
                f"Allowed parameters are: {[param.name for param in ParamType]}",
            )

        return v

    def add_to_parameters(
        self,
        name: str,
        value: float,
        unit: str,
        initial_value: Optional[float] = None,
        upper: Optional[float] = None,
        lower: Optional[float] = None,
        is_global: bool = False,
        stdev: Optional[float] = None,
        constant: bool = False,
        ontology: Optional[SBOTerm] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'KineticParameter' to attribute parameters

        Args:
            id (str): Unique identifier of the 'KineticParameter' object. Defaults to 'None'.
            name (): Name of the estimated parameter..
            value (): Numerical value of the estimated parameter..
            unit (): Unit of the estimated parameter..
            initial_value (): Initial value that was used for the parameter estimation.. Defaults to None
            upper (): Upper bound of the estimated parameter.. Defaults to None
            lower (): Lower bound of the estimated parameter.. Defaults to None
            is_global (): Specifies if this parameter is a global parameter.. Defaults to False
            stdev (): Standard deviation of the estimated parameter.. Defaults to None
            constant (): Specifies if this parameter is constant. Defaults to False
            ontology (): Type of the estimated parameter.. Defaults to None
        """

        params = {
            "name": name,
            "value": value,
            "unit": unit,
            "initial_value": initial_value,
            "upper": upper,
            "lower": lower,
            "is_global": is_global,
            "stdev": stdev,
            "constant": constant,
            "ontology": ontology,
        }

        if id is not None:
            params["id"] = id

        new_parameter = KineticParameter(**params)

        if any([parameter.name == new_parameter.name for parameter in self.parameters]):
            self.parameters = [
                new_parameter if parameter.name == new_parameter.name else parameter
                for parameter in self.parameters
            ]

            return new_parameter

        else:
            self.parameters.append(new_parameter)

            return new_parameter

    def get_parameter(self, param_name: str) -> KineticParameter:
        for parameter in self.parameters:
            if parameter.name == param_name:
                return parameter

        raise ValueError(f"Parameter '{param_name}' not found in kinetic model.")

    @property
    def function(self) -> callable:
        rate_expression = self.equation.split("=")[1]
        fun = sp.parse_expr(rate_expression, {"product": sp.Symbol("product")})

        return sp.lambdify(list(fun.free_symbols), fun)

    @property
    def eq_species(self) -> FrozenSet[str]:
        return set(
            symbol.name
            for symbol in self._equality.free_symbols
            if symbol.name in SPECIES
        )

    @property
    def _equality(self):
        sp_dict = {"product": sp.Symbol("product")}
        return sp.Equality(*[
            sp.parse_expr(side, sp_dict) for side in self.equation.split("=")
        ])

    @property
    def eq_parameters(self):
        return set(
            symbol.name
            for symbol in self._equality.free_symbols
            if symbol.name in [param.value for param in ParamType]
        )

    @property
    def _sp_rate_law(self):
        rate_law = self.equation.split("=")[1]
        return sp.parse_expr(rate_law, {"product": sp.Symbol("product")})

    @property
    def _sp_species_symbol(self) -> sp.Symbol:
        rate_law = self.equation.split("=")[0]
        species_eq = sp.parse_expr(rate_law, {"product": sp.Symbol("product")})

        if not len(species_eq.free_symbols) == 1:
            raise ValueError("Species equation must contain exactly one free symbol.")

        return species_eq

    @property
    def pretty_print(self):
        return pprint(self._equality)
