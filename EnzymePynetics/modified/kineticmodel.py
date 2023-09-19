import sdRDM
import sympy as sp

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .sboterm import SBOTerm
from .kineticparameter import KineticParameter


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

    def _set_bounds():
        pass

    @property
    def function(self) -> callable:
        rate_expression = self.equation.split("=")[1]
        fun = sp.parse_expr(rate_expression, {"product": sp.Symbol("product")})

        return sp.lambdify(list(fun.free_symbols), fun)

    @property
    def eq_species(self):
        pass

    @property
    def eq_parameters(self):
        pass
