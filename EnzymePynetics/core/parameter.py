import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .correlation import Correlation


@forge_signature
class Parameter(sdRDM.DataModel):

    """Defines a kinetic parameter."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("parameterINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the kinetic parameter.",
    )

    value: Optional[float] = Field(
        default=None,
        description="Value of the kinetic parameter.",
    )

    standard_deviation: Optional[float] = Field(
        default=None,
        description="1 sigma standard deviation of the kinetic parameter.",
    )

    correlations: List[Correlation] = Field(
        default_factory=ListPlus,
        multiple=True,
        descritpion="Correlation of parameter to other parameters of a model.",
    )

    upper_limit: Optional[float] = Field(
        default=None,
        description="Upper limit for parameter value.",
    )

    lower_limit: Optional[float] = Field(
        default=None,
        description="lower limit for parameter value.",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="5df8e2db86f7bbf45cece844377d5affbe7ec235"
    )

    def add_to_correlations(
        self,
        parameter: Optional[str] = None,
        value: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Correlation' to attribute correlations

        Args:
            id (str): Unique identifier of the 'Correlation' object. Defaults to 'None'.
            parameter (): Name of the parameter.. Defaults to None
            value (): Correlation value between -1 and 1.. Defaults to None
        """

        params = {
            "parameter": parameter,
            "value": value,
        }

        if id is not None:
            params["id"] = id

        self.correlations.append(Correlation(**params))
