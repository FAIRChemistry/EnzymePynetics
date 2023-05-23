import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .parameter import Parameter


@forge_signature
class ModelResult(sdRDM.DataModel):

    """Description of a kinetic model"""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("modelresultINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the kinetic model.",
    )

    equation: Optional[str] = Field(
        default=None,
        description="Equation of the kinetic model.",
    )

    parameters: List[Parameter] = Field(
        description="Kinetic parameters of the model.",
        default_factory=ListPlus,
        multiple=True,
    )

    AIC: Optional[float] = Field(
        default=None,
        description="Akaike information criterion.",
    )

    BIC: Optional[float] = Field(
        default=None,
        description="Bayesian information criterion.",
    )

    RMSD: Optional[float] = Field(
        default=None,
        description="Root mean square deviation between model and measurement data.",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="d5daf6f3f0aeef1866a8c6e759604dded4038194"
    )

    def add_to_parameters(
        self,
        name: Optional[str] = None,
        value: Optional[float] = None,
        standard_deviation: Optional[float] = None,
        upper_limit: Optional[float] = None,
        lower_limit: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Parameter' to attribute parameters

        Args:
            id (str): Unique identifier of the 'Parameter' object. Defaults to 'None'.
            name (): Name of the kinetic parameter.. Defaults to None
            value (): Value of the kinetic parameter.. Defaults to None
            standard_deviation (): 1 sigma standard deviation of the kinetic parameter.. Defaults to None
            upper_limit (): Upper limit for parameter value.. Defaults to None
            lower_limit (): lower limit for parameter value.. Defaults to None
        """

        params = {
            "name": name,
            "value": value,
            "standard_deviation": standard_deviation,
            "upper_limit": upper_limit,
            "lower_limit": lower_limit,
        }

        if id is not None:
            params["id"] = id

        self.parameters.append(Parameter(**params))