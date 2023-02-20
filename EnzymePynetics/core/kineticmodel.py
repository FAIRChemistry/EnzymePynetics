import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .parameter import Parameter


@forge_signature
class KineticModel(sdRDM.DataModel):
    """Description of a kinetic model"""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("kineticmodelINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(description="Name of the kinetic model.", default=None)

    equation: Optional[str] = Field(
        description="Equation of the kinetic model.", default=None
    )

    parameters: List[Parameter] = Field(
        description="Kinetic parameters of the model.", default_factory=ListPlus
    )

    AIC: Optional[float] = Field(
        description="Akaike information criterion.", default=None
    )

    BIC: Optional[float] = Field(
        description="Bayesian information criterion.", default=None
    )

    RMSD: Optional[float] = Field(
        description="Root mean square deviation between model and measurement data.",
        default=None,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="65530220022f81dc42567f7a1e75530dfdf77be4"
    )

    def add_to_parameters(
        self,
        name: Optional[str] = None,
        value: Optional[float] = None,
        standard_deviation: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Parameter' to the attribute 'parameters'.

        Args:


            id (str): Unique identifier of the 'Parameter' object. Defaults to 'None'.


            name (Optional[str]): Name of the kinetic parameter. Defaults to None


            value (Optional[float]): Value of the kinetic parameter. Defaults to None


            standard_deviation (Optional[float]): Standard deviation of the kinetic parameter. Defaults to None
        """

        params = {
            "name": name,
            "value": value,
            "standard_deviation": standard_deviation,
        }
        if id is not None:
            params["id"] = id
        parameters = [Parameter(**params)]
        self.parameters = self.parameters + parameters
