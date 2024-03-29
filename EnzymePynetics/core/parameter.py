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

    correlations: List[Correlation] = Field(
        default_factory=ListPlus,
        multiple=True,
        description="Correlation of parameter to other parameters of a model.",
    )
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="70285185b8d9c7baf61e12dd52d943624695a510"
    )

    def add_to_correlations(
        self,
        parameter_name: Optional[str] = None,
        value: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Correlation' to attribute correlations

        Args:
            id (str): Unique identifier of the 'Correlation' object. Defaults to 'None'.
            parameter_name (): Name of the parameter.. Defaults to None
            value (): Correlation value between -1 and 1.. Defaults to None
        """
        params = {"parameter_name": parameter_name, "value": value}
        if id is not None:
            params["id"] = id
        self.correlations.append(Correlation(**params))
        return self.correlations[-1]
