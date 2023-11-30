import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .correlation import Correlation
from .parameter import Parameter


@forge_signature
class ModelResult(sdRDM.DataModel):
    """Description of a kinetic model"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("modelresultINDEX"),
        xml="@id",
    )

    equations: List[str] = Field(
        description="Equation of the kinetic model.",
        default_factory=ListPlus,
        multiple=True,
    )

    parameters: List[Parameter] = Field(
        description="Kinetic parameters of the model.",
        default_factory=ListPlus,
        multiple=True,
    )

    fit_success: Optional[bool] = Field(
        default=None,
        description="Whether or not model fitting was possible.",
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
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="25b9eaa1ad02d290cf2a98a59b5a6f1730cb7652"
    )

    def add_to_parameters(
        self,
        name: Optional[str] = None,
        correlations: List[Correlation] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Parameter' to attribute parameters

        Args:
            id (str): Unique identifier of the 'Parameter' object. Defaults to 'None'.
            name (): Name of the kinetic parameter.. Defaults to None
            correlations (): Correlation of parameter to other parameters of a model.. Defaults to ListPlus()
        """
        params = {"name": name, "correlations": correlations}
        if id is not None:
            params["id"] = id
        self.parameters.append(Parameter(**params))
        return self.parameters[-1]
