import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .species import Species
from .measurement import Measurement
from .parameter import Parameter
from .modelresult import ModelResult


@forge_signature
class EnzymeKinetics(sdRDM.DataModel):
    """Base class, dealing with measurement data of an enzyme kinetics assay."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("enzymekineticsINDEX"),
        xml="@id",
    )

    title: Optional[str] = Field(
        default=None,
        description="Title of the kinetic experiment.",
    )

    model_results: List[ModelResult] = Field(
        description="Fitted kinetic models which were used for parameter estimation.",
        default_factory=ListPlus,
        multiple=True,
    )

    measurements: List[Measurement] = Field(
        description="Measurement data for a given initial substrate concentration.",
        default_factory=ListPlus,
        multiple=True,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="0cf682ecfef9dde7d67b9d03ea7fc30b2fca9ee4"
    )

    def add_to_model_results(
        self,
        name: Optional[str] = None,
        equations: List[str] = ListPlus(),
        parameters: List[Parameter] = ListPlus(),
        fit_success: Optional[bool] = None,
        AIC: Optional[float] = None,
        BIC: Optional[float] = None,
        RMSD: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'ModelResult' to attribute model_results

        Args:
            id (str): Unique identifier of the 'ModelResult' object. Defaults to 'None'.
            name (): Name of the kinetic model.. Defaults to None
            equations (): Equation of the kinetic model.. Defaults to ListPlus()
            parameters (): Kinetic parameters of the model.. Defaults to ListPlus()
            fit_success (): Whether or not model fitting was possible.. Defaults to None
            AIC (): Akaike information criterion.. Defaults to None
            BIC (): Bayesian information criterion.. Defaults to None
            RMSD (): Root mean square deviation between model and measurement data.. Defaults to None
        """

        params = {
            "name": name,
            "equations": equations,
            "parameters": parameters,
            "fit_success": fit_success,
            "AIC": AIC,
            "BIC": BIC,
            "RMSD": RMSD,
        }

        if id is not None:
            params["id"] = id

        self.model_results.append(ModelResult(**params))

    def add_to_measurements(
        self,
        species: List[Species] = ListPlus(),
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
        pH: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Measurement' to attribute measurements

        Args:
            id (str): Unique identifier of the 'Measurement' object. Defaults to 'None'.
            species (): Reactants of the reaction.. Defaults to ListPlus()
            temperature (): Temperature of the reaction.. Defaults to None
            temperature_unit (): Temperature unit.. Defaults to None
            pH (): pH of the reaction. Defaults to None
        """

        params = {
            "species": species,
            "temperature": temperature,
            "temperature_unit": temperature_unit,
            "pH": pH,
        }

        if id is not None:
            params["id"] = id

        self.measurements.append(Measurement(**params))
