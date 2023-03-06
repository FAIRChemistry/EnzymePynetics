import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .measurement import Measurement
from .timetypes import TimeTypes
from .kineticmodel import KineticModel
from .parameter import Parameter
from .species import Species


@forge_signature
class EnzymeKinetics(sdRDM.DataModel):
    """Base class, dealing with measurement data of an enzyme kinetics assay."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("enzymekineticsINDEX"),
        xml="@id",
    )

    measurements: List[Measurement] = Field(
        description="Measurement data for a given initial substrate concentration.",
        default_factory=ListPlus,
    )

    title: Optional[str] = Field(
        description="Title of the kinetic experiment.", default=None
    )

    kinetic_models: List[KineticModel] = Field(
        description="Kinetic moodels which were used for parameter estimation.",
        default_factory=ListPlus,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="f8440c4fbdc7c22ae6271ae9d554d1779f238938"
    )

    def add_to_kinetic_models(
        self,
        parameters: List[Parameter],
        name: Optional[str] = None,
        equation: Optional[str] = None,
        AIC: Optional[float] = None,
        BIC: Optional[float] = None,
        RMSD: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'KineticModel' to the attribute 'kinetic_models'.

        Args:


            id (str): Unique identifier of the 'KineticModel' object. Defaults to 'None'.


            parameters (List[Parameter]): Kinetic parameters of the model.


            name (Optional[str]): Name of the kinetic model. Defaults to None


            equation (Optional[str]): Equation of the kinetic model. Defaults to None


            AIC (Optional[float]): Akaike information criterion. Defaults to None


            BIC (Optional[float]): Bayesian information criterion. Defaults to None


            RMSD (Optional[float]): Root mean square deviation between model and measurement data. Defaults to None
        """

        params = {
            "parameters": parameters,
            "name": name,
            "equation": equation,
            "AIC": AIC,
            "BIC": BIC,
            "RMSD": RMSD,
        }
        if id is not None:
            params["id"] = id
        kinetic_models = [KineticModel(**params)]
        self.kinetic_models = self.kinetic_models + kinetic_models

    def add_to_measurements(
        self,
        species: List[Species],
        time: List[float],
        enzyme_conc: Optional[float] = None,
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
        pH: Optional[float] = None,
        time_unit: Optional[TimeTypes] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Measurement' to the attribute 'measurements'.

        Args:


            id (str): Unique identifier of the 'Measurement' object. Defaults to 'None'.


            species (List[Species]): Reactants of the reaction.


            time (List[float]): Time array corresponding to time-course data.


            enzyme_conc (Optional[float]): Enzyme concentration in the reaction. Defaults to None


            temperature (Optional[float]): Temperature of the reaction. Defaults to None


            temperature_unit (Optional[str]): Temperature unit. Defaults to None


            pH (Optional[float]): pH of the reaction. Defaults to None


            time_unit (Optional[TimeTypes]): Time data unit. Defaults to None
        """

        params = {
            "species": species,
            "time": time,
            "enzyme_conc": enzyme_conc,
            "temperature": temperature,
            "temperature_unit": temperature_unit,
            "pH": pH,
            "time_unit": time_unit,
        }
        if id is not None:
            params["id"] = id
        measurements = [Measurement(**params)]
        self.measurements = self.measurements + measurements
