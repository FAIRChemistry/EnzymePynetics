import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .concentrationtypes import ConcentrationTypes
from .measurement import Measurement
from .series import Series
from .timetypes import TimeTypes
from .kineticmodel import KineticModel
from .parameter import Parameter


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
        default="3de8cc7f43153d5cbb0cbfd736e91aca3ea2eab1"
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
        initial_substrate_conc: float,
        data: List[Series],
        data_conc_unit: ConcentrationTypes,
        time: List[float],
        time_unit: TimeTypes,
        enzyme_conc: Optional[float] = None,
        inhibitor_conc: Optional[float] = None,
        inhibitor_conc_unit: Optional[ConcentrationTypes] = None,
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
        pH: Optional[float] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Measurement' to the attribute 'measurements'.

        Args:


            id (str): Unique identifier of the 'Measurement' object. Defaults to 'None'.


            initial_substrate_conc (float): Initial substrate concentration of the measurement.


            data (List[Series]): One or multiple time-course concentration data arrays.


            data_conc_unit (ConcentrationTypes): Molar concentration unit of the measured data.


            time (List[float]): Time array corresponding to time-course data.


            time_unit (TimeTypes): Time data unit.


            enzyme_conc (Optional[float]): Enzyme concentration in the reaction. Defaults to None


            inhibitor_conc (Optional[float]): Inhibitor concentration, if applied to the reaction. Defaults to None


            inhibitor_conc_unit (Optional[ConcentrationTypes]): Inhibitor concentration in the reaction, if applied. Defaults to None


            temperature (Optional[float]): Temperature of the reaction. Defaults to None


            temperature_unit (Optional[str]): Temperature unit. Defaults to None


            pH (Optional[float]): pH of the reaction. Defaults to None
        """

        params = {
            "initial_substrate_conc": initial_substrate_conc,
            "data": data,
            "data_conc_unit": data_conc_unit,
            "time": time,
            "time_unit": time_unit,
            "enzyme_conc": enzyme_conc,
            "inhibitor_conc": inhibitor_conc,
            "inhibitor_conc_unit": inhibitor_conc_unit,
            "temperature": temperature,
            "temperature_unit": temperature_unit,
            "pH": pH,
        }
        if id is not None:
            params["id"] = id
        measurements = [Measurement(**params)]
        self.measurements = self.measurements + measurements
