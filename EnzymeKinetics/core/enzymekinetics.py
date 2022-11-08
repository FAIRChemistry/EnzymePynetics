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
from .stoichiometrytypes import StoichiometryTypes
from .inhibitor import Inhibitor
from .timetypes import TimeTypes


@forge_signature
class EnzymeKinetics(sdRDM.DataModel):
    """Base class, dealing with measurement data of an enzyme kinetics assay."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("enzymekineticsINDEX"),
        xml="@id",
    )

    data_conc_unit: ConcentrationTypes = Field(
        ..., description="Molar concentration unit of the measured data."
    )

    time_unit: TimeTypes = Field(..., description="Time data unit.")

    title: Optional[str] = Field(
        description="Title of the kinetic experiment", default=None
    )

    reactant_name: Optional[str] = Field(
        description="Name of the measured reactant.", default=None
    )

    measurements: List[Measurement] = Field(
        description="Measurement data for a given initial substrate concentration.",
        default_factory=ListPlus,
    )

    stoichiometry: Optional[StoichiometryTypes] = Field(
        description=(
            "Define whether 'substrate' or 'product' concentration was measured."
        ),
        default=None,
    )

    time: List[float] = Field(
        description="Time array corresponding to time-course data.",
        default_factory=ListPlus,
    )

    inhibitors: List[Inhibitor] = Field(
        description="Inhibitor data, if applied to the reaction.",
        default_factory=ListPlus,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/enzyme-kinetics-datamodel.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="618620591831f969dcda4fa6751bb740594ce6b4"
    )

    def add_to_measurements(
        self,
        initial_substrate_conc: float,
        data: List[Series],
        enzyme_conc: Optional[float] = None,
        inhibitor_conc_unit: Optional[ConcentrationTypes] = None,
        inhibitor: Optional[Inhibitor] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Measurement' to the attribute 'measurements'.

        Args:


            id (str): Unique identifier of the 'Measurement' object. Defaults to 'None'.


            initial_substrate_conc (float): Initial substrate concentration of the measurement.


            data (List[Series]): One or multiple time-course concentration data arrays.


            enzyme_conc (Optional[float]): Enzyme concentration in the reaction. Defaults to None


            inhibitor_conc_unit (Optional[ConcentrationTypes]): Inhibitor concentration in the reaction, if applied. Defaults to None


            inhibitor (Optional[Inhibitor]): Inhibitor instance, of the respective inhibitor applied to the reaction. Defaults to None
        """

        params = {
            "initial_substrate_conc": initial_substrate_conc,
            "data": data,
            "enzyme_conc": enzyme_conc,
            "inhibitor_conc_unit": inhibitor_conc_unit,
            "inhibitor": inhibitor,
        }
        if id is not None:
            params["id"] = id
        measurements = [Measurement(**params)]
        self.measurements = self.measurements + measurements

    def add_to_inhibitors(
        self,
        name: Optional[str] = None,
        concentration: Optional[float] = None,
        conconcentration_unit: Optional[ConcentrationTypes] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Inhibitor' to the attribute 'inhibitors'.

        Args:


            id (str): Unique identifier of the 'Inhibitor' object. Defaults to 'None'.


            name (Optional[str]): Name or ID of the inhibitor. Defaults to None


            concentration (Optional[float]): Inhibitor concentration in the reaction. Defaults to None


            conconcentration_unit (Optional[ConcentrationTypes]): Name or ID of the inhibitor. Defaults to None
        """

        params = {
            "name": name,
            "concentration": concentration,
            "conconcentration_unit": conconcentration_unit,
        }
        if id is not None:
            params["id"] = id
        inhibitors = [Inhibitor(**params)]
        self.inhibitors = self.inhibitors + inhibitors
