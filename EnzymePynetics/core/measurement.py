import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .concentrationtypes import ConcentrationTypes
from .series import Series
from .timetypes import TimeTypes


@forge_signature
class Measurement(sdRDM.DataModel):
    """A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("measurementINDEX"),
        xml="@id",
    )

    initial_substrate_conc: float = Field(
        ..., description="Initial substrate concentration of the measurement."
    )

    enzyme_conc: Optional[float] = Field(
        description="Enzyme concentration in the reaction.", default=None
    )

    data: List[Series] = Field(
        description="One or multiple time-course concentration data arrays.",
        default_factory=ListPlus,
    )

    inhibitor_conc: Optional[float] = Field(
        description="Inhibitor concentration, if applied to the reaction.", default=None
    )

    inhibitor_conc_unit: Optional[ConcentrationTypes] = Field(
        description="Inhibitor concentration in the reaction, if applied.", default=None
    )

    data_conc_unit: ConcentrationTypes = Field(
        ..., description="Molar concentration unit of the measured data."
    )

    time_unit: TimeTypes = Field(..., description="Time data unit.")

    temperature: Optional[float] = Field(
        description="Temperature of the reaction.", default=None
    )

    temperature_unit: Optional[str] = Field(
        description="Temperature unit.", default=None
    )

    pH: Optional[float] = Field(description="pH of the reaction", default=None)

    time: List[float] = Field(
        description="Time array corresponding to time-course data.",
        default_factory=ListPlus,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="3de8cc7f43153d5cbb0cbfd736e91aca3ea2eab1"
    )

    def add_to_data(self, values: List[float], id: Optional[str] = None) -> None:
        """
        Adds an instance of 'Series' to the attribute 'data'.

        Args:


            id (str): Unique identifier of the 'Series' object. Defaults to 'None'.


            values (List[float]): Time-course data of an individual reaction.
        """

        params = {"values": values}
        if id is not None:
            params["id"] = id
        data = [Series(**params)]
        self.data = self.data + data
