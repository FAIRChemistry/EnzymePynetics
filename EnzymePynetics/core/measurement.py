import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .timetypes import TimeTypes
from .species import Species


@forge_signature
class Measurement(sdRDM.DataModel):
    """A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("measurementINDEX"),
        xml="@id",
    )

    enzyme_conc: Optional[float] = Field(
        description="Enzyme concentration in the reaction.", default=None
    )

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

    time_unit: Optional[TimeTypes] = Field(description="Time data unit.", default=None)

    species: Optional[Species] = Field(
        description="Reactants of the reaction.", default=None
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="577e0e8515c62e37c47732400090bb756ba93616"
    )
