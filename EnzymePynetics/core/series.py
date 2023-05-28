import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .timetypes import TimeTypes
from .concentrationtypes import ConcentrationTypes


@forge_signature
class Series(sdRDM.DataModel):

    """Time-course data of an individual reaction."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("seriesINDEX"),
        xml="@id",
    )

    values: List[float] = Field(
        description="Time-course data of an individual reaction.",
        default_factory=ListPlus,
        multiple=True,
    )

    value_unit: Optional[ConcentrationTypes] = Field(
        default=None,
        description="Unit of the measurement data.",
    )

    time: List[float] = Field(
        description="Time array corresponding to time-course data.",
        default_factory=ListPlus,
        multiple=True,
    )

    time_unit: Optional[TimeTypes] = Field(
        default=None,
        description="Time data unit.",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="fadd40dc58d78832f5ebbb3627bf1c09494e86ca"
    )
