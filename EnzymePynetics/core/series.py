import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


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

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="f55399675de0cb2f40c1c7d175da4e51d3cdfdd2"
    )
