import sdRDM

from typing import Optional, Union
from typing import List
from pydantic import PrivateAttr
from pydantic import Field
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
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="bc9021f3deca4945ad07d6e331ee26510a76f7ca"
    )
