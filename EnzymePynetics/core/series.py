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
        default="f8a545552dffd4a2f6848b989b2f6e3061d8e08f"
    )
