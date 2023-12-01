import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Correlation(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("correlationINDEX"),
        xml="@id",
    )

    parameter_name: Optional[str] = Field(
        default=None,
        description="Name of the parameter.",
    )

    value: Optional[float] = Field(
        default=None,
        description="Correlation value between -1 and 1.",
    )
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="32f5c786a31e457bebca5c7812c1c6df2c1a9857"
    )
