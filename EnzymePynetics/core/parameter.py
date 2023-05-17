import sdRDM

from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Parameter(sdRDM.DataModel):

    """Defines a kinetic parameter."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("parameterINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the kinetic parameter.",
    )

    value: Optional[float] = Field(
        default=None,
        description="Value of the kinetic parameter.",
    )

    standard_deviation: Optional[float] = Field(
        default=None,
        description="Standard deviation of the kinetic parameter.",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="5dd38e6325dc3e396867b90d4d975126b5f12c6b"
    )
