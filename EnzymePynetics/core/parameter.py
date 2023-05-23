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
        description="1 sigma standard deviation of the kinetic parameter.",
    )

    upper_limit: Optional[float] = Field(
        default=None,
        description="Upper limit for parameter value.",
    )

    lower_limit: Optional[float] = Field(
        default=None,
        description="lower limit for parameter value.",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="27257696b2eb57a59c4f86615479b4cf40300291"
    )
