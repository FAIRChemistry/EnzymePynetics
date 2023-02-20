import sdRDM

from typing import Optional, Union
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Parameter(sdRDM.DataModel):
    """Defines a kinetic parameter"""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("parameterINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        description="Name of the kinetic parameter", default=None
    )

    value: Optional[float] = Field(
        description="Value of the kinetic parameter.", default=None
    )

    standard_deviation: Optional[float] = Field(
        description="Standard deviation of the kinetic parameter.", default=None
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="3de8cc7f43153d5cbb0cbfd736e91aca3ea2eab1"
    )
