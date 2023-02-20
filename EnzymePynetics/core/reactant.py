import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field
from typing import Optional

from .reactanttypes import ReactantTypes


@forge_signature
class Reactant(sdRDM.DataModel):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("reactantINDEX"),
        xml="@id",
    )
    name: Optional[str] = Field(
        description="name of the reactant.",
        default=None,
    )

    reactant_type: Optional[ReactantTypes] = Field(
        description=(
            "Define whether 'substrate' or 'product' concentration was measured."
        ),
        default=None,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="3de8cc7f43153d5cbb0cbfd736e91aca3ea2eab1"
    )
