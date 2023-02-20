from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field
from typing import Optional

from .abstractspecies import AbstractSpecies


@forge_signature
class Inhibitor(AbstractSpecies):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("inhibitorINDEX"),
        xml="@id",
    )
    ergerg: Optional[str] = Field(
        description="ergerge",
        default=None,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="65530220022f81dc42567f7a1e75530dfdf77be4"
    )
