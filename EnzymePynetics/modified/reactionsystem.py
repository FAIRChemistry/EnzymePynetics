import sdRDM

from typing import Optional
from pydantic import Field
from sdRDM.base.utils import forge_signature, IDGenerator


from .modelresult import ModelResult
from .reaction import Reaction


@forge_signature
class ReactionSystem(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("reactionsystemINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the reaction system",
    )

    reactions: Optional[Reaction] = Field(
        default=None,
    )

    result: Optional[ModelResult] = Field(
        default=ModelResult(),
        description="Result of the kinetic model fitting.",
    )
