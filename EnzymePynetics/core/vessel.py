import sdRDM

from typing import Optional
from pydantic import Field, PositiveFloat, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Vessel(sdRDM.DataModel):
    """This object describes vessels in which the experiment has been carried out. These can include any type of vessel used in biocatalytic experiments."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("vesselINDEX"),
        xml="@id",
    )

    name: str = Field(
        ...,
        description="Name of the used vessel.",
        template_alias="Name",
    )

    volume: PositiveFloat = Field(
        ...,
        description="Volumetric value of the vessel.",
        template_alias="Volume value",
    )

    unit: str = Field(
        ...,
        description="Volumetric unit of the vessel.",
        template_alias="Volume unit",
    )

    constant: bool = Field(
        description="Whether the volume of the vessel is constant or not.",
        default=True,
    )

    uri: Optional[str] = Field(
        default=None,
        description="URI of the vessel.",
    )

    creator_id: Optional[str] = Field(
        default=None,
        description="Unique identifier of the author.",
    )
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="70285185b8d9c7baf61e12dd52d943624695a510"
    )
