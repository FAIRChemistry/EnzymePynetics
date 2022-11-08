import sdRDM

from typing import Optional, Union
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from .concentrationtypes import ConcentrationTypes


@forge_signature
class Inhibitor(sdRDM.DataModel):
    """Parameters of the inhibitor, if inhibitor was applied to the enzyme reaction."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("inhibitorINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(description="Name or ID of the inhibitor", default=None)

    concentration: Optional[float] = Field(
        description="Inhibitor concentration in the reaction.", default=None
    )

    conconcentration_unit: Optional[ConcentrationTypes] = Field(
        description="Name or ID of the inhibitor", default=None
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/enzyme-kinetics-datamodel.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="00e6ba1119bf53a6e603b00423e11858c9469283"
    )
