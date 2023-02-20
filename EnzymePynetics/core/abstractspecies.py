import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field

from .concentrationtypes import ConcentrationTypes


@forge_signature
class AbstractSpecies(sdRDM.DataModel):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("abstractspeciesINDEX"),
        xml="@id",
    )
    name: str = Field(
        ...,
        description="name of the reactant.",
    )

    conc_unit: ConcentrationTypes = Field(
        ...,
        description="Concentration unit of the measurement data.",
    )

    initial_conc: float = Field(
        ...,
        description="Initial concentration of the reactant.",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="65530220022f81dc42567f7a1e75530dfdf77be4"
    )
