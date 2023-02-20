import sdRDM

from typing import Optional, Union
from typing import Optional
from typing import List
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .reactanttypes import ReactantTypes
from .abstractspecies import AbstractSpecies
from .series import Series


@forge_signature
class Reactant(sdRDM.DataModel):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("reactantINDEX"),
        xml="@id",
    )

    reactant_type: Optional[ReactantTypes] = Field(
        description=(
            "Define whether 'substrate' or 'product' concentration was measured."
        ),
        default=None,
    )

    data: List[Series] = Field(
        description="One or multiple time-course measurement data arrays.",
        default_factory=ListPlus,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="65530220022f81dc42567f7a1e75530dfdf77be4"
    )

    def add_to_data(self, values: List[float], id: Optional[str] = None) -> None:
        """
        Adds an instance of 'Series' to the attribute 'data'.

        Args:


            id (str): Unique identifier of the 'Series' object. Defaults to 'None'.


            values (List[float]): Time-course data of an individual reaction.
        """

        params = {"values": values}
        if id is not None:
            params["id"] = id
        data = [Series(**params)]
        self.data = self.data + data
