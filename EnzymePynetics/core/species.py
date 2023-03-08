import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .concentrationtypes import ConcentrationTypes
from .reactanttypes import ReactantTypes
from .series import Series


@forge_signature
class Species(sdRDM.DataModel):
    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("speciesINDEX"),
        xml="@id",
    )

    name: str = Field(..., description="name of the reactant.")

    conc_unit: ConcentrationTypes = Field(
        ..., description="Concentration unit of the measurement data."
    )

    initial_conc: float = Field(
        ..., description="Initial concentration of the reactant."
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
        default="aeef6068d7875e75bfc81dacd796d9320df104e3"
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
