import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .series import Series
from .concentrationtypes import ConcentrationTypes
from .speciestypes import SpeciesTypes


@forge_signature
class Species(sdRDM.DataModel):

    """"""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("speciesINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="name of the reactant.",
    )

    conc_unit: Optional[ConcentrationTypes] = Field(
        default=None,
        description="Concentration unit of the measurement data.",
    )

    initial_conc: Optional[float] = Field(
        default=None,
        description="Initial concentration of the reactant.",
    )

    species_type: Optional[SpeciesTypes] = Field(
        default=None,
        description="Define the role of the species in the reaction.",
    )

    data: List[Series] = Field(
        description="One or multiple time-course measurement data arrays.",
        default_factory=ListPlus,
        multiple=True,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="e6b0d4e43f04fb348415726dbcaa82e8ea0c63ad"
    )

    def add_to_data(
        self, values: List[float] = ListPlus(), id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'Series' to attribute data

        Args:
            id (str): Unique identifier of the 'Series' object. Defaults to 'None'.
            values (): Time-course data of an individual reaction.. Defaults to ListPlus()
        """

        params = {
            "values": values,
        }

        if id is not None:
            params["id"] = id

        self.data.append(Series(**params))
