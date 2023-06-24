import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .timetypes import TimeTypes
from .series import Series
from .concentrationtypes import ConcentrationTypes
from .speciestypes import SpeciesTypes


@forge_signature
class Species(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
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

    time_unit: Optional[TimeTypes] = Field(
        default=None,
        description="Time data unit.",
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
        default="df3e2f4bc8efec727aab1644f342fefba9029ee7"
    )

    def add_to_data(
        self,
        values: List[float] = ListPlus(),
        time: List[float] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Series' to attribute data

        Args:
            id (str): Unique identifier of the 'Series' object. Defaults to 'None'.
            values (): Time-course data of an individual reaction.. Defaults to ListPlus()
            time (): Time array corresponding to time-course data.. Defaults to ListPlus()
        """

        params = {
            "values": values,
            "time": time,
        }

        if id is not None:
            params["id"] = id

        self.data.append(Series(**params))
