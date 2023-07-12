import sdRDM

from typing import List, Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .speciestypes import SpeciesTypes
from .concentrationtypes import ConcentrationTypes
from .timetypes import TimeTypes
from .series import Series
from .species import Species


@forge_signature
class Measurement(sdRDM.DataModel):
    """A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("measurementINDEX"),
        xml="@id",
    )

    species: List[Species] = Field(
        description="Reactants of the reaction.",
        default_factory=ListPlus,
        multiple=True,
    )

    temperature: Optional[float] = Field(
        default=None,
        description="Temperature of the reaction.",
    )

    temperature_unit: Optional[str] = Field(
        default=None,
        description="Temperature unit.",
    )

    pH: Optional[float] = Field(
        default=None,
        description="pH of the reaction",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="de7db67cfa26a2116c3dfd86376b03ccadf2cacf"
    )

    def add_to_species(
        self,
        name: Optional[str] = None,
        conc_unit: Optional[ConcentrationTypes] = None,
        time_unit: Optional[TimeTypes] = None,
        initial_conc: Optional[float] = None,
        species_type: Optional[SpeciesTypes] = None,
        data: List[Series] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Species' to attribute species

        Args:
            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.
            name (): name of the reactant.. Defaults to None
            conc_unit (): Concentration unit of the measurement data.. Defaults to None
            time_unit (): Time data unit.. Defaults to None
            initial_conc (): Initial concentration of the reactant.. Defaults to None
            species_type (): Define the role of the species in the reaction.. Defaults to None
            data (): One or multiple time-course measurement data arrays.. Defaults to ListPlus()
        """

        params = {
            "name": name,
            "conc_unit": conc_unit,
            "time_unit": time_unit,
            "initial_conc": initial_conc,
            "species_type": species_type,
            "data": data,
        }

        if id is not None:
            params["id"] = id

        self.species.append(Species(**params))
