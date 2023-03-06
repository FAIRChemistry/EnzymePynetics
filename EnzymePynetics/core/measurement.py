import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .timetypes import TimeTypes
from .species import Species
from .concentrationtypes import ConcentrationTypes
from .reactanttypes import ReactantTypes
from .series import Series


@forge_signature
class Measurement(sdRDM.DataModel):
    """A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("measurementINDEX"),
        xml="@id",
    )

    enzyme_conc: Optional[float] = Field(
        description="Enzyme concentration in the reaction.", default=None
    )

    temperature: Optional[float] = Field(
        description="Temperature of the reaction.", default=None
    )

    temperature_unit: Optional[str] = Field(
        description="Temperature unit.", default=None
    )

    pH: Optional[float] = Field(description="pH of the reaction", default=None)

    time: List[float] = Field(
        description="Time array corresponding to time-course data.",
        default_factory=ListPlus,
    )

    time_unit: Optional[TimeTypes] = Field(description="Time data unit.", default=None)

    species: List[Species] = Field(
        description="Reactants of the reaction.", default_factory=ListPlus
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="f8440c4fbdc7c22ae6271ae9d554d1779f238938"
    )

    def add_to_species(
        self,
        name: str,
        conc_unit: ConcentrationTypes,
        initial_conc: float,
        data: List[Series],
        reactant_type: Optional[ReactantTypes] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Species' to the attribute 'species'.

        Args:


            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.


            name (str): name of the reactant.


            conc_unit (ConcentrationTypes): Concentration unit of the measurement data.


            initial_conc (float): Initial concentration of the reactant.


            data (List[Series]): One or multiple time-course measurement data arrays.


            reactant_type (Optional[ReactantTypes]): Define whether "substrate" or "product" concentration was measured. Defaults to None
        """

        params = {
            "name": name,
            "conc_unit": conc_unit,
            "initial_conc": initial_conc,
            "data": data,
            "reactant_type": reactant_type,
        }
        if id is not None:
            params["id"] = id
        species = [Species(**params)]
        self.species = self.species + species
