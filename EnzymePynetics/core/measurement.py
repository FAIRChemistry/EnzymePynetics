import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from .concentrationtypes import ConcentrationTypes
from .series import Series
from .timetypes import TimeTypes
from .inhibitor import Inhibitor
from .reactant import Reactant
from .reactanttypes import ReactantTypes


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

    time_unit: TimeTypes = Field(..., description="Time data unit.")

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

    reactants: List[Reactant] = Field(
        description="Reactants of the reaction.", default_factory=ListPlus
    )

    inhibitor: Optional[Inhibitor] = Field(
        description="Inhibitor applied to the reaction.", default=None
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/EnzymePynetics.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="65530220022f81dc42567f7a1e75530dfdf77be4"
    )

    def add_to_reactants(
        self,
        name: str,
        conc_unit: ConcentrationTypes,
        initial_conc: float,
        data: List[Series],
        reactant_type: Optional[ReactantTypes] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Reactant' to the attribute 'reactants'.

        Args:


            id (str): Unique identifier of the 'Reactant' object. Defaults to 'None'.


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
        reactants = [Reactant(**params)]
        self.reactants = self.reactants + reactants
