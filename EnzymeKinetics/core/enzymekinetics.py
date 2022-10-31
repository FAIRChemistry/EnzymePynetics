import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field
from typing import List
from typing import Optional

from .concentrationtypes import ConcentrationTypes
from .measurement import Measurement
from .series import Series
from .stoichiometrytypes import StoichiometryTypes
from .timetypes import TimeTypes


@forge_signature
class EnzymeKinetics(sdRDM.DataModel):

    """Base class, dealing with measurement data of an enzyme kinetics assay."""

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("enzymekineticsINDEX"),
        xml="@id",
    )
    data_conc_unit: ConcentrationTypes = Field(
        ...,
        description="Molar concentration unit of the measured data.",
    )

    time_unit: TimeTypes = Field(
        ...,
        description="Time data unit.",
    )

    title: Optional[str] = Field(
        description="Title of the kinetic experiment",
        default=None,
    )

    reactant_name: Optional[str] = Field(
        description="Name of the measured reactant.",
        default=None,
    )

    measurements: List[Measurement] = Field(
        description="Measurement data for a given initial substrate concentration.",
        default_factory=ListPlus,
    )

    stoichiometry: Optional[StoichiometryTypes] = Field(
        description=(
            "Define whether 'substrate' or 'product' concentration was measured."
        ),
        default=None,
    )

    time: List[float] = Field(
        description="Time array corresponding to time-course data.",
        default_factory=ListPlus,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/enzyme-kinetics-datamodel.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="e2fad0d34e83be28e91495256d6f0dc4f7c278e8"
    )

    def add_to_measurements(
        self,
        initial_substrate_conc: float,
        enzyme_conc: List[float],
        data: List[Series],
        values: List[float],
        inhibitor_conc: Optional[float] = None,
        inhibitor_conc_unit: Optional[str] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        Adds an instance of 'Measurement' to the attribute 'measurements'.

        Args:
            id (str): Unique identifier of the 'Measurement' object. Defaults to 'None'.
            initial_substrate_conc (float): Initial substrate concentration of the measurement.
            enzyme_conc (List[float]): Enzyme concentration in the measurement.
            data (List[Series]): One or multiple time-course concentration data arrays.
            values (List[float]): Time-course data of an individual reaction.
            inhibitor_conc (Optional[float]): inhibitor concentration for the measurement, if inhibitor was present. Defaults to None
            inhibitor_conc_unit (Optional[str]): Concentration unit of the inhibitior. Defaults to None
        """

        params = {
            "initial_substrate_conc": initial_substrate_conc,
            "enzyme_conc": enzyme_conc,
            "data": data,
            "values": values,
            "inhibitor_conc": inhibitor_conc,
            "inhibitor_conc_unit": inhibitor_conc_unit,
        }

        if id is not None:
            params["id"] = id

        measurements = [Measurement(**params)]

        self.measurements = self.measurements + measurements
