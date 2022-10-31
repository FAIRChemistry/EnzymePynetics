import sdRDM

from typing import Optional, Union
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import Field
from typing import List
from typing import Optional

from .series import Series


@forge_signature
class Measurement(sdRDM.DataModel):

    """A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling.
    """

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("measurementINDEX"),
        xml="@id",
    )
    initial_substrate_conc: float = Field(
        ...,
        description="Initial substrate concentration of the measurement.",
    )

    enzyme_conc: List[float] = Field(
        description="Enzyme concentration in the measurement.",
        default_factory=ListPlus,
    )

    inhibitor_conc: Optional[float] = Field(
        description=(
            "inhibitor concentration for the measurement, if inhibitor was present."
        ),
        default=None,
    )

    inhibitor_conc_unit: Optional[str] = Field(
        description="Concentration unit of the inhibitior.",
        default=None,
    )

    data: List[Series] = Field(
        description="One or multiple time-course concentration data arrays",
        default_factory=ListPlus,
    )

    values: List[float] = Field(
        description="Time-course data of an individual reaction.",
        default_factory=ListPlus,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/enzyme-kinetics-datamodel.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="1ded8d01b49366f4dd52e8773b83ded9b364db93"
    )
