import sdRDM

from typing import Optional, Union
from typing import List
from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


@forge_signature
class Measurement(sdRDM.DataModel):
    """A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling.
    """

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("measurementINDEX"),
        xml="@id",
    )

    initial_conc: Optional[float] = Field(
        description="Initial substrate concentration of the measurement.", default=None
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/haeussma/enzyme-kinetics-datamodel.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="e2fad0d34e83be28e91495256d6f0dc4f7c278e8"
    )
