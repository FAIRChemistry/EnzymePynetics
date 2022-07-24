import sdRDM


from typing import Optional
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from pydantic import Field


class Parameter(sdRDM.DataModel):

    """This is another object used to describe the parameters of given dataset. As a final note, it is important to use the description of an object to its fullest. As you might noticed, the space in between the object definition ```###``` can be freely used to describe what this object is actually about. Ultimately, this gives you the opportunity to ensure users completely understand what the intention and use case of this object is in a readable way."""

    key: str = Field(
        ...,
        description="Name of the parameter",
        xml="@key",
        dataverse="pyDaRUS.Process.method_parameters.name",
    )

    value: float = Field(
        ...,
        description="Respective value of a parameter",
        xml="@value",
        dataverse="pyDaRUS.Process.method_parameters.value",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/JR-1991/sdrdm-template.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="6aa415cb09360648ff18aa0d9e977b017dd79c9e"
    )
