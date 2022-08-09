import sdRDM


from typing import Optional
from pydantic import PrivateAttr
from sdRDM.base.listplus import ListPlus
from pydantic import Field
from typing import Optional


class Author(sdRDM.DataModel):

    """This is another object that represents the author of the dataset. Please note, that the options here contain all required fields but also custom ones. In this example, the ```Dataverse``` option specifies where each field should be mapped, when exported to a Dataverse format. Hence, these options allow you to link your dataset towards any other data model without writing code by yourself."""

    name: str = Field(
        ...,
        description="Full name including given and family name",
        dataverse="pyDaRUS.Citation.author.name",
    )

    affiliation: Optional[str] = Field(
        description="To which organization the author is affiliated to",
        dataverse="pyDaRUS.Citation.author.affiliation",
        default=None,
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/JR-1991/sdrdm-template.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="2e6c6d2488a46547b281b6243038ee0eb29c0b97"
    )
