import sdRDM

from typing import Optional
from typing import List
from typing import Optional, Union
from pydantic import PrivateAttr
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from .author import Author
from .parameter import Parameter


class Root(sdRDM.DataModel):
    """This is the root of the data model and contains all objects defined in this example. While its good practice to have a single root, you can define as many roots as you like. Furthermore, the name does not have to be ```Root``` and can be any other name.
    """

    description: str = Field(
        ...,
        description="Describes the content of the dataset.",
        dataverse="pyDaRUS.Citation.description.text",
    )

    title: str = Field(
        ..., description="Title of the work", dataverse="pyDaRUS.Citation.title"
    )

    subject: str = Field(
        ...,
        description="Subject of matter linked to the dataset",
        dataverse="pyDaRUS.Citation.subject",
    )

    authors: List[Author] = Field(
        description="Authors of this dataset.", default_factory=ListPlus
    )

    parameters: List[Parameter] = Field(
        description="Parameters to start and configure some process",
        default_factory=ListPlus,
    )

    id: str = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("rootINDEX"),
        xml="@id",
    )

    __repo__: Optional[str] = PrivateAttr(
        default="git://github.com/JR-1991/sdrdm-template.git"
    )

    __commit__: Optional[str] = PrivateAttr(
        default="2e6c6d2488a46547b281b6243038ee0eb29c0b97"
    )

    def add_to_authors(
        self, name: str, affiliation: Optional[str] = None, id: Optional[str] = None
    ) -> None:
        """
        Adds an instance of 'Author' to the attribute 'authors'.

        Args:


            id (str): Unique identifier of the 'Author' object. Defaults to 'None'.


            name (str): Full name including given and family name.


            affiliation (Optional[str]): To which organization the author is affiliated to. Defaults to None
        """

        params = {"name": name, "affiliation": affiliation}
        if id is not None:
            params["id"] = id
        authors = [Author(**params)]
        self.authors = self.authors + authors

    def add_to_parameters(
        self, key: str, value: float, id: Optional[str] = None
    ) -> None:
        """
        Adds an instance of 'Parameter' to the attribute 'parameters'.

        Args:


            id (str): Unique identifier of the 'Parameter' object. Defaults to 'None'.


            key (str): Name of the parameter.


            value (float): Respective value of a parameter.
        """

        params = {"key": key, "value": value}
        if id is not None:
            params["id"] = id
        parameters = [Parameter(**params)]
        self.parameters = self.parameters + parameters
