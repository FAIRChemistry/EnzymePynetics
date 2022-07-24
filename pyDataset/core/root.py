import sdRDM


from typing import Optional
from pydantic import PrivateAttr
from pydantic import Field
from typing import List
from .author import Author
from .parameter import Parameter


class Root(sdRDM.DataModel):

    __url__: Optional[str] = PrivateAttr(
        default="git://github.com/JR-1991/sdrdm-template.git"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="3ef119eb64e4d575b5ec473e67b3e34fd1f8f09f"
    )

    """This is the root of the data model and contains all objects defined in this example. While its good practice to have a single root, you can define as many roots as you like. Furthermore, the name does not have to be ```Root``` and can be any other name.
"""

    description: str = Field(
        ...,
        description="Describes the content of the dataset.",
        dataverse="pyDaRUS.Citation.description.text",
    )

    title: str = Field(
        ...,
        description="Title of the work",
        dataverse="pyDaRUS.Citation.title",
    )

    subject: List[str] = Field(
        description="Subject of matter linked to the dataset",
        dataverse="pyDaRUS.Citation.subject",
        default_factory=list,
    )

    authors: List[Author] = Field(
        description="Authors of this dataset.",
        default_factory=list,
    )

    parameters: List[Parameter] = Field(
        description="Parameters to start and configure some process",
        default_factory=list,
    )

    def add_to_authors(
        self,
        name: str,
        affiliation: str,
    ) -> None:
        """
        Adds an instance of 'Author' to the attribute 'authors'.

        Args:
            name (str): Full name including given and family name.
            affiliation (str): To which organization the author is affiliated to.
        """

        self.authors.append(
            Author(
                name=name,
                affiliation=affiliation,
            )
        )

    def add_to_parameters(
        self,
        key: str,
        value: float,
    ) -> None:
        """
        Adds an instance of 'Parameter' to the attribute 'parameters'.

        Args:
            key (str): Name of the parameter.
            value (float): Respective value of a parameter.
        """

        self.parameters.append(
            Parameter(
                key=key,
                value=value,
            )
        )
