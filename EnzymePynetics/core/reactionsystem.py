import sdRDM

from typing import List, Optional
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .modelresult import ModelResult
from .reaction import Reaction
from .kineticmodel import KineticModel
from .sboterm import SBOTerm
from .reactionelement import ReactionElement


@forge_signature
class ReactionSystem(sdRDM.DataModel):

    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("reactionsystemINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Name of the reaction system",
    )

    reactions: List[Reaction] = Field(
        description="Reactions of the reaction system",
        default_factory=ListPlus,
        multiple=True,
    )

    result: Optional[ModelResult] = Field(
        default=ModelResult(),
        description="Result of the kinetic model fitting.",
    )

    def add_to_reactions(
        self,
        name: str,
        reversible: bool = False,
        temperature: Optional[float] = None,
        temperature_unit: Optional[str] = None,
        ph: Optional[float] = None,
        ontology: SBOTerm = SBOTerm.BIOCHEMICAL_REACTION,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        model: Optional[KineticModel] = None,
        educts: List[ReactionElement] = ListPlus(),
        products: List[ReactionElement] = ListPlus(),
        modifiers: List[ReactionElement] = ListPlus(),
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'Reaction' to attribute reactions

        Args:
            id (str): Unique identifier of the 'Reaction' object. Defaults to 'None'.
            name (): Name of the reaction..
            reversible (): Whether the reaction is reversible or irreversible. Defaults to False
            temperature (): Numeric value of the temperature of the reaction.. Defaults to None
            temperature_unit (): Unit of the temperature of the reaction.. Defaults to None
            ph (): PH value of the reaction.. Defaults to None
            ontology (): Ontology defining the role of the given species.. Defaults to SBOTerm.BIOCHEMICAL_REACTION
            uri (): URI of the reaction.. Defaults to None
            creator_id (): Unique identifier of the author.. Defaults to None
            model (): Kinetic model decribing the reaction.. Defaults to None
            educts (): List of educts containing ReactionElement objects.. Defaults to ListPlus()
            products (): List of products containing ReactionElement objects.. Defaults to ListPlus()
            modifiers (): List of modifiers (Proteins, snhibitors, stimulators) containing ReactionElement objects.. Defaults to ListPlus()
        """

        params = {
            "name": name,
            "reversible": reversible,
            "temperature": temperature,
            "temperature_unit": temperature_unit,
            "ph": ph,
            "ontology": ontology,
            "uri": uri,
            "creator_id": creator_id,
            "model": model,
            "educts": educts,
            "products": products,
            "modifiers": modifiers,
        }

        if id is not None:
            params["id"] = id

        self.reactions.append(Reaction(**params))

        return self.reactions[-1]
