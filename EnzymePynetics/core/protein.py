from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .abstractspecies import AbstractSpecies
from .sboterm import SBOTerm


@forge_signature
class Protein(AbstractSpecies):
    """This objects describes the proteins that were used or produced in the course of the experiment."""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("proteinINDEX"),
        xml="@id",
    )

    sequence: str = Field(
        ...,
        description="Amino acid sequence of the protein",
        template_alias="Sequence",
    )

    ecnumber: Optional[str] = Field(
        default=None,
        description="EC number of the protein.",
        pattern="(\\d+.)(\\d+.)(\\d+.)(\\d+)",
        template_alias="EC Number",
    )

    organism: Optional[str] = Field(
        default=None,
        description="Organism the protein was expressed in.",
        template_alias="Source organism",
    )

    organism_tax_id: Optional[str] = Field(
        default=None,
        description="Taxonomy identifier of the expression host.",
    )

    uniprotid: Optional[str] = Field(
        default=None,
        description=(
            "Unique identifier referencing a protein entry at UniProt. Use this"
            " identifier to initialize the object from the UniProt database."
        ),
        template_alias="UniProt ID",
    )

    ontology: SBOTerm = Field(
        description="None",
        default=SBOTerm.PROTEIN,
    )
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="70285185b8d9c7baf61e12dd52d943624695a510"
    )
