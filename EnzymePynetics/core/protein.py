
from typing import Optional
from pydantic import Field, PrivateAttr
from sdRDM.base.utils import forge_signature, IDGenerator
from .sboterm import SBOTerm
from .abstractspecies import AbstractSpecies


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
        default=SBOTerm.CATALYST,
    )
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="1eb6fe3fd3d8f9a2eb20911141e5b01d10a9bd57"
    )
