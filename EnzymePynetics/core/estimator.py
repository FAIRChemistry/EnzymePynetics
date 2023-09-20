import sdRDM

from typing import Optional, Union, List
from pydantic import PrivateAttr, Field, validator
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from pydantic import StrictBool
from pydantic import PositiveFloat

from .modelresult import ModelResult
from .correlation import Correlation
from .reactionsystem import ReactionSystem
from .measurementdata import MeasurementData
from .protein import Protein
from .reaction import Reaction
from .abstractspecies import AbstractSpecies
from .reactionelement import ReactionElement
from .measurement import Measurement
from .parameter import Parameter
from .reactant import Reactant
from .kineticmodel import KineticModel
from .replicate import Replicate
from .kineticparameter import KineticParameter
from .datatypes import DataTypes
from .vessel import Vessel
from .sboterm import SBOTerm


@forge_signature
class Estimator(sdRDM.DataModel):

                
    """"""
    
    id: Optional[str] = Field(
            description="Unique identifier of the given object.",
            default_factory=IDGenerator("estimatorINDEX"),
            xml="@id"
    )
    
    name: Optional[str] = Field(
    
    
    default=None,
    
    description="Title of the kinetic experiment",
    
    )

    measured_reactant:: Optional[Reactant] = Field(
    
    
    default=None,
    
    description="Reactant that is measured in the experiment",
    
    )

    reaction_systems: List[ReactionSystem] = Field(
    
    
    description="Reactions of multiple species",
    
    default_factory=ListPlus,
    
    multiple=True,
    
    )

    species: List[Union[AbstractSpecies, Protein, Reactant]] = Field(
    
    
    description="Reactants, Inhibitor, Activators and Catalysts of the reaction",
    
    default_factory=ListPlus,
    
    multiple=True,
    
    )

    reaction: Optional[Reaction] = Field(
    
    
    default=None,
    
    description="Reaction proceeding in measurements",
    
    )

    models: List[KineticModel] = Field(
    
    
    description="Kinetic model options used to fit to measurement data.",
    
    default_factory=ListPlus,
    
    multiple=True,
    
    )

    measurements: List[Measurement] = Field(
    
    
    description="Measurement data for a given initial substrate concentration.",
    
    default_factory=ListPlus,
    
    multiple=True,
    
    )

    
    
    def add_to_reaction_systems(
        self,
        name: Optional[str] = None ,
        reactions: List[Reaction] = ListPlus() ,
        result: Optional[ModelResult] = None ,
        id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'ReactionSystem' to attribute reaction_systems

        Args:
            id (str): Unique identifier of the 'ReactionSystem' object. Defaults to 'None'.
            name (): Name of the reaction system. Defaults to None
            reactions (): Reactions of the reaction system. Defaults to ListPlus()
            result (): Result of the kinetic model fitting.. Defaults to None
        """

        params = {
            "name": name,
            "reactions": reactions,
            "result": result,
        }

        if id is not None:
            params["id"] = id

        self.reaction_systems.append(ReactionSystem(**params))

        return self.reaction_systems[-1]

    def add_abstract_species_to_species(
        self,
        name: str ,
        vessel_id: Vessel ,
        constant: StrictBool ,
        init_conc: Optional[float] = None ,
        unit: Optional[str] = None ,
        uri: Optional[str] = None ,
        creator_id: Optional[str] = None ,
        id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'AbstractSpecies' to attribute species

        Args:
            id (str): Unique identifier of the 'AbstractSpecies' object. Defaults to 'None'.
            name (): None.
            vessel_id (): None.
            constant (): None.
            init_conc (): None. Defaults to None
            unit (): None. Defaults to None
            uri (): None. Defaults to None
            creator_id (): None. Defaults to None
        """

        params = {
            "name": name,
            "vessel_id": vessel_id,
            "constant": constant,
            "init_conc": init_conc,
            "unit": unit,
            "uri": uri,
            "creator_id": creator_id,
        }

        if id is not None:
            params["id"] = id

        self.species.append(AbstractSpecies(**params))

        return self.species[-1]

    def add_protein_to_species(
        self,
        sequence: str ,
        name: str ,
        vessel_id: Vessel ,
        constant: StrictBool ,
        ecnumber: Optional[str] = None ,
        organism: Optional[str] = None ,
        organism_tax_id: Optional[str] = None ,
        uniprotid: Optional[str] = None ,
        ontology: SBOTerm = SBOTerm.CATALYST ,
        init_conc: Optional[float] = None ,
        unit: Optional[str] = None ,
        uri: Optional[str] = None ,
        creator_id: Optional[str] = None ,
        id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'Protein' to attribute species

        Args:
            id (str): Unique identifier of the 'Protein' object. Defaults to 'None'.
            sequence (): Amino acid sequence of the protein.
            name (): None.
            vessel_id (): None.
            constant (): None.
            ecnumber (): EC number of the protein.. Defaults to None
            organism (): Organism the protein was expressed in.. Defaults to None
            organism_tax_id (): Taxonomy identifier of the expression host.. Defaults to None
            uniprotid (): Unique identifier referencing a protein entry at UniProt. Use this identifier to initialize the object from the UniProt database.. Defaults to None
            ontology (): None. Defaults to SBOTerm.CATALYST
            init_conc (): None. Defaults to None
            unit (): None. Defaults to None
            uri (): None. Defaults to None
            creator_id (): None. Defaults to None
        """

        params = {
            "sequence": sequence,
            "name": name,
            "vessel_id": vessel_id,
            "constant": constant,
            "ecnumber": ecnumber,
            "organism": organism,
            "organism_tax_id": organism_tax_id,
            "uniprotid": uniprotid,
            "ontology": ontology,
            "init_conc": init_conc,
            "unit": unit,
            "uri": uri,
            "creator_id": creator_id,
        }

        if id is not None:
            params["id"] = id

        self.species.append(Protein(**params))

        return self.species[-1]

    def add_reactant_to_species(
        self,
        name: str ,
        vessel_id: Vessel ,
        constant: StrictBool ,
        smiles: Optional[str] = None ,
        inchi: Optional[str] = None ,
        chebi_id: Optional[str] = None ,
        ontology: SBOTerm = SBOTerm.SMALL_MOLECULE ,
        init_conc: Optional[float] = None ,
        unit: Optional[str] = None ,
        uri: Optional[str] = None ,
        creator_id: Optional[str] = None ,
        id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'Reactant' to attribute species

        Args:
            id (str): Unique identifier of the 'Reactant' object. Defaults to 'None'.
            name (): None.
            vessel_id (): None.
            constant (): None.
            smiles (): Simplified Molecular Input Line Entry System (SMILES) encoding of the reactant.. Defaults to None
            inchi (): International Chemical Identifier (InChI) encoding of the reactant.. Defaults to None
            chebi_id (): Unique identifier of the CHEBI database. Use this identifier to initialize the object from the CHEBI database.. Defaults to None
            ontology (): None. Defaults to SBOTerm.SMALL_MOLECULE
            init_conc (): None. Defaults to None
            unit (): None. Defaults to None
            uri (): None. Defaults to None
            creator_id (): None. Defaults to None
        """

        params = {
            "name": name,
            "vessel_id": vessel_id,
            "constant": constant,
            "smiles": smiles,
            "inchi": inchi,
            "chebi_id": chebi_id,
            "ontology": ontology,
            "init_conc": init_conc,
            "unit": unit,
            "uri": uri,
            "creator_id": creator_id,
        }

        if id is not None:
            params["id"] = id

        self.species.append(Reactant(**params))

        return self.species[-1]

    def add_to_models(
        self,
        name: str ,
        equation: str ,
        parameters: List[KineticParameter] = ListPlus() ,
        ontology: Optional[SBOTerm] = None ,
        id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'KineticModel' to attribute models

        Args:
            id (str): Unique identifier of the 'KineticModel' object. Defaults to 'None'.
            name (): Name of the kinetic law..
            equation (): Equation for the kinetic law..
            parameters (): List of estimated parameters.. Defaults to ListPlus()
            ontology (): Type of the estimated parameter.. Defaults to None
        """

        params = {
            "name": name,
            "equation": equation,
            "parameters": parameters,
            "ontology": ontology,
        }

        if id is not None:
            params["id"] = id

        self.models.append(KineticModel(**params))

        return self.models[-1]

    def add_to_measurements(
        self,
        name: str ,
        temperature: float ,
        temperature_unit: str ,
        ph: float ,
        global_time_unit: str ,
        species: List[MeasurementData] = ListPlus() ,
        global_time: List[float] = ListPlus() ,
        uri: Optional[str] = None ,
        creator_id: Optional[str] = None ,
        id: Optional[str] = None
    ) -> None:
        """
        This method adds an object of type 'Measurement' to attribute measurements

        Args:
            id (str): Unique identifier of the 'Measurement' object. Defaults to 'None'.
            name (): Name of the measurement.
            temperature (): Numeric value of the temperature of the reaction..
            temperature_unit (): Unit of the temperature of the reaction..
            ph (): PH value of the reaction..
            global_time_unit (): Unit of the global time..
            species (): Species of the measurement.. Defaults to ListPlus()
            global_time (): Global time of the measurement all replicates agree on.. Defaults to ListPlus()
            uri (): URI of the reaction.. Defaults to None
            creator_id (): Unique identifier of the author.. Defaults to None
        """

        params = {
            "name": name,
            "temperature": temperature,
            "temperature_unit": temperature_unit,
            "ph": ph,
            "global_time_unit": global_time_unit,
            "species": species,
            "global_time": global_time,
            "uri": uri,
            "creator_id": creator_id,
        }

        if id is not None:
            params["id"] = id

        self.measurements.append(Measurement(**params))

        return self.measurements[-1]
