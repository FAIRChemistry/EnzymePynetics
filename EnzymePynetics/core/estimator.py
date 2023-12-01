import sdRDM

import numpy as np
import pandas as pd
import plotly.express as px
from typing import Optional, Union, List
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator
from itertools import combinations
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from IPython.display import display
from tqdm import tqdm
from EnzymePynetics.ioutils import parse_enzymeml, _to_enzymeml
from lmfit import report_fit
from .modelresult import ModelResult
from .abstractspecies import AbstractSpecies
from .reactionelement import ReactionElement
from .reactionsystem import ReactionSystem
from .kineticparameter import KineticParameter
from .sboterm import SBOTerm
from .vessel import Vessel
from .reaction import Reaction
from .protein import Protein
from .reactant import Reactant
from .kineticmodel import KineticModel
from .measurementdata import MeasurementData
from .measurement import Measurement
from .paramtype import ParamType


@forge_signature
class Estimator(sdRDM.DataModel):
    """"""

    id: Optional[str] = Field(
        description="Unique identifier of the given object.",
        default_factory=IDGenerator("estimatorINDEX"),
        xml="@id",
    )

    name: Optional[str] = Field(
        default=None,
        description="Title of the kinetic experiment",
    )

    measured_reactant: Optional[Reactant] = Field(
        default=None,
        description="Reactant that is measured in the experiment",
    )

    reaction_systems: List[ReactionSystem] = Field(
        description="Reactions of multiple species",
        default_factory=ListPlus,
        multiple=True,
    )

    species: List[Union[Protein, Reactant]] = Field(
        description="Reactants, Inhibitor, Activators and Catalysts of the reaction",
        default_factory=ListPlus,
        multiple=True,
    )

    reactions: List[Reaction] = Field(
        description="Reactions proceeding in measurements",
        default_factory=ListPlus,
        multiple=True,
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
    __repo__: Optional[str] = PrivateAttr(
        default="https://github.com/haeussma/EnzymePynetics"
    )
    __commit__: Optional[str] = PrivateAttr(
        default="70285185b8d9c7baf61e12dd52d943624695a510"
    )

    def add_to_reaction_systems(
        self,
        name: Optional[str] = None,
        reactions: List[Reaction] = ListPlus(),
        result: Optional[ModelResult] = None,
        id: Optional[str] = None,
    ) -> None:
        """
        This method adds an object of type 'ReactionSystem' to attribute reaction_systems

        Args:
            id (str): Unique identifier of the 'ReactionSystem' object. Defaults to 'None'.
            name (): Name of the reaction system. Defaults to None
            reactions (): Reactions of the reaction system. Defaults to ListPlus()
            result (): Result of the kinetic model fitting.. Defaults to None
        """
        params = {"name": name, "reactions": reactions, "result": result}
        if id is not None:
            params["id"] = id
        self.reaction_systems.append(ReactionSystem(**params))
        return self.reaction_systems[-1]

    def add_protein_to_species(
        self,
        sequence: str,
        name: str,
        vessel_id: Vessel,
        constant: bool,
        ecnumber: Optional[str] = None,
        organism: Optional[str] = None,
        organism_tax_id: Optional[str] = None,
        uniprotid: Optional[str] = None,
        ontology: SBOTerm = SBOTerm.CATALYST,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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
        name: str,
        vessel_id: Vessel,
        constant: bool,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        chebi_id: Optional[str] = None,
        ontology: SBOTerm = SBOTerm.SMALL_MOLECULE,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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

    def add_to_models(
        self,
        name: str,
        equation: str,
        parameters: List[KineticParameter] = ListPlus(),
        ontology: Optional[SBOTerm] = None,
        id: Optional[str] = None,
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
        name: str,
        temperature: float,
        temperature_unit: str,
        ph: float,
        global_time_unit: str,
        species: List[MeasurementData] = ListPlus(),
        global_time: List[float] = ListPlus(),
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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

    def add_protein_to_species(
        self,
        sequence: str,
        name: str,
        vessel_id: Vessel,
        constant: bool,
        ecnumber: Optional[str] = None,
        organism: Optional[str] = None,
        organism_tax_id: Optional[str] = None,
        uniprotid: Optional[str] = None,
        ontology: SBOTerm = SBOTerm.CATALYST,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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
        name: str,
        vessel_id: Vessel,
        constant: bool,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        chebi_id: Optional[str] = None,
        ontology: SBOTerm = SBOTerm.SMALL_MOLECULE,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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

    def add_to_models(
        self,
        name: str,
        equation: str,
        parameters: List[KineticParameter] = ListPlus(),
        ontology: Optional[SBOTerm] = None,
        id: Optional[str] = None,
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
        name: str,
        temperature: float,
        temperature_unit: str,
        ph: float,
        global_time_unit: str,
        species: List[MeasurementData] = ListPlus(),
        global_time: List[float] = ListPlus(),
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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

    def add_abstract_species_to_species(
        self,
        name: str,
        vessel_id: Vessel,
        constant: bool,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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
        sequence: str,
        name: str,
        vessel_id: Vessel,
        constant: bool,
        ecnumber: Optional[str] = None,
        organism: Optional[str] = None,
        organism_tax_id: Optional[str] = None,
        uniprotid: Optional[str] = None,
        ontology: SBOTerm = SBOTerm.CATALYST,
        init_conc: Optional[float] = None,
        unit: Optional[str] = None,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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

    def add_to_models(
        self,
        name: str,
        equation: str,
        parameters: List[KineticParameter] = ListPlus(),
        ontology: Optional[SBOTerm] = None,
        id: Optional[str] = None,
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
        name: str,
        temperature: float,
        temperature_unit: str,
        ph: float,
        global_time_unit: str,
        species: List[MeasurementData] = ListPlus(),
        global_time: List[float] = ListPlus(),
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        id: Optional[str] = None,
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

    def add_reaction(
        self,
        id: str,
        name: str,
        reversible: bool = False,
        ontology: SBOTerm = SBOTerm.BIOCHEMICAL_REACTION,
        uri: Optional[str] = None,
        creator_id: Optional[str] = None,
        model: Optional[KineticModel] = None,
        educt: Reactant = None,
        product: Reactant = None,
        catalyst: Protein = None,
        inhibitor: Reactant = None,
    ) -> None:
        """
        This method adds an object of type 'Reaction' to attribute reaction

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

        # add educt and product to ReactionElements
        if educt:
            educt_reaction_element = [
                ReactionElement(
                    species_id=educt.id,
                    constant=educt.constant,
                    ontology=SBOTerm.SUBSTRATE,
                )
            ]
        else:
            educt_reaction_element = ListPlus()

        if product:
            product_reaction_element = [
                ReactionElement(
                    species_id=product.id,
                    constant=product.constant,
                    ontology=SBOTerm.PRODUCT,
                )
            ]
        else:
            product_reaction_element = ListPlus()

        # Add modifiers
        modifiers = []

        if not catalyst and len(self.enzymes) == 1:
            catalyst = self.enzymes[0]

        if catalyst:
            modifiers.append(
                ReactionElement(
                    species_id=catalyst.id,
                    constant=catalyst.constant,
                    ontology=SBOTerm.CATALYST,
                )
            )
        else:
            raise ValueError("No protein defined as enzyme. Use 'add_protein' first.")

        if inhibitor:
            modifiers.append(
                ReactionElement(
                    species_id=inhibitor.id,
                    constant=inhibitor.constant,
                    ontology=SBOTerm.INHIBITOR,
                )
            )

        params = {
            "name": name,
            "reversible": reversible,
            "temperature": self.temperature,
            "temperature_unit": self.temperature_unit,
            "ph": self.ph,
            "ontology": ontology,
            "uri": uri,
            "creator_id": creator_id,
            "model": model,
            "educts": educt_reaction_element,
            "products": product_reaction_element,
            "modifiers": modifiers,
        }

        if id is not None:
            params["id"] = id

        new_reaction = Reaction(**params)

        if any([reaction.id == new_reaction.id for reaction in self.reactions]):
            self.reactions = [
                new_reaction if reaction.id == new_reaction.id else reaction
                for reaction in self.reactions
            ]

            return new_reaction

        else:
            self.reactions.append(new_reaction)

            return new_reaction

    def add_model(
        self,
        id: str,
        name: str,
        equation: str,
        parameters: List[KineticParameter] = ListPlus(),
        ontology: Optional[SBOTerm] = None,
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
            "id": id,
            "name": name,
            "equation": equation,
            "parameters": parameters,
            "ontology": ontology,
        }

        new_model = KineticModel(**params)
        self._init_parameters(new_model)

        if any([model.id == new_model.id for model in self.models]):
            self.models = [
                new_model if model.id == new_model.id else model
                for model in self.models
            ]

            return new_model

        else:
            self.models.append(new_model)

            return new_model

    def _create_model_combinations(self):
        if len(self.reactions) > 1:
            raise ValueError(
                "Currently only one reaction per reaction system is supported."
            )
        self.reaction_systems = ListPlus()
        for substrate_model in self.substrate_models:
            # Add different substrate models to reaction
            substrate_reaction = Reaction(**self.reactions[0].to_dict())
            substrate_reaction.model = KineticModel(**substrate_model.to_dict())

            self.reaction_systems.append(
                ReactionSystem(
                    name=f"{substrate_model.name}",
                    reactions=[substrate_reaction],
                )
            )

            if self.enzyme_models:
                inactivation_reaction = self._create_inactivation_reaction()
                # create reaction system with enzyme models
                for enzyme_model in self.enzyme_models:
                    substrate_reaction = Reaction(**substrate_reaction.to_dict())
                    new_inactivation = Reaction(**inactivation_reaction.to_dict())
                    new_inactivation.model = KineticModel(**enzyme_model.to_dict())

                    self.reaction_systems.append(
                        ReactionSystem(
                            name=f"{substrate_model.name} with {enzyme_model.name}",
                            reactions=[substrate_reaction, new_inactivation],
                        )
                    )

    def _create_inactivation_reaction(self) -> Reaction:
        # Make copies of enzyme species
        if self.enzyme.ontology == SBOTerm.SMALL_MOLECULE.value:
            Species = Reactant
        else:
            Species = Protein
        active = Species(**self.enzyme.to_dict())
        active.name += " (active)"
        active.constant = False

        inactive = Species(**self.enzyme.to_dict())
        inactive.name += " (inactive)"
        inactive.id += "_inactive"
        inactive.constant = False
        self.species.append(inactive)

        # Make enzyme inactivation reaction
        inactivation = Reaction(**self.reactions[0].to_dict())
        inactivation.name = "enzyme inactivation"
        inactivation.reversible = False
        inactivation.educts = [
            ReactionElement(
                species_id=active.id,
                constant=False,
                ontology=SBOTerm.CATALYST,
            )
        ]
        inactivation.products = [
            ReactionElement(
                species_id=inactive.id,
                constant=False,
                ontology=SBOTerm.PROTEIN,  # , since inactive
            )
        ]

        return inactivation

    def _remove_nans(self):
        """Removes all samples which contain nan values."""

        nan_mask = np.isnan(self.substrate_data).any(axis=1)

        substrate_data = self.substrate_data[~nan_mask]
        product_data = self.product_data[~nan_mask]
        enzyme_data = self.enzyme_data[~nan_mask]
        time_data = self.time_data[~nan_mask]

        return [substrate_data, enzyme_data, product_data, time_data]

    def _subset_time(
        self,
        min_time: float,
        max_time: float,
        substrate_data: np.ndarray,
        enzyme_data: np.ndarray,
        product_data: np.ndarray,
        time_data: np.ndarray,
    ):
        if min_time is None and max_time is not None:
            subset_slice = time_data < max_time
        elif min_time is not None and max_time is None:
            subset_slice = time_data > min_time
        elif min_time is not None and max_time is not None:
            subset_slice = (time_data > min_time) & (time_data < max_time)
        else:
            raise ValueError("Either min_time or max_time must be specified.")

        return (
            substrate_data[subset_slice].reshape(time_data.shape[0], -1),
            enzyme_data[subset_slice].reshape(time_data.shape[0], -1),
            product_data[subset_slice].reshape(time_data.shape[0], -1),
            time_data[subset_slice].reshape(time_data.shape[0], -1),
        )

    def fit_models_fixed_params(
        self,
        model: KineticModel,
        max_time: float = None,
        fixed_params: List[str] = None,
    ):
        if isinstance(model, str):
            model = [
                system for system in self.reaction_systems if system.name == model
            ][0]

        for system in self.reaction_systems:
            if system.name != model.name:
                for reaction in system.reactions:
                    for param in reaction.model.parameters:
                        if param.name in fixed_params:
                            param.value = model.get_parameter(param.name).value

        substrate, enzyme, product, time = self._remove_nans()

        if max_time:
            substrate, enzyme, product, time = self._subset_time(
                max_time, substrate, enzyme, product, time
            )

        systems = tqdm(self.reaction_systems)
        for system in systems:
            systems.set_description(desc=f"Fitting {system.name} model")
            report = system.fit(
                substrate_data=substrate,
                enzyme_data=enzyme,
                product_data=product,
                times=time,
                fixed_params=fixed_params,
            )
            # print(report_fit(report))

        self.reaction_systems.sort(key=lambda x: x.result.AIC)

        display(self.fit_statistics())

    def fit_models(self, min_time: float = None, max_time: float = None):
        self._create_model_combinations()

        substrate, enzyme, product, time = self._remove_nans()

        if min_time or max_time:
            substrate, enzyme, product, time = self._subset_time(
                min_time, max_time, substrate, enzyme, product, time
            )

        systems = tqdm(self.reaction_systems)
        for system in systems:
            systems.set_description(desc=f"Fitting {system.name} model")
            report = system.fit(
                substrate_data=substrate,
                enzyme_data=enzyme,
                product_data=product,
                times=time,
            )
            print(report_fit(report))

        self.reaction_systems.sort(key=lambda x: x.result.AIC)

        display(self.fit_statistics())

    def fit_statistics(self):
        header = np.array(
            [
                ["Model", ""],
                ["AIC", ""],
                [ParamType.K_CAT.value, f"1 / {self.time_unit}"],
                [ParamType.K_M.value, self.substrate_unit],
                [
                    f"{ParamType.K_CAT.value} / {ParamType.K_M.value}",
                    f"{self.substrate_unit} / {self.time_unit}",
                ],
                [ParamType.K_IC.value, self.substrate_unit],
                [ParamType.K_IU.value, self.substrate_unit],
                [ParamType.K_IE.value, f"1 / {self.time_unit}"],
            ]
        )

        entries = []
        for system in self.reaction_systems:
            entry = dict.fromkeys(header[:, 0])
            entry["Model"] = system.name
            entry["AIC"] = system.result.AIC

            for reaction in system.reactions:
                for param in reaction.model.parameters:
                    entry[param.name] = param.value

            entries.append(entry)

        decimal_formatting = {
            "AIC": "{:.0f}",
            ParamType.K_M.value: "{3f}",
            ParamType.K_CAT.value: "{:.3f}",
            ParamType.K_IE.value: "{:.3f}",
            ParamType.K_IU.value: "{:.3f}",
            ParamType.K_IC.value: "{:.3f}",
        }

        df = pd.DataFrame(entries).set_index("Model").sort_values("AIC")
        df.columns = pd.MultiIndex.from_arrays(header[1:, :].T)

        df[f"{ParamType.K_CAT.value} / {ParamType.K_M.value}"] = (
            df[ParamType.K_CAT.value].values / df[ParamType.K_M.value].values
        )

        return (
            df.style.format("{:.3f}", na_rep="")
            .format("{:.0f}", subset=["AIC"], na_rep="failed")
            .background_gradient(cmap="Blues", subset=["AIC"], gmap=-df["AIC"])
        )

    def _format_html(self, param: str):
        if param == f"{ParamType.K_CAT.value} {ParamType.K_M.value}":
            param1, param2 = param.split()
            return f"{param1.replace('_', '<sub>') + '</sub>'} {param2.replace('_', '<sub>') + '</sub>'}<sup>-1</sup>"

        return param.replace("_", "<sub>") + "</sub>"

    def model_table(self, round_digits: int = 3, return_fig: bool = False):
        table_traces = self._get_table_traces(round_digits)

        fig = go.Figure(data=table_traces)
        if return_fig:
            return fig

        fig.update_layout(
            title=(
                f"{self.name} at {self.temperature} {self.temperature_unit} and pH"
                f" {self.ph}"
            ),
        )

        fig.show()

    def _get_table_traces(self, round_digits: int = 3):
        # Get column labels

        # Iterate
        parameters = {}
        for system in self.reaction_systems:
            for reaction in system.reactions:
                for param in reaction.model.parameters:
                    parameters[param.name] = param.unit

        columns = ["Model", "AIC"] + list(parameters.keys())
        columns.append(f"{ParamType.K_CAT.value} {ParamType.K_M.value}")

        headers = []
        for column in columns:
            if column == "Model":
                col = ["", f"<b>{column}</b>"]
                headers.append(col)
            elif column == "AIC":
                col = ["", f"<b>{column}</b>"]
                headers.append(col)
            elif column == f"{ParamType.K_CAT.value} {ParamType.K_M.value}":
                col = [
                    (
                        f"<b>{self._format_html(column)}</b><br>{ReactionSystem._format_unit(parameters[ParamType.K_CAT.value])} {ReactionSystem._format_unit(parameters[ParamType.K_M.value])}"
                    ),
                    "",
                ]
                headers.append(col)
            else:
                col = [
                    f"<b>{self._format_html(column)}</b><br>{ReactionSystem._format_unit(parameters[column])}",
                    "",
                ]
                headers.append(col)

        entries = []
        for system in self.reaction_systems:
            data = []
            for row_name in columns:
                if row_name == "Model":
                    data.append(system.name)
                    continue
                if row_name == "AIC":
                    data.append(round(system.result.AIC))
                    continue

                if row_name not in [param.name for param in system.result.parameters]:
                    data.append(float("nan"))
                    continue

                for reaction in system.reactions:
                    for parameter in reaction.model.parameters:
                        if parameter.name != row_name:
                            continue
                        percent_error = parameter.stdev / parameter.value * 100
                        if percent_error > 100:
                            percent_error = f"± >100 %"
                        else:
                            percent_error = (
                                f"± {percent_error:.0f} %" if percent_error else ""
                            )
                        data.append(
                            f"{round(parameter.value, round_digits)}<br>{percent_error}"
                        )
            entries.append(data)
        entries = np.array(entries).T
        entries[entries == "nan"] = ""

        return go.Table(
            header=dict(values=headers, align="center"),
            cells=dict(values=entries, align="center"),
        )

    def _get_correlation_traces(self):
        # Format model names
        model_names = [system.name for system in self.reaction_systems]
        for name_id, name in enumerate(model_names):
            if len(name.split()) > 3:
                model_names[name_id] = (
                    " ".join(name.split()[:3]) + "<br>" + " ".join(name.split()[3:])
                )

            # Format correlation labels
            param_name_map = {
                ParamType.K_CAT.value: "k<sub>cat</sub>",
                ParamType.K_M.value: "K<sub>m</sub>",
                ParamType.K_IE.value: "k<sub>ie</sub>",
                ParamType.K_IC.value: "K<sub>ic</sub>",
                ParamType.K_IU.value: "K<sub>iu</sub>",
            }

        params = [param.value for param in ParamType]
        unique_combinations = list(combinations(params, 2))

        # Add data to heatmap
        correlations = []
        for system in self.reaction_systems:
            system_entry = []
            for label in unique_combinations:
                key_label, value_label = label

                if key_label not in [param.name for param in system.result.parameters]:
                    system_entry.append(float("nan"))
                    continue

                for parameter in system.result.parameters:
                    if key_label != parameter.name:
                        continue

                    if value_label not in [
                        corr.parameter_name for corr in parameter.correlations
                    ]:
                        system_entry.append(float("nan"))
                        continue

                    for other_correlation in parameter.correlations:
                        if other_correlation.parameter_name != value_label:
                            continue
                        system_entry.append(other_correlation.value)

                    if parameter.name != key_label:
                        continue

                    for correlation in parameter.correlations:
                        if correlation.parameter_name != key_label:
                            continue

                        system_entry.append(correlation.value)

            correlations.append(system_entry)

        correlations = np.array(correlations)[::-1]
        model_names = list(reversed(model_names))
        text_values = np.around(correlations, decimals=2).astype(str)
        text_values[text_values == "nan"] = ""

        for label_id, label in enumerate(unique_combinations):
            key, value = label
            unique_combinations[
                label_id
            ] = f"{param_name_map[key]}<br>{param_name_map[value]}"

        return go.Heatmap(
            z=correlations,
            x=unique_combinations,
            y=model_names,
            colorscale=px.colors.diverging.balance,
            zmin=-1,
            zmax=1,
            text=text_values,
            texttemplate="%{text}",
        )

    def correlations(self, return_fig: bool = False):
        fig = go.Figure(data=self._get_correlation_traces())

        fig.update_layout(
            template="simple_white",
            xaxis_title="Parameter correlations",
            yaxis_title="Model",
            yaxis=dict(showgrid=False),
            xaxis=dict(showgrid=False),
            title=(
                "Correlations of estimated parameters at"
                f" {self.temperature} {self.temperature_unit} and pH {self.ph}"
            ),
        )

        if return_fig:
            return fig

        fig.show()

    def _init_parameters(self, model: KineticModel):
        value = float("nan")

        for param in model.eq_parameters:
            if param == ParamType.K_CAT.value:
                ontology = SBOTerm.K_CAT
                initial_value = self._init_kcat
                unit = f"1 / {self.time_unit}"
                upper = initial_value * 10
                lower = initial_value / 100

            elif param == ParamType.K_M.value:
                ontology = SBOTerm.K_M
                initial_value = self._init_km
                unit = self.substrate_unit
                upper = initial_value * 10
                lower = initial_value / 1000

            elif param == ParamType.K_IC.value:
                initial_value = self._init_km
                unit = self.substrate_unit
                ontology = None
                upper = initial_value * 5
                lower = initial_value / 100

            elif param == ParamType.K_IU.value:
                initial_value = self._init_km
                unit = self.substrate_unit
                ontology = None
                upper = initial_value * 10
                lower = initial_value / 100

            elif param == ParamType.K_IE.value:
                if self.time_unit == "s" or self.time_unit == "sec":
                    initial_value = np.log(2) / 60 / 60  # 60 min half life
                if self.time_unit == "min":
                    initial_value = np.log(2) / 60  # 60 min half life
                if self.time_unit == "h":
                    initial_value = np.log(2)
                unit = f"1 / {self.time_unit}"
                ontology = None
                upper = initial_value * 50
                lower = initial_value / 50

            else:
                raise ValueError(f"Parameter '{param}' not recognized.")

            model.add_to_parameters(
                name=param,
                initial_value=initial_value,
                value=value,
                unit=unit,
                ontology=ontology,
                upper=upper,
                lower=lower,
            )

    def to_enzymeml(
        self,
        enzymeml: "EnzymeML.EnzymeMLDocument",
        reaction_system: ReactionSystem = None,
        out_path: str = None,
    ) -> "EnzymeML.EnzymeMLDocument":
        return _to_enzymeml(enzymeml, reaction_system, out_path)

    def remove_replicate(self, replicate_id: str):
        for measurement in self.measurements:
            for species in measurement.species:
                if not species.replicates:
                    continue
                for replicate in species.replicates:
                    if replicate.id == replicate_id:
                        species.replicates.remove(replicate)

    @property
    def ph(self):
        if not all(
            [
                measurement.ph == self.measurements[0].ph
                for measurement in self.measurements
            ]
        ):
            raise ValueError("Measurements have inconsistent pH values.")
        return self.measurements[0].ph

    @property
    def temperature(self):
        if not all(
            [
                measurement.temperature == self.measurements[0].temperature
                for measurement in self.measurements
            ]
        ):
            raise ValueError("Measurements have inconsistent temperature values.")
        return self.measurements[0].temperature

    @property
    def temperature_unit(self):
        if not all(
            [
                measurement.temperature_unit == self.measurements[0].temperature_unit
                for measurement in self.measurements
            ]
        ):
            raise ValueError("Measurements have inconsistent temperature unit values.")
        return self.measurements[0].temperature_unit

    @property
    def time_unit(self):
        if not all(
            [
                measurement.global_time_unit == self.measurements[0].global_time_unit
                for measurement in self.measurements
            ]
        ):
            raise ValueError("Measurements have inconsistent time units.")
        return self.measurements[0].global_time_unit

    @property
    def substrate_unit(self):
        return self._get_consistent_unit(self.substrate)

    @property
    def product_unit(self):
        return self._get_consistent_unit(self.product)

    @property
    def enzyme_unit(self):
        return self._get_consistent_unit(self.enzyme)

    @property
    def inhibitor_unit(self):
        return self._get_consistent_unit(self.inhibitor)

    def _get_consistent_unit(self, species: AbstractSpecies) -> None:
        units = [measurement.unit for measurement in self._get_species_data(species)]
        if not all([unit == units[0] for unit in units]):
            raise ValueError("Measurements have inconsistent substrate units.")

        return units[0]

    @property
    def reactants(self):
        return [species for species in self.species if isinstance(species, Reactant)]

    @property
    def modifiers(self):
        return [species for species in self.species if species.constant == True]

    @property
    def enzymes(self):
        return [
            species
            for species in self.species
            if species.ontology == SBOTerm.CATALYST.value
        ]

    @property
    def substrate_models(self):
        return [
            model for model in self.models if model.equation.startswith("substrate")
        ]

    @property
    def enzyme_models(self):
        return [model for model in self.models if model.equation.startswith("catalyst")]

    @property
    def measured_reactant_role(self) -> SBOTerm:
        for reaction in self.reactions:
            for educt in reaction.educts:
                if educt.species_id == self.measured_reactant.id:
                    return SBOTerm(educt.ontology)

            for product in reaction.products:
                if product.species_id == self.measured_reactant.id:
                    return SBOTerm(product.ontology)

        raise ValueError(
            f"Measured reactant '{self.measured_reactant}' not found in the defined"
            " reaction."
        )

    @property
    def substrate(self):
        return self._get_species_of_role(SBOTerm.SUBSTRATE)

    @property
    def product(self):
        return self._get_species_of_role(SBOTerm.PRODUCT)

    @property
    def enzyme(self):
        return self._get_species_of_role(SBOTerm.CATALYST)

    @property
    def inhibitor(self):
        return self._get_species_of_role(SBOTerm.INHIBITOR)

    @property
    def init_substrate_data(self):
        init_substrates = []
        for n_replicates, measurement in zip(
            self._measurement_replicates, self.measurements
        ):
            for data in measurement.species:
                if data.species_id == self.substrate.id:
                    init_substrates.append([data.init_conc] * n_replicates)

        flat_list = np.array([item for sublist in init_substrates for item in sublist])

        return flat_list

    @property
    def substrate_data(self):
        if self.measured_reactant_role == SBOTerm.SUBSTRATE:
            return self._get_measured_data(self.substrate)
        else:
            return self._calculate_missing_reactant(self.product)

    @property
    def product_data(self):
        if self.measured_reactant_role == SBOTerm.PRODUCT:
            return self._get_measured_data(self.product)
        else:
            return self._calculate_missing_reactant(self.substrate)

    def _get_measured_data(self, reactant: Reactant):
        measurement_data = []
        for measurement in self._get_species_data(reactant):
            for replicate in measurement.replicates:
                measurement_data.append(replicate.data)

        return np.array(measurement_data).reshape(sum(self._measurement_replicates), -1)

    def _calculate_missing_reactant(self, existing_reactant: Reactant):
        # calculate_product
        if existing_reactant == self.substrate:
            return self.init_substrate_data[:, None] - self.substrate_data

        # calculate substrate
        else:
            return self.init_substrate_data[:, None] - self.product_data

    @property
    def _measurement_replicates(self) -> List[int]:
        measurement_replicates = []
        for measurement in self.measurements:
            for data in measurement.species:
                if data.species_id == self.measured_reactant.id:
                    measurement_replicates.append(len(data.replicates))

        return measurement_replicates

    @property
    def time_data(self):
        time_data = []
        for measurement in self._get_species_data(self.measured_reactant):
            for replicate in measurement.replicates:
                time_data.append(replicate.time)

        return np.array(time_data).reshape(sum(self._measurement_replicates), -1)

    @property
    def enzyme_data(self):
        enzyme_data = []
        for n_reps, measurement in zip(self._measurement_replicates, self.measurements):
            for data in measurement.species:
                if not data.species_id == self.enzyme.id:
                    continue

                if not data.replicates:
                    enzyme_data.append([data.init_conc] * n_reps)

        flat_list = np.array([item for sublist in enzyme_data for item in sublist])

        return np.repeat(np.array(flat_list), self.time_data.shape[1]).reshape(
            self.time_data.shape
        )

    @property
    def inhibitor_data(self):
        inhibitor_data = []
        for n_reps, measurement in zip(self._measurement_replicates, self.measurements):
            for data in measurement.species:
                if not data.species_id == self.inhibitor.id:
                    continue

                if not data.replicates:
                    inhibitor_data.append([data.init_conc] * n_reps)

        return np.repeat(np.array(inhibitor_data), self.time_data.shape[1]).reshape(
            self.time_data.shape
        )

    def _get_species_of_role(self, role: SBOTerm):
        for reaction in self.reactions:
            for educt in reaction.educts:
                if educt.ontology == role.value:
                    species_id = educt.species_id

            for product in reaction.products:
                if product.ontology == role.value:
                    species_id = product.species_id

            for modifier in reaction.modifiers:
                if modifier.ontology == role.value:
                    species_id = modifier.species_id

        return self._get_species(species_id)

    def _get_species(self, species_id: str):
        for species in self.species:
            if species.id == species_id:
                return species

    def _get_substrate_rates(self):
        substrsate_rates = np.diff(self.substrate_data)
        time_rates = np.diff(self.time_data)

        return np.abs(substrsate_rates / time_rates)

    @property
    def _init_kcat(self):
        substrate_rates = self._get_substrate_rates()
        normalized_rates = (
            substrate_rates / self.enzyme_data[:, : substrate_rates.shape[1]]
        )

        return np.nanmax(normalized_rates)

    @property
    def _init_km(self):
        return max(self.init_substrate_data) / 2

    @property
    def substrate_unit(self):
        for measurment in self.measurements:
            for species in measurment.species:
                if species.species_id == self.substrate.id:
                    return species.unit

    def _get_species_data(self, species: AbstractSpecies) -> MeasurementData:
        species_measurements = []
        for measurment in self.measurements:
            for data in measurment.species:
                if data.species_id == species.id:
                    species_measurements.append(data)
        return species_measurements

    @classmethod
    def from_enzymeml(
        cls,
        enzymeml: Union[str, "EnzymeMLDocument"],
        measured_reactant: Union[Reactant, str],
    ):
        return parse_enzymeml(cls, enzymeml, measured_reactant)

    def get_reaction_system(self, system_name: str) -> ReactionSystem:
        for system in self.reaction_systems:
            if system.name == system_name:
                return system

        raise ValueError(f"Reaction system '{system_name}' not found.")

    def visualize_data(self):
        n_cols = int(np.ceil(np.sqrt(len(self.measurements))))

        fig = make_subplots(
            rows=n_cols,
            cols=n_cols,
            shared_yaxes=True,
            subplot_titles=[f"{measurement.name}" for measurement in self.measurements],
            y_title=(
                f"{self.measured_reactant.name} {self._format_unit(self.substrate_unit)}"
            ),
            x_title=f"Time / {self.time_unit}",
        )

        row = (0,)
        col = (0,)
        for meas_id, measurement in enumerate(self.measurements):
            row = meas_id // n_cols
            col = meas_id % n_cols
            for species in measurement.species:
                if not species.replicates:
                    continue

                n_colors = len(species.replicates)
                colors = px.colors.sample_colorscale(
                    "Purpor_r", [n / n_colors for n in range(n_colors)]
                )

                for replicate, color in zip(species.replicates, colors):
                    fig.add_trace(
                        go.Scatter(
                            x=np.array(replicate.time),
                            y=np.array(replicate.data),
                            mode="markers",
                            marker=dict(color=color),
                            name=replicate.id,
                        ),
                        row=row + 1,
                        col=col + 1,
                    )

        init_substrate_conc = [
            f"{self.substrate.name} {substrate.init_conc} /"
            f" {self._format_unit(substrate.unit)}"
            for substrate in self._get_species_data(self.substrate)
        ]

        for annotation, subtitle in zip(fig.layout.annotations, init_substrate_conc):
            annotation.update(text=subtitle)
            annotation["font"].update(size=7)

        fig.update_yaxes(
            # title=f"{self.measured_reactant.name} {self.measured_reactant.unit}",
            ticks="outside",
            tickwidth=1,
            tickcolor="black",
            tickfont=dict(size=7),
            title_font=dict(size=7),
        )

        fig.update_xaxes(
            # title=f"Time {self.time_unit}",
            ticks="outside",
            tickwidth=1,
            tickcolor="black",
            tickfont=dict(size=8),
            title_font=dict(size=8),
        )

        fig.update_layout(
            title=f"{self.name}",
            template="simple_white",
            hoverlabel_align="right",
            legend_title=f"Measurement IDs",
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False),
        )
        fig.show()

    def visualize(self, min_time: float = None, max_time: float = None):
        # Initialize figure

        if min_time is None:
            min_time = self.time_data.min()
            index = 0
        else:
            min_time = self.time_data[0][index]
            index = min(np.where(self.time_data[0] > min_time)[0])

        reactant = self.measured_reactant

        if self.measured_reactant_role == SBOTerm.SUBSTRATE:
            vismode = 0
        else:
            vismode = 2

        n_colors = len(self.measurements)
        colors = px.colors.sample_colorscale(
            "turbo", [n / (n_colors - 1) for n in range(n_colors)]
        )

        fig = go.Figure()

        annotations = []
        steps = []

        # Add measured data
        new_colors = []
        show_legend = True
        for measurement, substrate, color in zip(
            self._get_species_data(reactant),
            self._get_species_data(self.substrate),
            colors,
        ):
            if show_legend:
                show_legend = True
            else:
                if old_init_conc != substrate.init_conc:
                    show_legend = True

            for replicate in measurement.replicates:
                if any(np.isnan(replicate.data)):
                    continue
                fig.add_trace(
                    go.Scatter(
                        x=replicate.time,
                        y=replicate.data,
                        mode="markers",
                        customdata=["measured"],
                        name=(
                            f"{substrate.init_conc} {self._format_unit(measurement.unit)}"
                        ),
                        marker=dict(color=color),
                        hovertemplate=f"Well ID: {replicate.id}",
                        showlegend=show_legend,
                    )
                )
                old_init_conc = measurement.init_conc
                show_legend = False
            new_colors.append(color)

        # Add annotation for raw data
        annotations.append(
            go.layout.Annotation(
                font=dict(color="black", size=10),
                x=0,
                y=0,
                showarrow=False,
                text=f"",
                textangle=0,
                xref="x",
                yref="paper",
                xanchor="left",
            )
        )

        # Add simulated data for each model
        systems = [
            system for system in self.reaction_systems if system.result.fit_success
        ]
        systems = sorted(systems, key=lambda system: system.result.AIC)

        dense_time = np.linspace(min_time, max(self.time_data.flatten()), 100)

        for system in systems:
            substrate_data = self._get_species_data(self.substrate)
            enzyme_data = self._get_species_data(self.enzyme)
            product_data = self._get_species_data(self.product)

            for meas_id, (substrate, enzyme, product, color) in enumerate(
                zip(substrate_data, enzyme_data, product_data, colors)
            ):
                init_conditions = np.zeros((1, 3))
                if substrate.replicates:
                    init_conditions[0, 0] = np.nanmean(
                        [rep.data[index] for rep in substrate.replicates]
                    )
                else:
                    if product.replicates:
                        init_conditions[0, 0] = substrate.init_conc - np.nanmean(
                            [rep.data[index] for rep in product.replicates]
                        )
                    else:
                        continue

                # enzyme
                init_conditions[0, 1] = enzyme.init_conc

                # product
                if product.replicates:
                    init_conditions[0, 2] = np.nanmean(
                        [rep.data[index] for rep in product.replicates]
                    )
                else:
                    init_conditions[0, 2] = product.init_conc - np.nanmean(
                        [rep.data[index] for rep in substrate.replicates]
                    )

                simulated_substrates = system.simulate(
                    [dense_time], init_conditions, system.fitted_params_dict
                )[:, :, vismode]

                # Add data traces for each model
                fig.add_trace(
                    go.Scatter(
                        x=dense_time,
                        y=simulated_substrates[0],
                        name=f"{system.name}",
                        mode="lines",
                        marker=dict(color=color),
                        customdata=[f"{system.name}"],
                        hoverinfo="skip",
                        showlegend=False,
                        visible=False,
                    )
                )

            # Add annotations for each model
            label_pos = -0.5
            annotations.append(
                go.layout.Annotation(
                    font=dict(color="black", size=10),
                    x=0,
                    y=label_pos,
                    showarrow=False,
                    text=f"{system._style_parameters()}",
                    textangle=0,
                    xref="paper",
                    yref="paper",
                    xanchor="left",
                )
            )

        # Add step for raw data
        steps.append(
            dict(
                method="update",
                args=[
                    dict(
                        visible=self._visibility_mask(
                            visible_traces=["measured"], fig_data=fig.data
                        )
                    ),
                    {
                        "title.text": (
                            f"{self.name} at"
                            f" {self.temperature} {self.temperature_unit} and pH"
                            f" {self.ph}"
                        )
                    },
                    {"annotations": [annotations[0]]},
                ],
                label=f"-",
            )
        )

        for annotation, system in zip(annotations[1:], systems):
            step = dict(
                method="update",
                args=[
                    dict(
                        visible=self._visibility_mask(
                            visible_traces=["measured", system.name], fig_data=fig.data
                        )
                    ),
                    dict(annotations=[annotation]),
                    {
                        "title.text": (
                            f"{self.name} at"
                            f" {self.temperature} {self.temperature_unit} and pH"
                            f" {self.ph}"
                        )
                    },
                ],
                label=f"{system.name}",
            )

            steps.append(step)

        # Add Slider
        sliders = [
            dict(
                active=0,
                currentvalue=dict(prefix="Model: ", font=dict(color="black")),
                tickcolor="white",
                tickwidth=0,
                font=dict(color="white"),
                pad={"t": 50},
                steps=steps,
            )
        ]

        fig.update_layout(
            title=(
                f"{self.name} at {self.temperature} {self.temperature_unit} and pH"
                f" {self.ph}"
            ),
            sliders=sliders,
            updatemenus=[
                dict(
                    type="buttons",
                    direction="right",
                    x=0.7,
                    y=1.3,
                    showactive=True,
                )
            ],
            template="simple_white",
            hoverlabel_align="right",
            legend_title=f"Initial {self.substrate.name}",
            xaxis_title=f"time / {self.time_unit}",
            yaxis_title=f"{reactant.name} / {self._format_unit(measurement.unit)}",
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False),
        )

        fig.show()

    @staticmethod
    def hex_to_rgba(hex: str, alpha: float = 1) -> str:
        alpha = int(alpha * 255)
        rgb = tuple(int(hex.strip("#")[i : i + 2], 16) for i in (0, 2, 4))
        return f"rgba{rgb + (alpha,)}"

    @staticmethod
    def _format_unit(unit: str) -> str:
        unit = unit.replace(" / l", " L<sup>-1</sup>")
        unit = unit.replace("1 / s", "s<sup>-1</sup>")
        unit = unit.replace("1 / min", "min<sup>-1</sup>")
        unit = unit.replace("umol", "µmol")
        unit = unit.replace("ug", "µg")
        return unit

    def _get_conditions(self, time: float):
        init_conditions = np.zeros((len(self.measurements), 3))

        index = self.time_data[0].tolist().index(time)

        replicates = self._measurement_replicates

        start_slice = 0
        end_slice = 0
        for measurement_id, rep in enumerate(replicates):
            end_slice += rep
            data_substrate = self.substrate_data[start_slice:end_slice]
            data_enzyme = self.enzyme_data[start_slice:end_slice]
            data_product = self.product_data[start_slice:end_slice]
            start_slice += rep

            init_conditions[measurement_id, 0] = np.nanmean(data_substrate[:, index])
            init_conditions[measurement_id, 1] = np.nanmean(data_enzyme[:, index])
            init_conditions[measurement_id, 2] = np.nanmean(data_product[:, index])

        nan_rows_mask = np.isnan(init_conditions).any(axis=1)

        # Use the mask to select rows without NaN values in the first dimension
        filtered_matrix = init_conditions[~nan_rows_mask]

        return filtered_matrix

    @staticmethod
    def _visibility_mask(visible_traces: list, fig_data: list) -> list:
        return [
            any(fig["customdata"][0] == trace for trace in visible_traces)
            for fig in fig_data
        ]
