import sdRDM
import sympy as sp

from typing import List, Optional, Union
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator


from .abstractspecies import AbstractSpecies
from .protein import Protein
from .reactant import Reactant
from .reaction import Reaction
from .reactionelement import ReactionElement
from .measurement import Measurement
from .sboterm import SBOTerm
from .kineticmodel import KineticModel
from .measurementdata import MeasurementData
from .kineticparameter import KineticParameter
from EnzymePynetics.models import K_CAT, K_M, K_IC, K_IU
from EnzymePynetics.ioutils import parse_enzymeml


SPECIES_ROLES = ["substrate", "product", "enzyme", "inhibitor"]
PARAMETERS = [K_CAT, K_M, K_IC, K_IU]


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

    species: List[AbstractSpecies] = Field(
        description="Reactants, Inhibitor, Activators and Catalysts of the reaction",
        default_factory=ListPlus,
        multiple=True,
    )

    reactions: List[Reaction] = Field(
        description="Reaction proceeding in measurements",
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

    def _add_to_species(self, new_species: AbstractSpecies) -> AbstractSpecies:
        """
        This method adds an object of type 'Species' to attribute species

        Args:
            id (str): Unique identifier of the 'Species' object. Defaults to 'None'.
            type (): Type of the species.
            name (): Name of the species. Defaults to None
        """

        if any([species.id == new_species.id for species in self.species]):
            self.species = [new_species if species.id ==
                            new_species.id else species for species in self.species]

            return new_species

        else:
            self.species.append(new_species)

            return new_species

    def add_protein(
            self,
            id: str,
            name: str,
            constant: bool,
            sequence: str,
            **kwargs
    ):

        # define abstract Vessel object
        vessel = self._define_dummy_vessel()

        params = {
            "id": id,
            "name": name,
            "constant": constant,
            "sequence": sequence,
            "vessel_id": vessel.id,
            **kwargs
        }

        return self._add_to_species(Protein(**params))

    def add_reactant(
            self,
            id: str,
            name: str,
            constant: bool,
            **kwargs
    ):

        # define abstract Vessel object
        vessel = self._define_dummy_vessel()

        params = {
            "id": id,
            "name": name,
            "constant": constant,
            "vessel_id": vessel.id,
            **kwargs
        }

        return self._add_to_species(Reactant(**params))

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
        enzyme: Protein = None,
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
                    ontology=SBOTerm.SUBSTRATE
                )
            ]
        else:
            educt_reaction_element = ListPlus()

        if product:
            product_reaction_element = [
                ReactionElement(
                    species_id=product.id,
                    constant=product.constant,
                    ontology=SBOTerm.PRODUCT
                )
            ]
        else:
            product_reaction_element = ListPlus()

        # Add modifiers
        modifiers = []

        if not enzyme and len(self.enzymes) == 1:
            enzyme = self.enzymes[0]

        if enzyme:
            modifiers.append(
                ReactionElement(
                    species_id=enzyme.id,
                    constant=enzyme.constant,
                    ontology=SBOTerm.CATALYST
                )
            )
        else:
            raise ValueError(
                "No protein defined as enzyme. Use 'add_protein' first.")

        if inhibitor:
            modifiers.append(
                ReactionElement(
                    species_id=inhibitor.id,
                    constant=inhibitor.constant,
                    ontology=SBOTerm.INHIBITOR
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
            self.reactions = [new_reaction if reaction.id ==
                              new_reaction.id else reaction for reaction in self.reactions]

            return new_reaction

        else:
            self.reactions.append(new_reaction)

            return new_reaction

    def add_model(
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

        if not parameters:
            self._extract_parameters(equation)

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

    def _extract_parameters(self, equation: str) -> List[KineticParameter]:

        # local_sympy_dict ensures that 'product' in the string expression is treated as
        # a symbol instead of a function
        sympy_dict = {'product': sp.Symbol('product')}
        expr_substrate = sp.parse_expr(equation, sympy_dict)
        free_symbols = list(expr_substrate.free_symbols)
        param_names = [
            symbol.name for symbol in free_symbols if symbol.name.lower() not in SPECIES_ROLES]

        for param_name in param_names:
            if param_name not in [param.name for param in PARAMETERS]:
                raise ValueError(
                    f"Parameter '{param_name}' is not a valid parameter: ({PARAMETERS})"
                )

        parameters = []
        for parameter in re.findall(r"([a-zA-Z]+)", equation):
            if parameter not in parameters:
                parameters.append(parameter)

        # define abstract Vessel object
        vessel = self._define_dummy_vessel()

        # define parameters
        parameters = [
            KineticParameter(
                name=parameter,
                value=1.0,
                unit="",
                constant=False,
                is_global=False,
                vessel_id=vessel.id,
            )
            for parameter in parameters
        ]

        return parameters

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

    @property
    def ph(self):
        if not all([measurement.ph == self.measurements[0].ph for measurement in self.measurements]):
            raise ValueError("Measurements have inconsistent pH values.")
        return self.measurements[0].ph

    @property
    def temperature(self):
        if not all([measurement.temperature == self.measurements[0].temperature for measurement in self.measurements]):
            raise ValueError(
                "Measurements have inconsistent temperature values.")
        return self.measurements[0].temperature

    @property
    def temperature_unit(self):
        if not all([measurement.temperature_unit == self.measurements[0].temperature_unit for measurement in self.measurements]):
            raise ValueError(
                "Measurements have inconsistent temperature unit values.")
        return self.measurements[0].temperature_unit

    @property
    def reactants(self):
        return [species for species in self.species if species.constant == False]

    @property
    def modifiers(self):
        return [species for species in self.species if species.constant == True]

    @property
    def enzymes(self):
        return [species for species in self.species if species.ontology == SBOTerm.CATALYST.value]

    @classmethod
    def from_enzymeml(cls, enzymeml_doc: Union[str, "EnzymeMLDocument"]):

        return parse_enzymeml(cls, enzymeml_doc)
