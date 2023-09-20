import copy
from pyexpat import model
import time
import numpy as np
import sdRDM
import sympy as sp

from typing import List, Optional, Union
from pydantic import Field, PrivateAttr
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from .abstractspecies import AbstractSpecies
from .protein import Protein
from .reactant import Reactant
from .reaction import Reaction
from .reactionelement import ReactionElement
from .reactionsystem import ReactionSystem
from .measurement import Measurement
from .sboterm import SBOTerm, ParamType
from .kineticmodel import KineticModel
from .measurementdata import MeasurementData
from .kineticparameter import KineticParameter
from EnzymePynetics.ioutils import parse_enzymeml


SPECIES_ROLES = ["substrate", "product", "enzyme", "inhibitor"]


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

    measured_reactant: Reactant = Field(
        default=None,
        description="Reactant that is measured in the experiment",
    )

    reaction_systems: Optional[ReactionSystem] = Field(
        description="Reactions of multiple species",
        default_factory=ListPlus,
        multiple=True,
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

        if not enzyme and len(self.enzymes) == 1:
            enzyme = self.enzymes[0]

        if enzyme:
            modifiers.append(
                ReactionElement(
                    species_id=enzyme.id,
                    constant=enzyme.constant,
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

    def _validate_units(self):
        if not self.substrate_unit == self.product_unit:
            raise ValueError("Substrate and product have different units.")

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
        parameterized_model = self._init_parameters(new_model)

        if any([model.id == new_model.id for model in self.models]):
            self.models = [
                new_model if model.id == new_model.id else model
                for model in self.models
            ]

            return new_model

        else:
            self.models.append(new_model)

            return new_model

    def _init_parameters(self, model: KineticModel) -> KineticModel:
        value = float("nan")

        for param in model.eq_parameters:
            if param == ParamType.K_CAT.value:
                ontology = SBOTerm.K_CAT
                initial_value = self._init_kcat
                unit = f"1 / {self.time_unit}"
                upper = initial_value * 100
                lower = initial_value * 0.001

            elif param == ParamType.K_M.value:
                ontology = SBOTerm.K_M
                initial_value = self._init_km
                unit = self.substrate_unit
                upper = initial_value * 100
                lower = initial_value * 0.001

            elif param == ParamType.K_IC.value:
                initial_value = self._init_km
                unit = self.substrate_unit
                ontology = None
                upper = initial_value * 100
                lower = initial_value * 0.001

            elif param == ParamType.K_IU.value:
                initial_value = self._init_km
                unit = self.substrate_unit
                ontology = None
                upper = initial_value * 100
                lower = initial_value * 0.001

            elif param == ParamType.K_IE.value:
                if self.time_unit == "s":
                    initial_value = np.log(2) / 90 * 60  # 90 min half life
                if self.time_unit == "min":
                    initial_value = np.log(2) / 90  # 90 min half life
                unit = self.time_unit
                ontology = None
                upper = initial_value * 100
                lower = initial_value * 0.001

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
        return self._get_consisten_unit(self.substrate)

    @property
    def product_unit(self):
        return self._get_consisten_unit(self.product)

    @property
    def enzyme_unit(self):
        return self._get_consisten_unit(self.enzyme)

    @property
    def inhibitor_unit(self):
        return self._get_consisten_unit(self.inhibitor)

    def _get_consisten_unit(self, species: AbstractSpecies) -> None:
        units = [measurement.unit for measurement in self._get_species_data(species)]
        if not all([unit == units[0] for unit in units]):
            raise ValueError("Measurements have inconsistent substrate units.")

        return units[0]

    @property
    def reactants(self):
        return [species for species in self.species if species.constant == False]

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
            model
            for model in self.models
            if model.equation.startswith(SPECIES_ROLES[0])
        ]

    @property
    def enzyme_models(self):
        return [
            model
            for model in self.models
            if model.equation.startswith(SPECIES_ROLES[2])
        ]

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
            f"Measured reactant '{self.measured_reactant}' not found in the defined reaction."
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

        return np.array(init_substrates).flatten()

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
        for n_reps, measurement in zip(self._measurement_replicates, self.measurements):
            time_data.append([measurement.global_time] * n_reps)

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

        return np.repeat(np.array(enzyme_data), self.time_data.shape[1]).reshape(
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
        for measurment in self.measurements:
            for data in measurment.species:
                if data.species_id == species.id:
                    yield data

    @classmethod
    def from_enzymeml(
        cls,
        enzymeml_doc: Union[str, "EnzymeMLDocument"],
        measured_reactant: Union[Reactant, str],
    ):
        return parse_enzymeml(cls, enzymeml_doc, measured_reactant)
