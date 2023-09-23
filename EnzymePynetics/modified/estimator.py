import copy
from pyexpat import model
import time
import numpy as np
import pandas as pd
import sdRDM
import sympy as sp
import plotly.express as px
from plotly import graph_objects as go
from plotly.graph_objs import Layout

from typing import List, Optional, Union
from pydantic import Field
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature, IDGenerator

from .abstractspecies import AbstractSpecies
from .modelresult import ModelResult
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

    reaction_systems: List[ReactionSystem] = Field(
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
        new_model.pretty_print
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

    def _create_model_combinations(self):
        inactivation_reaction = self._create_inactivation_reaction()

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
        active = Protein(**self.enzyme.to_dict())
        active.name += " (active)"
        active.constant = False

        inactive = Protein(**self.enzyme.to_dict())
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
        # remove nans

        nan_mask = np.isnan(self.substrate_data).any(axis=1)

        substrate_data = self.substrate_data[~nan_mask]
        product_data = self.product_data[~nan_mask]
        enzyme_data = self.enzyme_data[~nan_mask]
        time_data = self.time_data[~nan_mask]

        return [substrate_data, enzyme_data, product_data, time_data]

    def fit_models(self):
        self._create_model_combinations()

        substrate, enzyme, product, time = self._remove_nans()

        for system in self.reaction_systems:
            print(f"Fitting {system.name}")
            system.fit(
                substrate_data=substrate,
                enzyme_data=enzyme,
                product_data=product,
                times=time,
            )

        return self.fit_statistics()

    def fit_statistics(self):
        header = np.array(
            [
                ["Model", ""],
                ["AIC", ""],
                [ParamType.K_CAT.value, f"1 / {self.time_unit}"],
                [ParamType.K_M.value, self.substrate_unit],
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
            ParamType.K_M.value: "{:.2f}",
            ParamType.K_CAT.value: "{:.2f}",
            ParamType.K_IE.value: "{:.2f}",
            ParamType.K_IU.value: "{:.2f}",
            ParamType.K_IC.value: "{:.2f}",
        }

        df = pd.DataFrame(entries).set_index("Model").sort_values("AIC")
        df.columns = pd.MultiIndex.from_arrays(header[1:, :].T)

        return (
            df.style.format("{:.2f}", na_rep="")
            .format("{:.0f}", subset=["AIC"], na_rep="failed")
            .background_gradient(cmap="Blues", subset=["AIC"])
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

        params = {
            "name": name,
            "reactions": reactions,
            "result": result,
        }

        if id is not None:
            params["id"] = id

        self.reaction_systems.append(ReactionSystem(**params))

        return self.reaction_systems[-1]

    def _init_parameters(self, model: KineticModel) -> KineticModel:
        value = float("nan")

        for param in model.eq_parameters:
            if param == ParamType.K_CAT.value:
                ontology = SBOTerm.K_CAT
                initial_value = self._init_kcat
                unit = f"1 / {self.time_unit}"
                upper = initial_value * 1000
                lower = initial_value * 0.001

            elif param == ParamType.K_M.value:
                ontology = SBOTerm.K_M
                initial_value = self._init_km
                unit = self.substrate_unit
                upper = initial_value * 1000
                lower = initial_value * 0.001

            elif param == ParamType.K_IC.value:
                initial_value = self._init_km / 10
                unit = self.substrate_unit
                ontology = None
                upper = initial_value * 1000
                lower = initial_value * 0.001

            elif param == ParamType.K_IU.value:
                initial_value = self._init_km / 10
                unit = self.substrate_unit
                ontology = None
                upper = initial_value * 1000
                lower = initial_value * 0.001

            elif param == ParamType.K_IE.value:
                if self.time_unit == "s":
                    initial_value = np.log(2) / 90 * 60  # 90 min half life
                if self.time_unit == "min":
                    initial_value = np.log(2) / 90  # 90 min half life
                unit = f"1 / {self.time_unit}"
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

    def get_reaction_system(self, system_name: str) -> ReactionSystem:
        for system in self.reaction_systems:
            if system.name == system_name:
                return system

        raise ValueError(f"Reaction system '{system_name}' not found.")

    def visualize(self, reaction_system: ReactionSystem = None):
        # Initialize figure

        reactant = self.measured_reactant

        colors = px.colors.qualitative.Plotly

        fig = go.Figure()

        annotations = []
        steps = []

        # Add measured data
        new_colors = []
        show_legend = True
        for measurement, color in zip(self._get_species_data(reactant), colors):
            if show_legend:
                show_legend = True
            else:
                if old_init_conc != measurement.init_conc:
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
                        name=f"{measurement.init_conc} {self._format_unit(measurement.unit)}",
                        marker=dict(color=self.hex_to_rgba(color)),
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

        dense_time = np.linspace(0, max(self.time_data.flatten()), 100)

        for system in systems:
            init_conditions = self._get_init_conditions

            times = np.tile(dense_time, init_conditions.shape[0]).reshape(
                init_conditions.shape[0], -1
            )

            simulated_substrates = system.simulate(
                times, init_conditions, system.fitted_params_dict
            )[:, :, 0]

            for sim_substrate, time, color in zip(
                simulated_substrates, times, reversed(new_colors)
            ):
                # Add data traces for each model
                fig.add_trace(
                    go.Scatter(
                        x=time,
                        y=sim_substrate,
                        name=f"{system.name}",
                        mode="lines",
                        marker=dict(color="#000000"),
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
                        "title.text": f"{self.name} at {self.temperature} {self.temperature_unit} and pH {self.ph}"
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
                        "title.text": f"{self.name} at {self.temperature} {self.temperature_unit} and pH {self.ph}"
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
            title=f"{self.name} at {self.temperature} {self.temperature_unit} and pH {self.ph}",
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
            legend_title=f"Initial {reactant.name}",
            xaxis_title=f"time / {self.time_unit}",
            yaxis_title=f"{reactant.name} / {self._format_unit(measurement.unit)}",
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False),
        )

        fig.show()

    @staticmethod
    def hex_to_rgba(hex: str) -> str:
        rgb = tuple(int(hex.strip("#")[i : i + 2], 16) for i in (0, 2, 4))
        return f"rgba{rgb + (255,)}"

    @staticmethod
    def _format_unit(unit: str) -> str:
        unit = unit.replace(" / l", " L<sup>-1</sup>")
        unit = unit.replace("1 / s", "s<sup>-1</sup>")
        unit = unit.replace("1 / min", "min<sup>-1</sup>")
        unit = unit.replace("umol", "µmol")
        unit = unit.replace("ug", "µg")
        return unit

    @property
    def _get_init_conditions(self):
        init_conditions = np.zeros((len(self.measurements), 3))

        for i, measurement in enumerate(self.measurements):
            for data in measurement.species:
                if data.species_id == self.substrate.id:
                    if not data.replicates:
                        init_conditions[i, 0] = data.init_conc
                    else:
                        init_conditions[i, 0] = np.mean(
                            [replicate.data[0] for replicate in data.replicates]
                        )

                if data.species_id == self.enzyme.id:
                    if not data.replicates:
                        init_conditions[i, 1] = data.init_conc
                    else:
                        init_conditions[i, 1] = np.mean(
                            [replicate.data[0] for replicate in data.replicates]
                        )

                if data.species_id == self.product.id:
                    if not data.replicates:
                        init_conditions[i, 2] = data.init_conc
                    else:
                        init_conditions[i, 2] = np.mean(
                            [replicate.data[0] for replicate in data.replicates]
                        )

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
