from typing import List, Dict, Optional, Literal, Union

from pyenzyme import EnzymeMLDocument
from EnzymePynetics.core.enzymekinetics import EnzymeKinetics
from EnzymePynetics.core.speciestypes import SpeciesTypes
from EnzymePynetics.core.series import Series
from EnzymePynetics.core.measurement import Measurement
from EnzymePynetics.core.species import Species
from EnzymePynetics.tools.kineticmodel import KineticModel
from EnzymePynetics.tools.rate_equations import *

import numpy as np
import pandas as pd
from pandas import DataFrame
from pathlib import Path
from scipy.integrate import odeint
from lmfit import report_fit
from IPython.display import display
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import matplotlib.colors
from matplotlib.pyplot import cm

_SPECIES_TYPES = Literal["substrate", "product"]


class ParameterEstimator:
    def __init__(self, data: EnzymeKinetics):
        self.data = data
        self.models: Dict[str, KineticModel] = None
        self._measured_species = None

        (
            self.substrate,
            self.initial_substrate,
            self.product,
            self.enzyme,
            self.inhibitor,
            self.time,
        ) = self._initialize_measurement_data()
        self.initial_kcat = self._calculate_kcat()
        self.initial_Km = self._calculate_Km()

        # TODO shapcheck function to check for consistent array lengths

    def fit_models(
        self,
        initial_substrate_concs: list = None,
        start_time_index: int = None,
        stop_time_index: int = None,
        only_irrev_MM: bool = False,
        display_output: bool = True,
    ) -> None:
        """Fits the measurement data to a set of kinetic models.

        Args:
            initial_substrate_concs (list, optional): Enables to subset the measurement data by choosing one ore multiple initial substrate concentrations. Defaults to None.
            start_time_index (int, optional): Subset the data by choosing the index of thee first measurement point which should be condidered. Defaults to None.
            stop_time_index (int, optional): Choose last measurement point which should be considered. Defaults to None.
            enzyme_inactivation (bool, optional): _description_. Defaults to False.

        Prints DataFrame of all kinetic parameters of all fitted models sorted by Akaike information criterion.
        """

        # Subset data if one or multiple attributes are passed to the function
        if np.any([initial_substrate_concs, start_time_index, stop_time_index]):
            (
                self.subset_substrate,
                self.subset_product,
                self.subset_enzyme,
                self.subset_initial_substrate,
                self.subset_time,
                self.subset_inhibitor,
            ) = self._subset_data(
                initial_substrates=initial_substrate_concs,
                start_time_index=start_time_index,
                stop_time_index=stop_time_index,
            )
        else:
            self.subset_substrate = self.substrate
            self.subset_product = self.product
            self.subset_enzyme = self.enzyme
            self.subset_initial_substrate = self.initial_substrate
            self.subset_time = self.time
            self.subset_inhibitor = self.inhibitor

        # Initialize kinetics models
        self.models = self._initialize_models(
            substrate=self.subset_substrate,
            product=self.subset_product,
            enzyme=self.subset_enzyme,
            inhibitor=self.subset_inhibitor,
            only_irrev_MM=only_irrev_MM,
        )

        self._run_minimization(display_output)

        # Set units in ModelResults object
        for model_name, model in self.models.items():
            for p, parameter in enumerate(model.result.parameters):
                if parameter.name == "k_cat" or parameter.name == "K_ie":
                    self.models[model_name].result.parameters[
                        p
                    ].unit = f"1 / {self._time_unit}"
                elif parameter.name == "Km":
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self._substrate_unit
                elif np.any(self.subset_inhibitor != 0):
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self._inhibitor_unit
                else:
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self._substrate_unit

        self.result_dict = self._result_overview()
        if display_output:
            display(self.result_dict)

    def get_model_results(self, model: str = None):
        if model == None:
            model = self.result_dict.index[0]

        return self.models[model].result.params

    def _initialize_measurement_data(self):
        """
        Extracts data from data objects and reshapes it for fitting.
        """

        substrate_data = []
        initial_substrate_data = []
        product_data = []
        enzyme_data = []
        inhibitor_data = []
        time_data = []

        # extract data from each measurement and species
        for measurement in self.data.measurements:
            for species in measurement.species:
                if species.species_type == SpeciesTypes.SUBSTRATE.value:
                    initial_substrate_data.append(species.initial_conc)
                    self._substrate_unit = species.conc_unit
                    if len(species.data) != 0:
                        self._measured_species = species.species_type
                        for replicate in species.data:
                            substrate_data.append(replicate.values)
                            self._time_unit = replicate.time_unit
                            time_data.append(replicate.time)

                if species.species_type == SpeciesTypes.PRODUCT.value:
                    self._product_unit = species.conc_unit
                    if len(species.data) != 0:
                        self._measured_species = species.species_type
                        for replicate in species.data:
                            product_data.append(replicate.values)
                            self._time_unit = replicate.time_unit
                            time_data.append(replicate.time)

                if species.species_type == SpeciesTypes.INHIBITOR.value:
                    self._inhibitor_unit = species.conc_unit
                    if len(species.data) != 0:
                        for replicate in species.data:
                            inhibitor_data.append(replicate.values)
                    else:
                        inhibitor_data.append(species.initial_conc)

                if species.species_type == SpeciesTypes.ENZYME.value:
                    self._enzyme_unit = species.conc_unit
                    if len(species.data) != 0:
                        for replicate in species.data:
                            enzyme_data.append(replicate.values)
                    else:
                        enzyme_data.append(species.initial_conc)

        # np arrays
        substrate_array = np.array(substrate_data)
        initial_substrate_array = np.array(initial_substrate_data)
        product_array = np.array(product_data)
        enzyme_array = np.array(enzyme_data)
        inhibitor_array = np.array(inhibitor_data)
        time_array = np.array(time_data)

        # get shape information of data
        n_measurements = len(self.data.measurements)
        data_shape = (
            substrate_array.shape if substrate_array.size > 0 else product_array.shape
        )
        n_replicates = int(data_shape[0] / n_measurements)

        # adjust data arrays for inhibitor and enzyme according to number of replicates of each measurement
        enzyme_array = np.repeat(enzyme_array, n_replicates)
        enzyme_array = np.repeat(enzyme_array, data_shape[1]).reshape(data_shape)
        initial_substrate_array = np.repeat(initial_substrate_array, n_replicates)
        if len(inhibitor_array) == 0:
            inhibitor_array = np.zeros(data_shape)
        else:
            inhibitor_array = np.repeat(inhibitor_array, n_replicates)
            inhibitor_array = np.repeat(inhibitor_array, data_shape[1]).reshape(
                data_shape
            )

        # calcualte missing species based on specified initial substrate concentration
        if len(product_data) == 0:
            product_array = np.array(
                self._calculate_product(substrate_data, initial_substrate_array)
            )
        if len(substrate_data) == 0:
            substrate_array = np.array(
                self._calculate_substrate(product_data, initial_substrate_array)
            )
            self._substrate_unit = self._product_unit

        return (
            substrate_array,
            initial_substrate_array,
            product_array,
            enzyme_array,
            inhibitor_array,
            time_array,
        )

    def _calculate_substrate(
        self, product_data: List[List], initial_substrates: List
    ) -> List[List]:
        substrate = []

        for product_measurment, initial_substrate in zip(
            product_data, initial_substrates
        ):
            substrate.append(
                [initial_substrate - value for value in product_measurment]
            )

        return substrate

    def _calculate_product(
        self, substrate_data: List[List], initial_substrates: List
    ) -> List[List]:
        product = []
        for substrate_measurement, initial_substrate in zip(
            substrate_data, initial_substrates
        ):
            product.append(
                [initial_substrate - value for value in substrate_measurement]
            )

        return product

    @staticmethod
    def _get_y0s(substrate, enzyme, product, inhibitor):
        y0s = []
        for s, e, p, i in zip(substrate, enzyme, product, inhibitor):
            y0s.append((s[0], e[0], p[0], i[0]))

        return np.array(y0s)

    def _calculate_rates(self):
        """
        Calculates the change per time-unit between all measurement points.
        """
        concentration_intervals = np.diff(self.substrate)
        time_intervals = np.diff(self.time)
        rates = abs(concentration_intervals / time_intervals)
        return rates

    def _calculate_kcat(self) -> float:
        rates = self._calculate_rates()
        kcat = np.nanmax(rates.T / self.enzyme[:, 0])
        return kcat

    def _calculate_Km(self):
        return np.nanmax(self._calculate_rates() / 2)

    def _subset_data(
        self,
        initial_substrates: list = None,
        start_time_index: int = None,
        stop_time_index: int = None,
    ) -> tuple:
        """This function allows to subset the actual measurement data. Thereby, measurements of specific initial substrate concentrations
        can be specified. Additionally, the time-frame can be specified by defining the index of the first and last measurement time-point.

        Args:
            initial_substrates (list, optional). Defaults to None.
            start_time_index (int, optional). Defaults to None.
            stop_time_index (int, optional). Defaults to None.

        Raises:
            ValueError: If concentrations are passed, which are not defined in "initial_substrate_concentration"

        Returns:
            tuple: Data subset.
        """
        idx = np.array([])
        if initial_substrates == None or len(initial_substrates) == 0:
            idx = np.arange(self.substrate.shape[0])
        else:
            for concentration in initial_substrates:
                if concentration not in self.initial_substrate:
                    raise ValueError(
                        f"{concentration} not found in initial substrate concentrations. \nInitial substrate concentrations are {list(np.unique(self.initial_substrate))}"
                    )
                else:
                    idx = np.append(
                        idx, np.where(self.initial_substrate == concentration)[0]
                    )
        idx = idx.astype(int)

        new_substrate = self.substrate[idx, start_time_index:stop_time_index]
        new_product = self.product[idx, start_time_index:stop_time_index]
        new_enzyme = self.enzyme[idx, start_time_index:stop_time_index]
        new_initial_substrate = self.initial_substrate[idx]
        new_time = self.time[start_time_index:stop_time_index]
        new_inhibitor = self.inhibitor[idx]

        return (
            new_substrate,
            new_product,
            new_enzyme,
            new_initial_substrate,
            new_time,
            new_inhibitor,
        )

    def get_parameter_dict(self):
        return self.models[self.result_dict.index[0]].result.params

    def _initialize_models(
        self, substrate, product, enzyme, inhibitor, only_irrev_MM: bool
    ) -> Dict[str, KineticModel]:
        """Initializer for all kinetic models. If an inhibitor is provided, inhibition models are initialized. If no inhibitor is specified,
        inhibitory models for substrate and product inhibition are initialized additionally to the irreversible Michaelis Menten model.
        """

        y0 = self._get_y0s(
            substrate=substrate, product=product, enzyme=enzyme, inhibitor=inhibitor
        )

        irreversible_Michaelis_Menten_inactivation = KineticModel(
            name="irreversible Michaelis Menten with enzyme inactivation",
            params=[],
            y0=y0,
            kcat_initial=self.initial_kcat,
            Km_initial=self.initial_Km,
            model=irreversible_model,
            enzyme_inactivation=True,
        )

        irreversible_Michaelis_Menten = KineticModel(
            name="irreversible Michaelis Menten",
            params=[],
            y0=y0,
            kcat_initial=self.initial_kcat,
            Km_initial=self.initial_Km,
            model=irreversible_model,
            enzyme_inactivation=False,
        )

        if only_irrev_MM:
            return {
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
            }

        if np.all(self.inhibitor == 0):
            ### Product inhibition models

            y0 = self._get_y0s(
                substrate=self.substrate,
                enzyme=self.enzyme,
                product=self.product,
                inhibitor=self.product,
            )

            competitive_product_inhibition = KineticModel(
                name="competitive product inhibition",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            competitive_product_inhibition_inactivation = KineticModel(
                name="competitive product inhibition with enzyme inactivation",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            uncompetitive_product_inhibition = KineticModel(
                name="uncompetitive product inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            uncompetitive_product_inhibition_inactivation = KineticModel(
                name="uncompetitive product inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            noncompetitive_product_inhibition = KineticModel(
                name="non-competitive product inhibition",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            noncompetitive_product_inhibition_inactivation = KineticModel(
                name="non-competitive product inhibition with enzyme inactivation",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            ### Substrate inhibition models

            y0 = self._get_y0s(
                substrate=self.substrate,
                enzyme=self.enzyme,
                product=self.product,
                inhibitor=self.substrate,
            )

            substrate_inhibition = KineticModel(
                name="substrate inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=substrate_inhibition_model,
                enzyme_inactivation=False,
            )
            substrate_inhibition_inactivation = KineticModel(
                name="substrate inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=substrate_inhibition_model,
                enzyme_inactivation=True,
            )

            model_dict = {
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
                competitive_product_inhibition.name: competitive_product_inhibition,
                competitive_product_inhibition_inactivation.name: competitive_product_inhibition_inactivation,
                uncompetitive_product_inhibition.name: uncompetitive_product_inhibition,
                uncompetitive_product_inhibition_inactivation.name: uncompetitive_product_inhibition_inactivation,
                uncompetitive_product_inhibition.name: uncompetitive_product_inhibition,
                uncompetitive_product_inhibition_inactivation.name: uncompetitive_product_inhibition_inactivation,
                noncompetitive_product_inhibition.name: noncompetitive_product_inhibition,
                noncompetitive_product_inhibition_inactivation.name: noncompetitive_product_inhibition_inactivation,
                substrate_inhibition.name: substrate_inhibition,
                substrate_inhibition_inactivation.name: substrate_inhibition_inactivation,
            }

            return model_dict

        else:
            y0 = self._get_y0s(
                self.substrate, self.enzyme, self.product, self.inhibitor
            )

            competitive_inhibition = KineticModel(
                name="competitive inhibition",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_model,
                enzyme_inactivation=False,
            )
            competitive_inhibition_inactivation = KineticModel(
                name="competitive inhibition with enzyme inactivation",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_model,
                enzyme_inactivation=True,
            )

            uncompetitive_inhibition = KineticModel(
                name="uncompetitive inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_inhibition_model,
                enzyme_inactivation=False,
            )
            uncompetitive_inhibition_inactivation = KineticModel(
                name="uncompetitive inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_inhibition_model,
                enzyme_inactivation=True,
            )

            noncompetitive_inhibition = KineticModel(
                name="non-competitive inhibition",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_inhibition_model,
                enzyme_inactivation=False,
            )
            noncompetitive_inhibition_inactivation = KineticModel(
                name="non-competitive inhibition with enzyme inactivation",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_inhibition_model,
                enzyme_inactivation=True,
            )

            partially_competitive_inhibition = KineticModel(
                name="partially competitive inhibition",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=partially_competitive_inhibition_model,
                enzyme_inactivation=False,
            )
            partially_competitive_inhibition_inactivation = KineticModel(
                name="partially competitive inhibition with enzyme inactivation",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=partially_competitive_inhibition_model,
                enzyme_inactivation=True,
            )

            return {
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
                competitive_inhibition.name: competitive_inhibition,
                competitive_inhibition_inactivation.name: competitive_inhibition_inactivation,
                uncompetitive_inhibition.name: uncompetitive_inhibition,
                uncompetitive_inhibition_inactivation.name: uncompetitive_inhibition_inactivation,
                noncompetitive_inhibition.name: noncompetitive_inhibition,
                noncompetitive_inhibition_inactivation.name: noncompetitive_inhibition_inactivation,
                partially_competitive_inhibition.name: partially_competitive_inhibition,
                partially_competitive_inhibition_inactivation.name: partially_competitive_inhibition_inactivation,
            }

    def _run_minimization(self, display_output: bool):
        """Performs non-linear least-squared minimization to fit the data to the kinetic
        models by adjusting the kinetic parameters of the models.

        Returns:
            DataFrame: Overview of the kinetic parameters of all kinetic models.
        """

        if display_output:
            print("Fitting data to:")
        for kineticmodel in self.models.values():
            if display_output:
                print(f" - {kineticmodel.name} model")

            kineticmodel.fit(self.subset_substrate, self.subset_time)
            if display_output:
                print(f" -- Fitting succeeded: {kineticmodel.result.fit_success}")

    def _result_overview(self) -> DataFrame:
        """
        Prettifies the results of all kinetic models and organized them in a pandas DataFrame.
        """

        if np.all(self.inhibitor == 0):
            inhibitor_unit = self._substrate_unit
        else:
            inhibitor_unit = self._inhibitor_unit

        parameter_mapper = {
            "k_cat": f"kcat [1/{self._time_unit}]",
            "Km": f"Km [{self._substrate_unit}]",
            "K_ic": f"Ki competitive [{inhibitor_unit}]",
            "K_iu": f"Ki uncompetitive [{inhibitor_unit}]",
            "K_ie": f"ki time-dep enzyme-inactiv. [1/{self._time_unit}]",
        }

        result_dict = {}
        for model in self.models.values():
            name = model.name
            if model.result.fit_success:
                aic = round(model.result.AIC)
                rmsd = model.result.RMSD

                parameter_dict = {}
                for parameter in model.result.parameters:
                    name = parameter_mapper[parameter.name]
                    value = parameter.value
                    stderr = parameter.standard_deviation

                    try:
                        percentual_stderr = stderr / value * 100
                    except TypeError:
                        percentual_stderr = float("nan")

                    if name.startswith("Ki time-dep"):
                        parameter_dict[
                            name
                        ] = f"{value:.4f} +/- {percentual_stderr:.2f}%"
                    else:
                        parameter_dict[
                            name
                        ] = f"{value:.4f} +/- {percentual_stderr:.2f}%"

                    if parameter.name == "k_cat":
                        kcat = parameter.value
                        kcat_stderr = parameter.standard_deviation
                    if parameter.name == "Km":
                        Km = parameter.value
                        Km_stderr = parameter.standard_deviation

                if Km_stderr is None or kcat_stderr is None:
                    kcat_Km_stderr = float("nan")
                    kcat_Km = float("nan")
                    percentual_kcat_Km_stderr = float("nan")
                else:
                    kcat_Km = kcat / Km
                    kcat_Km_stderr = (
                        (kcat_stderr / kcat) ** 2 + (Km_stderr / Km) ** 2
                    ) ** 0.5 * kcat_Km
                    percentual_kcat_Km_stderr = kcat_Km_stderr / kcat_Km * 100

                parameter_dict[
                    f"kcat / Km [1/{self._time_unit} * 1/{self._substrate_unit}]"
                ] = f"{kcat_Km:.3f} +/- {percentual_kcat_Km_stderr:.2f}%"

                result_dict[model.name] = {"AIC": aic, "RMSD": rmsd, **parameter_dict}

        df = DataFrame.from_dict(result_dict).T.sort_values("AIC", ascending=True)
        df.fillna("-", inplace=True)
        return df.style.background_gradient(cmap="Blues")

    def get_species(self, species_type: SpeciesTypes, measurement: int = 0) -> Species:
        return next(
            (
                x
                for x in self.data.measurements[measurement].species
                if x.species_type == species_type
            )
        )

    @staticmethod
    def _HEX_to_RGBA_string(color: list) -> str:
        return f"rgba{tuple(color)}"

    @staticmethod
    def _visibility_mask(visible_traces: list, fig_data: list) -> list:
        return [
            any(fig["customdata"][0] == trace for trace in visible_traces)
            for fig in fig_data
        ]

    def _restore_replicate_dims(
        self, initial_substrate: List[float], data_array: np.ndarray
    ):
        unique_initial_substrate = np.unique(initial_substrate)
        data_shape = data_array.shape
        n_replicates = int(data_shape[0] / len(unique_initial_substrate))

        if np.all(self.subset_inhibitor == 0):
            return data_array.reshape(len(unique_initial_substrate), n_replicates, -1)

        else:
            unique_inhibitor_concs = np.unique(self.subset_inhibitor)
            n_replicates = int(
                data_shape[0]
                / len(unique_initial_substrate)
                / len(unique_inhibitor_concs)
            )

            return data_array.reshape(
                len(unique_inhibitor_concs),
                len(unique_initial_substrate),
                n_replicates,
                -1,
            )

    def visualize_subplots(self, visualized_species: _SPECIES_TYPES = None):
        colors = matplotlib.colors.to_rgba_array(px.colors.qualitative.Plotly)

        # Select which species to plot
        if visualized_species is None:
            visualized_species = self._measured_species

        if visualized_species == "substrate":
            measurement_data = self.subset_substrate
            species_tuple = 0
            visualized_species = self._measured_species
        else:
            measurement_data = self.subset_product
            species_tuple = 2

        unique_initial_substrates = np.unique(self.subset_initial_substrate)
        unique_inhibitors = np.unique(self.inhibitor)

        # Reshape data and time array to contain inhibitor and replicate information
        datas = self._restore_replicate_dims(
            self.subset_initial_substrate, measurement_data
        )
        times = self._restore_replicate_dims(
            self.subset_initial_substrate, self.subset_time
        )

        # get measured species and substrate species
        species = self.get_species(visualized_species)
        substrate = self.get_species("substrate")
        inhibitor = self.get_species("inhibitor")

        subplot_titles = []
        for row, inhibitor_conc in enumerate(unique_inhibitors):
            if inhibitor_conc == 0:
                subplot_titles.append(f"without {inhibitor.name}")
            else:
                subplot_titles.append(
                    f"{inhibitor_conc} {inhibitor.conc_unit} {inhibitor.name}"
                )

        fig = make_subplots(
            rows=datas.shape[0],
            cols=1,
            shared_yaxes="all",
            y_title=f"Initial {substrate.name} ({self._substrate_unit})",
            x_title=f"time ({self._time_unit})",
            subplot_titles=subplot_titles,
            horizontal_spacing=0.05,
            vertical_spacing=0.05,
        )

        # Adjust syle of subplot label
        for count, annotation in enumerate(fig.layout["annotations"]):
            if inhibitor.name in annotation["text"]:
                fig.layout["annotations"][count]["x"] = 0
                fig.layout["annotations"][count]["font"]["size"] = 12
                fig.layout["annotations"][count]["xanchor"] = "left"

        for inhibitor_count, (inhibitor_data, inhibitor_time) in enumerate(
            zip(datas, times)
        ):
            for data, time, init_sub, color in zip(
                inhibitor_data, inhibitor_time, unique_initial_substrates, colors
            ):
                for replicate_count, (replicate_data, replicate_time) in enumerate(
                    zip(data, time)
                ):
                    show_legend = (
                        True if replicate_count == 0 and inhibitor_count == 0 else False
                    )
                    color[-1] = 1
                    fig.add_trace(
                        go.Scatter(
                            x=replicate_time,
                            y=replicate_data,
                            name=f"{init_sub}",
                            mode="markers",
                            marker=dict(color=self._HEX_to_RGBA_string(color)),
                            showlegend=show_legend,
                            customdata=["replicates"],
                            hoverinfo="skip",
                            visible=False,
                        ),
                        col=1,
                        row=inhibitor_count + 1,
                    )

                if datas.shape[2] > 1:
                    show_legend = True if inhibitor_count == 0 else False

                    mean = np.mean(data, axis=0)
                    std = np.std(data, axis=0)
                    color[-1] = 1
                    fig.add_trace(
                        go.Scatter(
                            x=time[0],
                            y=mean,
                            name=f"{init_sub}",
                            mode="markers",
                            marker=dict(color=self._HEX_to_RGBA_string(color)),
                            customdata=["mean"],
                            hoverinfo="skip",
                            showlegend=show_legend,
                            visible=True,
                        ),
                        col=1,
                        row=inhibitor_count + 1,
                    )

                    color[-1] = 0.25  # change opacity of color
                    fig.add_trace(
                        go.Scatter(
                            name=f"{init_sub}",
                            x=time[0],
                            y=mean + std,
                            mode="lines",
                            line=dict(width=0),
                            showlegend=False,
                            customdata=["std"],
                            visible=True,
                        ),
                        col=1,
                        row=inhibitor_count + 1,
                    )
                    fig.add_trace(
                        go.Scatter(
                            name=f"{init_sub}",
                            x=time[0],
                            y=mean - std,
                            line=dict(width=0),
                            mode="lines",
                            fillcolor=self._HEX_to_RGBA_string(color),
                            fill="tonexty",
                            showlegend=False,
                            customdata=["std"],
                            visible=True,
                        ),
                        col=1,
                        row=inhibitor_count + 1,
                    )

        # Integrate each successfully fitted model
        successfull_models = []
        steps = []
        for model in self.models.values():
            if model.result.fit_success:
                successfull_models.append(model)
                means_y0s = np.mean(
                    self._restore_replicate_dims(
                        self.subset_initial_substrate, model.y0
                    ),
                    axis=2,
                )
                for inhibitor_count, (inhibitor_y0s, time) in enumerate(
                    zip(means_y0s, times)
                ):
                    datas = model.integrate(
                        model._fit_result.params,
                        time[:, 0, :],
                        inhibitor_y0s,
                    )

                    for data, t, color, init_sub in zip(
                        datas[:, :, species_tuple],
                        time[:, 0, :],
                        colors,
                        unique_initial_substrates,
                    ):
                        color[-1] = 1
                        fig.add_trace(
                            go.Scatter(
                                x=t,
                                y=data,
                                name=model.name,
                                marker=dict(color=self._HEX_to_RGBA_string(color)),
                                customdata=[model.name],
                                showlegend=False,
                                hoverinfo="name",
                                visible=False,
                            ),
                            col=1,
                            row=inhibitor_count + 1,
                        )
        for model in successfull_models:
            steps.append(
                dict(
                    method="update",
                    args=[
                        {
                            "visible": self._visibility_mask(
                                visible_traces=["mean", "std", model.name],
                                fig_data=fig.data,
                            )
                        },
                        {"title": f"Data + Model"},
                    ],
                    label=f"{model.name}",
                )
            )

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

        max_data = np.nanmax([ys["y"] for ys in fig.__dict__["_data_objs"]])

        # # Buttons
        # buttons = []
        # print(self._visibility_mask(visible_traces=["mean", "std"], fig_data=fig.data))
        # buttons.append(
        #     dict(
        #         method="restyle",
        #         label="Averages",
        #         visible=True,
        #         args=[
        #             dict(
        #                 visible=self._visibility_mask(
        #                     visible_traces=["mean", "std"], fig_data=fig.data
        #                 )
        #             )
        #         ],
        #         args2=[
        #             dict(
        #                 visible=self._visibility_mask(
        #                     visible_traces=["replicates"], fig_data=fig.data
        #                 )
        #             )
        #         ],
        #     )
        # )

        fig.update_layout(
            sliders=sliders,
            updatemenus=[
                dict(
                    type="buttons",
                    direction="right",
                    x=0.7,
                    y=1.3,
                    showactive=True,
                    # buttons=buttons,
                )
            ],
            yaxis_range=[0 - 0.05 * max_data, max_data + 0.05 * max_data],
        )

        # Add title, legend...
        fig.update_layout(
            showlegend=True,
            title="Measured data",
            yaxis_title=f"{species.name} ({species.conc_unit})",
            xaxis_title=f"time ({self._time_unit})",
            hovermode="closest",
            legend_title_text=f"Initial {substrate.name} ({self._substrate_unit})",
            hoverlabel_namelength=-1,
        )

        fig.update_layout(height=400 + len(unique_inhibitors) * 150)

        return fig

    def visualize_model_overview(self, visualized_species: _SPECIES_TYPES = None):
        fig = go.Figure()
        colors = matplotlib.colors.to_rgba_array(px.colors.qualitative.T10_r)

        # Select which species to plot
        if visualized_species is None:
            visualized_species = self._measured_species

        if visualized_species == "substrate":
            measurement_data = self.subset_substrate
            species_tuple = 0
            visualized_species = self._measured_species
        else:
            measurement_data = self.subset_product
            species_tuple = 2

        datas = self._restore_replicate_dims(
            self.subset_initial_substrate, measurement_data
        )
        times = self._restore_replicate_dims(
            self.subset_initial_substrate, self.subset_time
        )
        initial_substrates = []
        for sub in self.subset_initial_substrate:
            if sub not in initial_substrates:
                initial_substrates.append(sub)

        # plot replicates

        for i, (data, time, init_sub) in enumerate(
            zip(datas, times, initial_substrates)
        ):
            for j, (replicate_data, replicate_time) in enumerate(zip(data, time)):
                show_legend = True if j == 0 else False
                fig.add_trace(
                    go.Scatter(
                        x=replicate_time,
                        y=replicate_data,
                        name=f"{init_sub}",
                        mode="markers",
                        marker=dict(color=self._HEX_to_RGBA_string(colors[i])),
                        showlegend=show_legend,
                        customdata=["replicates"],
                        hoverinfo="skip",
                        visible=False,
                    )
                )

        # plot means
        if len(datas.shape) == 3:
            means = np.mean(datas, axis=1)
            stds = np.std(datas, axis=1)

            for color, mean, std, init_sub, time in zip(
                colors, means, stds, initial_substrates, times[:, 0, :]
            ):
                fig.add_trace(
                    go.Scatter(
                        x=time,
                        y=mean,
                        name=f"{init_sub}",
                        mode="markers",
                        marker=dict(color=self._HEX_to_RGBA_string(color)),
                        customdata=["mean"],
                        hoverinfo="skip",
                        visible=True,
                    )
                )

                color[-1] = 0.25  # change opacity of color
                fig.add_trace(
                    go.Scatter(
                        name=f"{init_sub}",
                        x=time,
                        y=mean + std,
                        mode="lines",
                        line=dict(width=0),
                        showlegend=False,
                        customdata=["std"],
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        name=f"{init_sub}",
                        x=time,
                        y=mean - std,
                        line=dict(width=0),
                        mode="lines",
                        fillcolor=self._HEX_to_RGBA_string(color),
                        fill="tonexty",
                        showlegend=False,
                        customdata=["std"],
                    )
                )

        # Integrate each successfully fitted model
        successful_models = []
        steps = []
        annotations = []

        for model_name in self.result_dict.index:
            model = self.models[model_name]
            if model.result.fit_success:
                successful_models.append(model.name)
                mean_y0s = np.mean(
                    self._restore_replicate_dims(
                        self.subset_initial_substrate, model.y0
                    ),
                    axis=1,
                )
                datas = model.integrate(
                    model._fit_result.params,
                    self.subset_time,
                    mean_y0s,
                )

                for data, time, color, init_sub in zip(
                    datas[:, :, species_tuple],
                    times[:, 0, :],
                    colors,
                    initial_substrates,
                ):
                    color[-1] = 1
                    fig.add_trace(
                        go.Scatter(
                            x=time,
                            y=data,
                            mode="lines",
                            name=model.name,
                            marker=dict(color=self._HEX_to_RGBA_string(color)),
                            customdata=[f"{model.name}"],
                            showlegend=False,
                            hoverinfo="name",
                            visible=False,
                        )
                    )
                param_value_map = dict(
                    k_cat="<b><i>k</i><sub>cat</sub>:</b>",
                    Km="<b><i>K</i><sub>m</sub>:</b>",
                    K_ie="<b><i>K</i><sub>ie</sub>:</b>",
                    K_ic="<b><i>K</i><sub>ic</sub>:</b>",
                    K_iu="<b><i>K</i><sub>iu</sub>:</b>",
                )
                params = ""
                for parameter in model.result.parameters:
                    params = (
                        params
                        + f"{param_value_map[parameter.name]} "
                        + f"{parameter.value:.3f} "
                        + f"{parameter.unit} \n\n"
                    )
                annotations.append(
                    go.layout.Annotation(
                        font=dict(color="black"),
                        x=0,
                        y=-0.55,
                        showarrow=False,
                        text=f"<b>AIC:</b> {round(model.result.AIC)}\n\n {params}",
                        textangle=0,
                        xref="x",
                        yref="paper",
                        xanchor="left",
                    )
                )

        empty_annotation = go.layout.Annotation(
            font=dict(color="white"),
            x=0,
            y=-0.55,
            showarrow=False,
            text=f"",
            textangle=0,
            xref="x",
            yref="paper",
            xanchor="left",
        )

        steps.append(
            dict(
                method="update",
                args=[
                    dict(
                        visible=self._visibility_mask(
                            visible_traces=["mean", "std"], fig_data=fig.data
                        )
                    ),
                    dict(title=f"Measured data", annotations=[empty_annotation]),
                ],
                label="Measured data",
            )
        )
        for model, annotation in zip(successful_models, annotations):
            step = dict(
                method="update",
                args=[
                    dict(
                        visible=self._visibility_mask(
                            visible_traces=["mean", "std", model], fig_data=fig.data
                        )
                    ),
                    dict(title=f"Data + Model", annotations=[annotation]),
                ],
                label=f"{model}",
            )

            steps.append(step)

        sliders = [
            dict(
                active=0,
                currentvalue=dict(prefix="", font=dict(color="black")),
                tickcolor="white",
                tickwidth=0,
                font=dict(color="white"),
                pad={"t": 50},
                steps=steps,
            )
        ]

        max_data = np.nanmax([ys["y"] for ys in fig.__dict__["_data_objs"]])

        # # Buttons
        # buttons = []
        # buttons.append(
        #     dict(
        #         method="restyle",
        #         label="Averages",
        #         visible=True,
        #         args=[
        #             dict(
        #                 visible=self._visibility_mask(
        #                     visible_traces=["mean", "std"], fig_data=fig.data
        #                 )
        #             ),
        #         ],
        #         args2=[
        #             dict(
        #                 visible=self._visibility_mask(
        #                     visible_traces=["replicates"], fig_data=fig.data
        #                 ),
        #                 args=[{"annotations": annotation}],
        #             ),
        #         ],
        #     )
        # )

        fig.update_layout(
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
            yaxis_range=[0 - 0.05 * max_data, max_data + 0.05 * max_data],
        )

        # Add title, legend...
        species = self.get_species(visualized_species)
        substrate = self.get_species("substrate")
        fig.update_layout(
            showlegend=True,
            title="Measured data",
            yaxis_title=f"{species.name} ({species.conc_unit.replace(' / l', ' L<sup>-1</sup>')})",
            xaxis_title=f"time ({self._time_unit})",
            hovermode="closest",
            legend_title_text=f"Initial {substrate.name} <br><sub>({self._substrate_unit.replace(' / l', ' L<sup>-1</sup>')})</sub></br>",
            hoverlabel_namelength=-1,
        )

        return fig

    @classmethod
    def from_EnzymeML(
        cls,
        enzmldoc: Union[EnzymeMLDocument, Path],
        substrate_id: str = "s0",
        measured_species_id: str = "s0",
        protein_id: str = "p0",
        inhibitor_id: str = None,
    ):
        def get_unit(unit: str) -> str:
            if "mole" in unit:
                return unit.replace("mole", "mol")
            else:
                return unit

        if isinstance(enzmldoc, str) or isinstance(enzmldoc, Path):
            enzmldoc = EnzymeMLDocument.fromFile(enzmldoc)

        elif isinstance(enzmldoc, EnzymeMLDocument):
            enzmldoc = enzmldoc

        else:
            raise ValueError(
                f"enzmldoc is of type. Needs to be eighter EnzymeMLDocument or string-like path to the omex file."
            )

        pH = enzmldoc.getReaction("r0").ph
        temperature = enzmldoc.getReaction("r0").temperature
        temperature_unit = enzmldoc.getReaction("r0").temperature_unit

        measurements = []
        for measurement in enzmldoc.measurement_dict.values():
            substrate = measurement.getReactant(substrate_id)
            enzyme = measurement.getProtein(protein_id)
            measured_species = measurement.getReactant(measured_species_id)

            substrate_species = Species(
                id=substrate_id,
                name=enzmldoc.getReactant(substrate_id).name,
                initial_conc=substrate.init_conc,
                conc_unit=get_unit(substrate.unit),
                species_type=SpeciesTypes.SUBSTRATE.value,
            )

            replicates = [
                Series(
                    values=rep.data,
                    value_unit=get_unit(substrate.unit),
                    time=rep.time,
                    time_unit=rep.time_unit,
                )
                for rep in measured_species.replicates
            ]
            if len(replicates) == 0:
                raise ValueError(
                    f"Species {measured_species_id} does not contain measurement data. Specify the according 'measured_species_id' from the EnzymeMLDocument"
                )

            if substrate_id == measured_species_id:
                substrate_species.data = replicates
                product_species = []
            else:
                product_species = Species(
                    id=measured_species_id,
                    name=enzmldoc.getReactant(measured_species_id).name,
                    initial_conc=measured_species.init_conc,
                    conc_unit=get_unit(measured_species.unit),
                    species_type=SpeciesTypes.PRODUCT.value,
                    data=replicates,
                )

            enzyme_species = Species(
                id=protein_id,
                name=enzmldoc.getProtein(protein_id).name,
                initial_conc=enzyme.init_conc,
                species_type=SpeciesTypes.ENZYME.value,
            )

            species = [substrate_species, enzyme_species, product_species]

            if inhibitor_id != None:
                inhibitor = measurement.getReactant(inhibitor_id)

                inhibitor_species = Species(
                    id=inhibitor_id,
                    name=enzmldoc.getReactant(inhibitor_id).name,
                    initial_conc=inhibitor.init_conc,
                    conc_unit=get_unit(inhibitor.unit),
                    species_type=SpeciesTypes.INHIBITOR.value,
                )

                species = [
                    substrate_species,
                    enzyme_species,
                    product_species,
                    inhibitor_species,
                ]

            measurements.append(
                Measurement(
                    temperature=temperature,
                    temperature_unit=temperature_unit,
                    pH=pH,
                    species=species,
                )
            )

        enzyme_kinetics = EnzymeKinetics(name=enzmldoc.name, measurements=measurements)

        return cls(data=enzyme_kinetics)
