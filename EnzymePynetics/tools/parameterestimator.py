from ast import Raise, Sub
import time
from turtle import mode, position
from typing import List, Dict, Literal, Union
from matplotlib.pyplot import flag

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
from itertools import groupby
from pathlib import Path
from IPython.display import display
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import matplotlib.colors


class ParameterEstimator:
    def __init__(self, data: EnzymeKinetics, measured_species: SpeciesTypes = None):

        # Identify measured species (Substrate or Product)
        if measured_species:
            self._measured_species = measured_species
        else:
            self._measured_species = self._get_measured_species(
                data.measurements)

        # Get names of all species
        self.substrate_name = self._get_species_name(
            data, SpeciesTypes.SUBSTRATE)
        self.product_name = self._get_species_name(
            data, SpeciesTypes.PRODUCT)
        self.enzyme_name = self._get_species_name(
            data, SpeciesTypes.ENZYME)
        self.inhibitor_name = self._get_species_name(
            data, SpeciesTypes.INHIBITOR)
        self.measured_species_name = self._get_species_name(
            data, self._measured_species)

        # Get and verify units of all species
        (
            self.substrate_unit,
            self.product_unit,
            self.enzyme_unit,
            self.inhibitor_unit,
            self.time_unit
        ) = self._get_units(
            data)

        # Get measurement data
        self.data = self._create_dataframe(data)

        # Create model attribute
        self.models: Dict[str, KineticModel] = None

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
        # if np.any([initial_substrate_concs, start_time_index, stop_time_index]):
        #     (
        #         self.subset_substrate,
        #         self.subset_product,
        #         self.subset_enzyme,
        #         self.subset_initial_substrate,
        #         self.subset_time,
        #         self.subset_inhibitor,
        #     ) = self._subset_data(
        #         initial_substrates=initial_substrate_concs,
        #         start_time_index=start_time_index,
        #         stop_time_index=stop_time_index,
        #     )
        # else:
        #     self.subset_substrate = self.substrate
        #     self.subset_product = self.product
        #     self.subset_enzyme = self.enzyme
        #     self.subset_initial_substrate = self.initial_substrate
        #     self.subset_time = self.time
        #     self.subset_inhibitor = self.inhibitor

        fitting_time, fitting_data = self._prepare_fitting_data()
        substrate, enzyme, product, inhibitor = fitting_data

        y0s = fitting_data[:, :, 0]

        # Initialize kinetics models
        self.models = self._initialize_models(
            y0s=y0s,
            only_irrev_MM=only_irrev_MM,
            inhibitor_species=np.all(inhibitor == 0),
            init_kcat=self._calculate_kcat(substrate, enzyme, fitting_time),
            init_Km=np.nanmax(self._calculate_rates(
                substrate, fitting_time) / 2)
        )

        self._run_minimization(display_output, substrate, fitting_time)

        # Set units in ModelResults object
        for model_name, model in self.models.items():
            for p, parameter in enumerate(model.result.parameters):
                if parameter.name == "k_cat" or parameter.name == "K_ie":
                    self.models[model_name].result.parameters[
                        p
                    ].unit = f"1 / {self.time_unit}"
                elif parameter.name == "Km":
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self.substrate_unit
                elif np.all(inhibitor != 0):
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self.inhibitor_unit
                else:
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self.substrate_unit

        self.result_dict = self._result_overview(inhibitor)
        if display_output:
            display(self.result_dict)

    def get_model_results(self, model: str = None):
        if model == None:
            model = self.result_dict.index[0]

        return self.models[model].result.params

    def _get_init_conc(self, measurement: Measurement, species_type: SpeciesTypes) -> float:
        """Extracts the initial concentration of a species for a given measurement"""

        try:
            return self._get_species(measurement, species_type).initial_conc

        except StopIteration:
            if species_type is SpeciesTypes.INHIBITOR:
                return 0
            else:
                raise StopIteration()

    def _create_dataframe(self, kinetic_data: EnzymeKinetics) -> pd.DataFrame:
        """Extracts measurement data and conditions from ```EnzymeKinetics``` object.
        Returns DataFrame."""

        entries = []

        # Extract data of measured species and measurement conditions
        for measurement in kinetic_data.measurements:

            measured_species = self._get_species(
                measurement, self._measured_species)
            self._measured_species_name = measured_species.name

            inhibitor = self._get_init_conc(
                measurement, SpeciesTypes.INHIBITOR)

            init_substrate = self._get_init_conc(
                measurement, SpeciesTypes.SUBSTRATE)
            enzyme = self._get_init_conc(measurement, SpeciesTypes.ENZYME)

            for replicate_id, replicate in enumerate(measured_species.data):
                for time, concentration in zip(replicate.time, replicate.values):
                    entries.append(
                        {
                            SpeciesTypes.INHIBITOR.value: inhibitor,
                            "init_substrate": init_substrate,
                            "replicate": replicate_id,
                            "time": time,
                            SpeciesTypes.ENZYME.value: enzyme,
                            self._measured_species.value: concentration,
                        }
                    )
        dframe = pd.DataFrame.from_dict(entries)

        dframe = self._calculate_missing_species(
            dframe, self._measured_species)

        dframe = dframe.set_index(
            [SpeciesTypes.INHIBITOR.value, "init_substrate",
                "replicate"]
        ).sort_index()

        return dframe

        # return dframe

    def _get_species_name(self, kinetic_data: EnzymeKinetics, species: SpeciesTypes):
        """Returns the name of a species, if species is not defined,
        none is returned"""

        try:
            return self._get_species(kinetic_data.measurements[0], species).name

        except StopIteration:
            return None

    def visualize_data(self, species: SpeciesTypes = None):
        """Plots data of the measured species

        Args:
            species (SpeciesTypes, optional): Species to visualize. Defaults to None.

        Returns:
            px.Figure: Plot
        """
        if isinstance(species, str):
            species = SpeciesTypes(species.lower())

        if not species:
            species = self._measured_species

        dframe = self.data.sort_index().reset_index()

        fig = px.scatter(data_frame=dframe,
                         y=species.value,
                         x="time",
                         color=dframe["init_substrate"].astype(str),
                         animation_frame=SpeciesTypes.INHIBITOR.value)

        fig.update_layout(
            legend=dict(title=f"Initial {self.substrate_name} ({self.substrate_unit})"))
        fig.update_traces(customdata=["raw"])

        # Set axis labels with corresponding unit
        fig.update_xaxes(dict(title=f"time ({self.time_unit})"))

        if species == SpeciesTypes.SUBSTRATE:
            ylable = self.substrate_name
            unit = self.substrate_unit
        if species == SpeciesTypes.PRODUCT:
            ylable = self.product_name
            unit = self.product_unit
        if species == SpeciesTypes.ENZYME:
            ylable = self.enzyme_name
            unit = self.enzyme_unit
        if species == SpeciesTypes.INHIBITOR:
            ylable = self.inhibitor_name
            unit = self.inhibitor_unit

        ylable = self.substrate_name if species == SpeciesTypes.SUBSTRATE else self.product_name

        fig.update_yaxes(dict(title=f"{ylable} ({unit})"))

        # Set slider label
        fig.update_layout(
            sliders=[dict(currentvalue=dict(prefix=f"{self.inhibitor_name} ",
                                            suffix=f" {self.inhibitor_unit}"))])

        return fig

    def _get_units(self, data: EnzymeKinetics):
        """Extracts and verifies units for measurement data."""

        # Substrate unit
        try:
            substrate_units = [self._get_species(
                meas, SpeciesTypes.SUBSTRATE).conc_unit for meas in data.measurements]

            if self.all_equal(substrate_units):
                substrate_unit = substrate_units[0]
            else:
                raise ValueError(
                    f"{SpeciesTypes.SUBSTRATE.value} species unit is not identical for all measurments."
                )

        except StopIteration:
            raise ValueError(
                f"{SpeciesTypes.SUBSTRATE.value} species not found in measurements. "
                f"Add {SpeciesTypes.SUBSTRATE.value} species to each measurement with respective "
                "unit and initial concentration."
            )

        # Product unit
        try:
            product_units = [self._get_species(
                meas, SpeciesTypes.PRODUCT).conc_unit for meas in data.measurements]

            if self.all_equal(product_units):
                product_unit = product_units[0]
            else:
                raise ValueError(
                    f"{SpeciesTypes.SUBSTRATE.value} species unit is not identical for all measurments."
                )

        except StopIteration:
            product_unit = substrate_unit

        if substrate_unit != product_unit:
            raise ValueError(
                f"{SpeciesTypes.SUBSTRATE.value} unit ({substrate_unit}) and {SpeciesTypes.PRODUCT.value} "
                f"unit ({product_unit}) need to be identical."
            )

        # Enzyme unit
        try:
            enzyme_units = [self._get_species(
                meas, SpeciesTypes.ENZYME).conc_unit for meas in data.measurements]
            if self.all_equal(enzyme_units):
                enzyme_unit = enzyme_units[0]
            else:
                raise ValueError(
                    f"{SpeciesTypes.SUBSTRATE.value} species unit is not identical for all measurments."
                )

        except StopIteration:
            raise ValueError(
                f"{SpeciesTypes.ENZYME.value} species not found in measurements. "
                f"Add {SpeciesTypes.ENZYME.value} species to each measurement with respective "
                "unit and initial concentration."
            )

        # Inhibitor unit
        try:
            inhibitor_units = [self._get_species(
                meas, SpeciesTypes.INHIBITOR).conc_unit for meas in data.measurements]
            if self.all_equal(inhibitor_units):
                inhibitor_unit = inhibitor_units[0]
            else:
                raise ValueError(
                    f"{SpeciesTypes.SUBSTRATE.value} species unit is not identical for all measurments."
                )

        except StopIteration:
            inhibitor_unit = substrate_unit

        # Time unit
        time_units = []
        for meas in data.measurements:
            for species in meas.species:
                time_units.append(species.time_unit)

        time_units = list(filter(lambda item: item is not None, time_units))

        if self.all_equal(time_units):
            time_unit = time_units[0]
        else:
            raise ValueError(
                f"Time unit is not identical for all measurments.",
            )

        return (substrate_unit, product_unit, enzyme_unit, inhibitor_unit, time_unit)

    def _prepare_fitting_data(self) -> tuple:
        """Creates arrays for all substrate, enzyme, product, and inhibitor concentrations.
        Returns tuple(time_array, measurement_arrays)"""

        inhibitor_lvls, init_substrate_lvls, replicate_lvls = [
            lvl.values for lvl in self.data.index.levels]

        dframe = self.data.reset_index()

        fitting_data = []
        time_data = []
        for inhibitor_lvl in inhibitor_lvls:
            for init_substrate_lvl in init_substrate_lvls:
                for replicate_lvl in replicate_lvls:
                    entries = dframe.loc[(dframe[SpeciesTypes.INHIBITOR.value] == inhibitor_lvl) & (
                        dframe['init_substrate'] == init_substrate_lvl) & (dframe['replicate'] == replicate_lvl)].T

                    if not entries.empty:
                        time = entries.loc["time"].values
                        substrate = entries.loc[SpeciesTypes.SUBSTRATE.value].values
                        enzyme = entries.loc[SpeciesTypes.ENZYME.value].values
                        inhibitor = entries.loc[SpeciesTypes.INHIBITOR.value].values
                        product = entries.loc[SpeciesTypes.PRODUCT.value].values

                        fitting_data.append(
                            [substrate, enzyme, product, inhibitor])
                        time_data.append(time)

        time_data = np.array(time_data)
        fitting_data = np.array(fitting_data).swapaxes(
            0, 1)  # Substrate, Enzyme, Product, Inhibitor

        return (time_data, fitting_data)

    def _calculate_rates(self, substrate, time):
        """Calculates the change per time-unit between all measurement points."""
        concentration_intervals = np.diff(substrate)
        time_intervals = np.diff(time)
        rates = abs(concentration_intervals / time_intervals)
        return rates

    def _calculate_kcat(self, substrate, enzyme, time) -> float:
        rates = self._calculate_rates(substrate, time)
        kcat = np.nanmax(rates / enzyme[:, :rates.shape[-1]])
        return kcat

    def _calculate_Km(self, substrate, time):
        # print(self._calculate_rates(substrate, time))
        rates = self._calculate_rates(substrate, time)
        return np.nanmax(self._calculate_rates(substrate, time) / 2)

    # def _subset_data(
    #     self,
    #     initial_substrates: list = None,
    #     start_time_index: int = None,
    #     stop_time_index: int = None,
    # ) -> tuple:
    #     """This function allows to subset the actual measurement data. Thereby, measurements of specific initial substrate concentrations
    #     can be specified. Additionally, the time-frame can be specified by defining the index of the first and last measurement time-point.

    #     Args:
    #         initial_substrates (list, optional). Defaults to None.
    #         start_time_index (int, optional). Defaults to None.
    #         stop_time_index (int, optional). Defaults to None.

    #     Raises:
    #         ValueError: If concentrations are passed, which are not defined in "initial_substrate_concentration"

    #     Returns:
    #         tuple: Data subset.
    #     """
    #     idx = np.array([])
    #     if initial_substrates == None or len(initial_substrates) == 0:
    #         idx = np.arange(self.substrate.shape[0])
    #     else:
    #         for concentration in initial_substrates:
    #             if concentration not in self.initial_substrate:
    #                 raise ValueError(
    #                     f"{concentration} not found in initial substrate concentrations. \nInitial substrate concentrations are {list(np.unique(self.initial_substrate))}"
    #                 )
    #             else:
    #                 idx = np.append(
    #                     idx, np.where(self.initial_substrate ==
    #                                   concentration)[0]
    #                 )
    #     idx = idx.astype(int)

    #     new_substrate = self.substrate[idx, start_time_index:stop_time_index]
    #     new_product = self.product[idx, start_time_index:stop_time_index]
    #     new_enzyme = self.enzyme[idx, start_time_index:stop_time_index]
    #     new_initial_substrate = self.initial_substrate[idx]
    #     new_time = self.time[start_time_index:stop_time_index]
    #     new_inhibitor = self.inhibitor[idx]

    #     return (
    #         new_substrate,
    #         new_product,
    #         new_enzyme,
    #         new_initial_substrate,
    #         new_time,
    #         new_inhibitor,
    #     )

    def get_parameter_dict(self):
        return self.models[self.result_dict.index[0]].result.params

    def _initialize_models(
        self, y0s, init_Km: float, init_kcat: float, inhibitor_species: bool, only_irrev_MM: bool = False
    ) -> Dict[str, KineticModel]:
        """Initializer for all kinetic models. If an inhibitor is provided, inhibition models are initialized. If no inhibitor is specified,
        inhibitory models for substrate and product inhibition are initialized additionally to the irreversible Michaelis Menten model.
        """

        substrate, enzyme, product, inhibitor = y0s

        y0 = np.array([substrate, enzyme, product, inhibitor]).T

        irreversible_Michaelis_Menten_inactivation = KineticModel(
            name="irreversible Michaelis Menten with enzyme inactivation",
            params=[],
            y0=y0,
            kcat_initial=init_kcat,
            Km_initial=init_Km,
            model=irreversible_model,
            enzyme_inactivation=True,
        )

        irreversible_Michaelis_Menten = KineticModel(
            name="irreversible Michaelis Menten",
            params=[],
            y0=y0,
            kcat_initial=init_kcat,
            Km_initial=init_Km,
            model=irreversible_model,
            enzyme_inactivation=False,
        )

        if only_irrev_MM:
            return {
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
            }

        if inhibitor_species:
            # Product inhibition models

            y0 = np.array([substrate, enzyme, product, inhibitor]).T

            competitive_product_inhibition = KineticModel(
                name="competitive product inhibition",
                params=["K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=competitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            competitive_product_inhibition_inactivation = KineticModel(
                name="competitive product inhibition with enzyme inactivation",
                params=["K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=competitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            uncompetitive_product_inhibition = KineticModel(
                name="uncompetitive product inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=uncompetitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            uncompetitive_product_inhibition_inactivation = KineticModel(
                name="uncompetitive product inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=uncompetitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            noncompetitive_product_inhibition = KineticModel(
                name="non-competitive product inhibition",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=noncompetitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            noncompetitive_product_inhibition_inactivation = KineticModel(
                name="non-competitive product inhibition with enzyme inactivation",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=noncompetitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            # Substrate inhibition models

            y0 = np.array([substrate, enzyme, product, substrate]).T

            substrate_inhibition = KineticModel(
                name="substrate inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=substrate_inhibition_model,
                enzyme_inactivation=False,
            )
            substrate_inhibition_inactivation = KineticModel(
                name="substrate inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
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
            y0 = np.array([substrate, enzyme, product, inhibitor]).T

            competitive_inhibition = KineticModel(
                name="competitive inhibition",
                params=["K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=competitive_inhibition_model,
                enzyme_inactivation=False,
            )
            competitive_inhibition_inactivation = KineticModel(
                name="competitive inhibition with enzyme inactivation",
                params=["K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=competitive_inhibition_model,
                enzyme_inactivation=True,
            )

            uncompetitive_inhibition = KineticModel(
                name="uncompetitive inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=uncompetitive_inhibition_model,
                enzyme_inactivation=False,
            )
            uncompetitive_inhibition_inactivation = KineticModel(
                name="uncompetitive inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=uncompetitive_inhibition_model,
                enzyme_inactivation=True,
            )

            noncompetitive_inhibition = KineticModel(
                name="non-competitive inhibition",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=noncompetitive_inhibition_model,
                enzyme_inactivation=False,
            )
            noncompetitive_inhibition_inactivation = KineticModel(
                name="non-competitive inhibition with enzyme inactivation",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=noncompetitive_inhibition_model,
                enzyme_inactivation=True,
            )

            partially_competitive_inhibition = KineticModel(
                name="partially competitive inhibition",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=partially_competitive_inhibition_model,
                enzyme_inactivation=False,
            )
            partially_competitive_inhibition_inactivation = KineticModel(
                name="partially competitive inhibition with enzyme inactivation",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=partially_competitive_inhibition_model,
                enzyme_inactivation=True,
            )
            competitive_inhibition_with_substrate_inhibition_model

            comp_inhib_sub_inhib = KineticModel(
                name="competitive inhibition with uncompetitive substrate inhibition",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=competitive_inhibition_with_substrate_inhibition_model,
                enzyme_inactivation=False,
            )
            comp_inhib_sub_inhib_enz_inact = KineticModel(
                name="competitive inhibition with uncompetitive substrate inhibition with enzyme inactivation",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=init_kcat,
                Km_initial=init_Km,
                model=competitive_inhibition_with_substrate_inhibition_model,
                enzyme_inactivation=True,
            )

            return {
                comp_inhib_sub_inhib_enz_inact.name: comp_inhib_sub_inhib_enz_inact,
                comp_inhib_sub_inhib.name: comp_inhib_sub_inhib,
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
                competitive_inhibition.name: competitive_inhibition,
                competitive_inhibition_inactivation.name: competitive_inhibition_inactivation,
                uncompetitive_inhibition.name: uncompetitive_inhibition,
                uncompetitive_inhibition_inactivation.name: uncompetitive_inhibition_inactivation,
                # noncompetitive_inhibition.name: noncompetitive_inhibition,
                # noncompetitive_inhibition_inactivation.name: noncompetitive_inhibition_inactivation,
                partially_competitive_inhibition.name: partially_competitive_inhibition,
                partially_competitive_inhibition_inactivation.name: partially_competitive_inhibition_inactivation,
            }

    def _run_minimization(self, display_output: bool, substrate, time):
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

            kineticmodel.fit(substrate, time)
            if display_output:
                print(
                    f" -- Fitting succeeded: {kineticmodel.result.fit_success}")

    def _result_overview(self, inhibitor) -> pd.DataFrame:
        """
        Prettifies the results of all kinetic models and organized them in a pandas DataFrame.
        """

        if np.all(inhibitor == 0):
            inhibitor_unit = self.substrate_unit
        else:
            inhibitor_unit = self.inhibitor_unit

        parameter_mapper = {
            "k_cat": f"kcat [1/{self.time_unit}]",
            "Km": f"Km [{self.substrate_unit}]",
            "K_ic": f"Ki competitive [{inhibitor_unit}]",
            "K_iu": f"Ki uncompetitive [{inhibitor_unit}]",
            "K_ie": f"ki time-dep enzyme-inactiv. [1/{self.time_unit}]",
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
                    f"kcat / Km [1/{self.time_unit} * 1/{self.substrate_unit}]"
                ] = f"{kcat_Km:.3f} +/- {percentual_kcat_Km_stderr:.2f}%"

                result_dict[model.name] = {
                    "AIC": aic, "RMSD": rmsd, **parameter_dict}

        df = pd.DataFrame.from_dict(result_dict).T.sort_values(
            "AIC", ascending=True)
        df.fillna("-", inplace=True)
        return df.style.background_gradient(cmap="Blues")

    def _get_visualization_data(self):
        """Gets y0 tuples and respective time array for unique initial conditions in the data set"""

        inhibitor_lvls, init_substrate_lvls, replicate_lvls = [
            lvl.values for lvl in self.data.index.levels]

        dframe = self.data.reset_index()

        results = []

        for inhibitor_lvl in inhibitor_lvls:
            for init_substrate_lvl in init_substrate_lvls:
                substrates = []
                products = []
                for replicate_lvl in replicate_lvls:
                    entries = dframe.loc[(
                        dframe[SpeciesTypes.INHIBITOR.value] == inhibitor_lvl) & (
                        dframe['init_substrate'] == init_substrate_lvl) & (
                        dframe["replicate"] == replicate_lvl
                    )]

                    if not entries.empty:
                        substrates.append(
                            entries[SpeciesTypes.SUBSTRATE.value].iloc[0])
                        products.append(
                            entries[SpeciesTypes.PRODUCT.value].iloc[0])
                        enzyme = entries[SpeciesTypes.ENZYME.value].iloc[0]
                        inhibitor = entries[SpeciesTypes.INHIBITOR.value].iloc[0]
                        measurement_time = entries["time"].values

                mean_substrate = np.mean(substrates)
                mean_product = np.mean(products)

                y0 = (mean_substrate, enzyme,
                      mean_product, inhibitor)

                for model in self.models.values():
                    if model.result.fit_success:

                        dense_time = np.linspace(
                            min(measurement_time), max(measurement_time), 100)

                        datas = model.integrate(
                            model._fit_result.params,
                            [dense_time],
                            [y0]
                        )

                        for data, time in zip(datas[0], dense_time):
                            sub, enz, prod, inhib = data

                            results.append(
                                {
                                    "model": model.name,
                                    f"{SpeciesTypes.SUBSTRATE.value}_simulated": sub,
                                    f"{SpeciesTypes.ENZYME.value}_simulated": enz,
                                    f"{SpeciesTypes.PRODUCT.value}_simulated": prod,
                                    f"{SpeciesTypes.INHIBITOR.value}_simulated": inhib,
                                    "time_simulated": time,
                                    "init_substrate": init_substrate_lvl,
                                    "inhibitor": inhibitor_lvl,
                                }
                            )

        new_df = pd.DataFrame(results).sort_index()
        combined_df = pd.concat(
            [self.data.reset_index(), new_df], ignore_index=True).set_index([f"{SpeciesTypes.INHIBITOR.value}", "model", "init_substrate"]).sort_index()

        return combined_df

    def _style_parameters(self, model: KineticModel):
        param_name_map = dict(
            k_cat="<b><i>k</i><sub>cat</sub>:</b>",
            Km="<b><i>K</i><sub>M</sub>:</b>",
            K_ie="<b><i>K</i><sub>ie</sub>:</b>",
            K_ic="<b><i>K</i><sub>ic</sub>:</b>",
            K_iu="<b><i>K</i><sub>iu</sub>:</b>",
        )
        params = ""
        for parameter in model.result.parameters:
            params = (
                params
                + f"{param_name_map[parameter.name]} "
                + f"{parameter.value:.3f} "
                + f"{self._format_unit(parameter.unit)} \n\n"
            )
        return params

    def visualize_fit(self, visualized_species: SpeciesTypes = None):

        # Select which species to plot
        if visualized_species is None:
            visualized_species = self._measured_species
            visualized_species_name = self._measured_species_name

        if visualized_species == SpeciesTypes.SUBSTRATE.value:
            visualized_species = SpeciesTypes.SUBSTRATE
            visualized_species_name = self.substrate_name

        if visualized_species == SpeciesTypes.PRODUCT.value:
            visualized_species = SpeciesTypes.PRODUCT
            visualized_species_name = self.product_name

        # Initialize figure
        colors = px.colors.qualitative.Plotly
        dframe = self._get_visualization_data()

        # Add labels for subplots, displaying different inhibitor concentrations
        inhibitor_levels = dframe.index.get_level_values(0).unique().values
        model_levels = dframe.index.get_level_values(1).unique().values
        init_substrate_levels = dframe.index.get_level_values(
            2).unique().values

        subplot_titles = []
        if inhibitor_levels.size > 1:

            for row, inhibitor_conc in enumerate(inhibitor_levels):
                if inhibitor_conc == 0:
                    subplot_titles.append(f"without {self.inhibitor_name}")
                else:
                    subplot_titles.append(
                        f"{inhibitor_conc} {self._format_unit(self.inhibitor_unit)} {self.inhibitor_name}"
                    )

        fig = make_subplots(
            rows=inhibitor_levels.size, cols=1,
            y_title=f"{visualized_species_name} ({self._format_unit(self.substrate_unit)})",
            x_title=f"time ({self._format_unit(self.time_unit)})",
            subplot_titles=subplot_titles,
            vertical_spacing=0.05
        )

        inhibitor_annotations = []
        for count, annotation in enumerate(fig.layout["annotations"]):
            if annotation["y"] == 0:
                x_annotation = annotation
            elif annotation["x"] == 0:
                y_annotation = annotation
            else:
                fig.layout["annotations"][count]["x"] = 0
                fig.layout["annotations"][count]["font"]["size"] = 12
                fig.layout["annotations"][count]["xanchor"] = "left"
                inhibitor_annotations.append(fig.layout["annotations"][count])

        ####################
        ### Add traces to figure ###
        ####################

        annotations = []
        steps = []

        # Add measurement data
        for row, inhibitor_lvl in enumerate(inhibitor_levels):
            for init_substrate_lvl, color in zip(init_substrate_levels, colors):
                df_measurement_data = dframe.loc[inhibitor_lvl,
                                                 float("nan"), init_substrate_lvl]
                fig.add_trace(go.Scatter(
                    x=df_measurement_data["time"].values,
                    y=df_measurement_data[visualized_species.value].values,
                    name=f"{init_substrate_lvl}",
                    mode="markers",
                    marker=dict(color=self.hex_to_rgba(color)),
                    customdata=["raw"],
                    hoverinfo="skip",
                    visible=True
                ), row=row+1, col=1)

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

        # Add fitted models
        successful_models = [model for model in self.models.values(
        ) if model.result.fit_success == True]
        successful_models.sort(
            key=lambda model: model.result.AIC)

        for row, inhibitor_lvl in enumerate(inhibitor_levels):
            for model in successful_models:
                for init_substrate_lvl, color in zip(init_substrate_levels, colors):

                    df_subset = dframe.loc[inhibitor_lvl,
                                           model.name, init_substrate_lvl]
                    fig.add_trace(go.Scatter(
                        x=df_subset["time_simulated"].values,
                        y=df_subset[f"{visualized_species.value}_simulated"].values,
                        name=f"{model.name}",
                        mode="lines",
                        marker=dict(color=self.hex_to_rgba(color)),
                        customdata=[f"{model.name}"],
                        hoverinfo="skip",
                        showlegend=False,
                        visible=False
                    ), row=row+1, col=1)

        steps.append(
            dict(
                method="update",
                args=[
                    dict(
                        visible=self._visibility_mask(
                            visible_traces=["raw"], fig_data=fig.data
                        )
                    ),
                    dict(title=f"", annotations=[
                         annotations[0]] + [y_annotation] + [x_annotation] + inhibitor_annotations),
                ],
                label=f"-",
            )
        )

        label_pos = -0.13 if inhibitor_levels.size > 1 else -0.4

        for model in successful_models:
            annotations.append(
                go.layout.Annotation(
                    font=dict(color="black", size=10),
                    x=0,
                    y=label_pos,
                    showarrow=False,
                    text=f"<b>AIC:</b> {round(model.result.AIC)}\n\n {self._style_parameters(model)}",
                    textangle=0,
                    xref="paper",
                    yref="paper",
                    xanchor="left",
                )
            )

        for model, annotation in zip(successful_models, annotations[1:]):

            step = dict(
                method="update",
                args=[
                    dict(
                        visible=self._visibility_mask(
                            visible_traces=["raw", model.name], fig_data=fig.data
                        )
                    ),
                    dict(title=f"",
                         annotations=[annotation] + [y_annotation] + [x_annotation] + inhibitor_annotations),
                ],
                label=f"{model.name}",
            )

            steps.append(step)

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
        print(len(annotations))
        print(len(successful_models))

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
            # yaxis_range=[-0.05 * np.nanmax(ydata), 1.05 * np.nanmax(ydata)],
            # xaxis_range=[-0.05 * np.nanmax(xdata), 1.05 * np.nanmax(xdata)],
        )

        fig.update_layout(height=400 + len(inhibitor_levels) * 150)

        names = set()
        fig.for_each_trace(
            lambda trace:
                trace.update(showlegend=False)
                if (trace.name in names) else names.add(trace.name))

        fig.update_layout(
            showlegend=True,
            hovermode="closest",
            legend_title_text=f"Initial {self.substrate_name} ({self._format_unit(self.substrate_unit)})",
            hoverlabel_namelength=-1,
        )

        config = {
            "toImageButtonOptions": {
                "format": "svg",  # one of png, svg, jpeg, webp
                "filename": "custom_image",
                # "height": 600,
                # "width": 700,
                "scale": 1,  # Multiply title/legend/axis/canvas sizes by this factor
            }
        }

        return fig.show(config=config)

    @classmethod
    def from_EnzymeML(
        cls,
        enzmldoc: Union[EnzymeMLDocument, Path],
        substrate_id: str = "s0",
        product_id: str = "s1",
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
                f"enzmldoc is of type {type(enzmldoc)}. Needs to be eighter EnzymeMLDocument or string-like path to an omex archive."
            )

        pH = enzmldoc.getReaction("r0").ph
        temperature = enzmldoc.getReaction("r0").temperature
        temperature_unit = enzmldoc.getReaction("r0").temperature_unit

        if measured_species_id == substrate_id:
            measured_species_type = SpeciesTypes.SUBSTRATE

        else:
            measured_species_type = SpeciesTypes.PRODUCT

        measurements = []
        for measurement in enzmldoc.measurement_dict.values():
            substrate = measurement.getReactant(substrate_id)
            product = measurement.getReactant(product_id)
            enzyme = measurement.getProtein(protein_id)
            measured_species = measurement.getReactant(measured_species_id)

            substrate_species = Species(
                id=substrate_id,
                name=enzmldoc.getReactant(substrate_id).name,
                initial_conc=substrate.init_conc,
                conc_unit=get_unit(substrate.unit),
                species_type=SpeciesTypes.SUBSTRATE,
                time_unit=[
                    rep.time_unit for rep in measured_species.replicates][0]
            )

            replicates = [
                Series(values=rep.data, time=rep.time) for rep in measured_species.replicates
            ]
            if len(replicates) == 0:
                raise ValueError(
                    f"Species {measured_species_id} does not contain measurement data. Specify 'measured_species_id' according to according in the EnzymeMLDocument"
                )

            if substrate_id == measured_species_id:
                substrate_species.data = replicates
                product_species = Species(
                    id=product_id,
                    name=enzmldoc.getReactant(product_id).name,
                    initial_conc=product.init_conc,
                    conc_unit=get_unit(product.unit),
                    species_type=SpeciesTypes.PRODUCT,
                    data=[],
                )
            else:
                product_species = Species(
                    id=measured_species_id,
                    name=enzmldoc.getReactant(measured_species_id).name,
                    initial_conc=measured_species.init_conc,
                    conc_unit=get_unit(measured_species.unit),
                    species_type=SpeciesTypes.PRODUCT,
                    data=replicates,
                )

            enzyme_species = Species(
                id=protein_id,
                name=enzmldoc.getProtein(protein_id).name,
                initial_conc=enzyme.init_conc,
                species_type=SpeciesTypes.ENZYME,
                conc_unit=get_unit(enzyme.unit)
            )

            species = [substrate_species, enzyme_species, product_species]

            if inhibitor_id != None:
                inhibitor = measurement.getReactant(inhibitor_id)

                inhibitor_species = Species(
                    id=inhibitor_id,
                    name=enzmldoc.getReactant(inhibitor_id).name,
                    initial_conc=inhibitor.init_conc,
                    conc_unit=get_unit(inhibitor.unit),
                    species_type=SpeciesTypes.INHIBITOR,
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

        enzyme_kinetics = EnzymeKinetics(
            name=enzmldoc.name, measurements=measurements)

        return cls(data=enzyme_kinetics, measured_species=measured_species_type)

    @staticmethod
    def all_equal(iterable):
        g = groupby(iterable)
        return next(g, True) and not next(g, False)

    @staticmethod
    def _get_species(measurement: Measurement, species_type: SpeciesTypes) -> Species:
        """Returns the respective species of a measurement"""

        return next(species for species in measurement.species
                    if species.species_type == species_type.value)

    @staticmethod
    def _repeat_value(time_array, species_array):
        """Repeats initial conditions of a species for n measurements, according to time axis of 
        respective measurement."""

        for idx, time in enumerate(time_array):
            species_array[idx] = np.repeat(species_array[idx], len(time))

        return np.array(species_array)

    @staticmethod
    def _get_y0s(substrate, enzyme, product, inhibitor):
        y0s = []
        for s, e, p, i in zip(substrate, enzyme, product, inhibitor):
            y0s.append((s[0], e[0], p[0], i[0]))

        return np.array(y0s)

    @staticmethod
    def _HEX_to_RGBA_string(color: list) -> str:
        return f"rgba{tuple(color)}"

    @staticmethod
    def _visibility_mask(visible_traces: list, fig_data: list) -> list:
        return [
            any(fig["customdata"][0] == trace for trace in visible_traces)
            for fig in fig_data
        ]

    @staticmethod
    def _format_unit(unit: str) -> str:
        unit = unit.replace(" / l", " L<sup>-1</sup>")
        unit = unit.replace("1 / s", "s<sup>-1</sup>")
        unit = unit.replace("1 / min", "min<sup>-1</sup>")
        unit = unit.replace("umol", "µmol")
        unit = unit.replace("ug", "µg")
        return unit

    @staticmethod
    def group_lists(measurement_data: List[List[float]], init_inhibitors: List[float]):
        result = []
        temp = []
        prev_val = init_inhibitors[0]

        for i in range(len(measurement_data)):
            if init_inhibitors[i] != prev_val:
                result.append(temp)
                temp = []

            temp.append(measurement_data[i])
            prev_val = init_inhibitors[i]

        result.append(temp)
        return result

    @staticmethod
    def _get_measured_species(measurements: List[Measurement]) -> SpeciesTypes:
        """Checks if substrate or product data is provided."""

        measured_species = []
        for measurement in measurements:
            for species in measurement.species:
                if species.data:
                    measured_species.append(species.species_type)

        if len(measured_species) > 1:
            raise ValueError(
                "Data contains measurments for substrate and product. Define which species should be used for "
                "parameter estimation, by specifying 'measured_species'. Using substrate and product data "
                "is currently not supported."
            )

        if not measured_species:
            raise ValueError("Data contains no measurment data.")

        if len(measured_species) == 1:
            return measured_species[0]

    @staticmethod
    def _calculate_missing_species(dframe: pd.DataFrame, measured_species: SpeciesTypes) -> pd.DataFrame:
        """Calculates missing species (e.g. substrate if product was measured)
        and adds it to the ```DataFrame```."""

        if measured_species is SpeciesTypes.SUBSTRATE:
            dframe[SpeciesTypes.PRODUCT.value] = dframe["init_substrate"] - \
                dframe[SpeciesTypes.SUBSTRATE.value]

            return dframe

        if measured_species is SpeciesTypes.PRODUCT:
            dframe[SpeciesTypes.SUBSTRATE.value] = dframe["init_substrate"] - \
                dframe[SpeciesTypes.PRODUCT.value]

            return dframe

        raise ValueError(
            "Measured species not defined."
        )

    # staticmethod
    def get_first_n_elements(tup, time):
        new_tuple = ()
        for arr in tup:
            new_arr = arr[:n]
            new_tuple += (new_arr,)
        return new_tuple

    @staticmethod
    def hex_to_rgba(hex: str) -> str:
        rgb = tuple(int(hex.strip("#")[i:i+2], 16) for i in (0, 2, 4))
        return f"rgba{rgb + (255,)}"
