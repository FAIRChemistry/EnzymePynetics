from ast import Raise, Sub
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
        self.data = data
        self.models: Dict[str, KineticModel] = None

        if measured_species:
            self._measured_species = measured_species
        else:
            self._measured_species = self._get_measured_species(
                data.measurements)

        self.dframe = self._initialize_measurement_data()
        # self.initial_kcat = self._calculate_kcat()
        # self.initial_Km = self._calculate_Km()

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
        print(f"conc_unit: {self._conc_unit}")
        for model_name, model in self.models.items():
            for p, parameter in enumerate(model.result.parameters):
                if parameter.name == "k_cat" or parameter.name == "K_ie":
                    self.models[model_name].result.parameters[
                        p
                    ].unit = f"1 / {self._time_unit}"
                elif parameter.name == "Km":
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self._conc_unit
                elif np.any(self.subset_inhibitor != 0):
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self._inhibitor_conc_unit
                else:
                    self.models[model_name].result.parameters[
                        p
                    ].unit = self._conc_unit

        self.result_dict = self._result_overview()
        if display_output:
            display(self.result_dict)

    def get_model_results(self, model: str = None):
        if model == None:
            model = self.result_dict.index[0]

        return self.models[model].result.params

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

    def _initialize_dataframe(self) -> pd.DataFrame:
        """Extracts measurement data and conditions from ```EnzymeKinetics``` object.
        Returns DataFrame."""

        entries = []

        # Extract data of measured species and measurement conditions
        for measurement in self.data.measurements:

            measured_species = self._get_species(
                measurement, self._measured_species)
            inhibitor = self._get_init_conc(
                measurement, SpeciesTypes.INHIBITOR)
            init_substrate = self._get_init_conc(
                measurement, SpeciesTypes.SUBSTRATE)
            enzyme = self._get_init_conc(measurement, SpeciesTypes.ENZYME)

            for replicate in measured_species.data:
                for time, concentration in zip(replicate.time, replicate.values):
                    entries.append(
                        {
                            SpeciesTypes.INHIBITOR.value: inhibitor,
                            "init_substrate": init_substrate,
                            SpeciesTypes.ENZYME.value: enzyme,
                            "time": time,
                            self._measured_species.value: concentration,
                        }
                    )
        dframe = pd.DataFrame.from_dict(entries)

        dframe = self._calculate_missing_species(dframe)

        dframe = dframe.set_index(
            [SpeciesTypes.INHIBITOR.value, "init_substrate",
                SpeciesTypes.ENZYME.value, "time"]
        ).sort_index()

        return dframe

    def _calculate_missing_species(self, dframe: pd.DataFrame) -> pd.DataFrame:
        """Calculates missing species (e.g. substrate if product was measured)
        and adds it to the ```DataFrame```."""

        if self._measured_species is SpeciesTypes.SUBSTRATE:
            dframe[SpeciesTypes.PRODUCT.value] = dframe["init_substrate"] - \
                dframe[SpeciesTypes.SUBSTRATE.value]

            return dframe

        elif self._measured_species is SpeciesTypes.PRODUCT:
            dframe[SpeciesTypes.SUBSTRATE.value] = dframe["init_substrate"] - \
                dframe[SpeciesTypes.SUBSTRATE.value]

            return dframe

        else:
            raise ValueError(
                "Messured species not defined."
            )

    def _get_units(self, data: EnzymeKinetics, species: SpeciesTypes):

        # Substrate unit
        try:
            substrate_units = [self._get_species(
                meas, SpeciesTypes.SUBSTRATE) for meas in data.measurements]
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
                meas, SpeciesTypes.PRODUCT) for meas in data.measurements]

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
                meas, SpeciesTypes.ENZYME) for meas in data.measurements]
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
                meas, SpeciesTypes.INHIBITOR) for meas in data.measurements]
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
                if species.data:
                    for replicate in species.data:
                        time_units.append()

        # self._time_unit = replicate.time_unit
        # self._conc_unit = replicate.values_unit

        # try:
        #     self._inhibitor_conc_unit = self._get_species(
        #         measurement, SpeciesTypes.INHIBITOR).conc_unit
        # except StopIteration:
        #     self._inhibitor_conc_unit = None

    def all_equal(iterable):
        g = groupby(iterable)
        return next(g, True) and not next(g, False)

    @staticmethod
    def _get_species(measurement: Measurement, species_type: SpeciesTypes) -> Species:
        """Returns the respective species of a measurement"""

        return next(species for species in measurement.species
                    if species.species_type == species_type.value)

    def _get_init_conc(self, measurement: Measurement, species_type: SpeciesTypes) -> float:
        """Extracts the initial concentration of a species for a given measurement"""

        try:
            return self._get_species(measurement, species_type).initial_conc

        except StopIteration:
            if species_type is SpeciesTypes.INHIBITOR:
                return 0
            else:
                raise StopIteration()

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
                        idx, np.where(self.initial_substrate ==
                                      concentration)[0]
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
            # Product inhibition models

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

            # Substrate inhibition models

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
            competitive_inhibition_with_substrate_inhibition_model

            comp_inhib_sub_inhib = KineticModel(
                name="competitive inhibition with uncompetitive substrate inhibition",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_with_substrate_inhibition_model,
                enzyme_inactivation=False,
            )
            comp_inhib_sub_inhib_enz_inact = KineticModel(
                name="competitive inhibition with uncompetitive substrate inhibition with enzyme inactivation",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
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
                print(
                    f" -- Fitting succeeded: {kineticmodel.result.fit_success}")

    def _result_overview(self) -> pd.DataFrame:
        """
        Prettifies the results of all kinetic models and organized them in a pandas DataFrame.
        """

        if np.all(self.inhibitor == 0):
            inhibitor_unit = self._conc_unit
        else:
            inhibitor_unit = self._inhibitor_conc_unit

        parameter_mapper = {
            "k_cat": f"kcat [1/{self._time_unit}]",
            "Km": f"Km [{self._conc_unit}]",
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
                    f"kcat / Km [1/{self._time_unit} * 1/{self._conc_unit}]"
                ] = f"{kcat_Km:.3f} +/- {percentual_kcat_Km_stderr:.2f}%"

                result_dict[model.name] = {
                    "AIC": aic, "RMSD": rmsd, **parameter_dict}

        df = DataFrame.from_dict(result_dict).T.sort_values(
            "AIC", ascending=True)
        df.fillna("-", inplace=True)
        return df.style.background_gradient(cmap="Blues")

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

    def visualize_subplots(
        self, visualized_species: SpeciesTypes = None, save_path: str = None
    ):
        colors = matplotlib.colors.to_rgba_array(px.colors.qualitative.Plotly)

        # Select which species to plot
        if visualized_species is None:
            visualized_species = self._measured_species

        if visualized_species == SpeciesTypes.SUBSTRATE.value:
            measurement_data = self.subset_substrate
            species_tuple = 0
            visualized_species = self._measured_species
        else:
            measurement_data = self.subset_product
            species_tuple = 2

        # Get measured data for visualization
        measurement_data, measurement_time = self._get_timecourse_data(
            visualized_species)

        init_substrates = []
        init_enzymes = []
        init_inhibitors = []
        for measurement in self.data.measurements:
            init_substrates.append(self._get_init_conc(
                measurement, SpeciesTypes.SUBSTRATE))
            init_enzymes.append(self._get_init_conc(
                measurement, SpeciesTypes.ENZYME))
            init_inhibitors.append(self._get_init_conc(
                measurement, SpeciesTypes.INHIBITOR))

        # plot replicates
        unique_initial_substrates = np.unique(init_substrates)
        unique_inhibitors = np.unique(init_inhibitors)

        # Reshape data and time array to contain inhibitor and replicate information
        datas = self.group_lists(measurement_data, init_inhibitors)
        times = self.group_lists(measurement_time, init_inhibitors)

        # get measured species and substrate species
        species = self._get_species(measurement, visualized_species)
        substrate = self._get_species(measurement, SpeciesTypes.SUBSTRATE)
        inhibitor = self._get_species(measurement, SpeciesTypes.INHIBITOR)

        subplot_titles = []
        for row, inhibitor_conc in enumerate(unique_inhibitors):
            if inhibitor_conc == 0:
                subplot_titles.append(f"without {inhibitor.name}")
            else:
                subplot_titles.append(
                    f"{inhibitor_conc} {inhibitor.conc_unit} {inhibitor.name}"
                )

        fig = make_subplots(
            rows=len(datas),
            cols=1,
            shared_yaxes="all",
            y_title=f"Initial {substrate.name} ({self._conc_unit})",
            x_title=f"time ({self._time_unit})",
            subplot_titles=subplot_titles,
            horizontal_spacing=0.05,
            vertical_spacing=0.05,
        )

        # Adjust style of subplot label
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

        ydata = [ys["y"] for ys in fig.__dict__["_data_objs"]]
        xdata = [xs["x"] for xs in fig.__dict__["_data_objs"]]

        # Integrate each successfully fitted model

        successfull_models = []
        steps = []
        for model in self.models.values():
            if model.result.fit_success:
                successfull_models.append(model)
                print(model.y0)

                means_y0s = model.y0

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
                                marker=dict(
                                    color=self._HEX_to_RGBA_string(color)),
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
            # yaxis_range=[-0.05 * np.nanmax(ydata), 1.05 * np.nanmax(ydata)],
        )

        # Add title, legend...
        fig.update_layout(
            showlegend=True,
            title="Measured data",
            yaxis_title=f"{species.name} ({species.conc_unit})",
            xaxis_title=f"time ({self._time_unit})",
            hovermode="closest",
            legend_title_text=f"Initial {substrate.name} ({self._conc_unit})",
            hoverlabel_namelength=-1,
        )

        fig.update_layout(height=400 + len(unique_inhibitors) * 150)

        return fig

    def _get_timecourse_data(self, species_type: SpeciesTypes) -> List[List[List[float]]]:
        """Gets measurement data of a species for the data set.
        Returns List[List[List]] representing measurements[replicates[timecourse]] along with 
        time information with the same shape"""

        data = []
        time = []
        for measurement in self.data.measurements:
            species = self._get_species(measurement, species_type)

            replicates = []
            subset_time = []
            for replicate in species.data:
                replicates.append(replicate.values)
                subset_time.append(replicate.time)

            data.append(replicates)
            time.append(subset_time)

        return data, time

    def _get_timecourse_inhibitors(self, species_type: SpeciesTypes) -> List[List[List[float]]]:
        """Gets measurement data of a species for the data set.
        Returns List[List[List]] representing measurements[replicates[timecourse]] along with 
        time information with the same shape"""

        inhibitors = []
        old_inhib_conc = 0
        for measurement in self.data.measurements:
            data = []
            inhib_conc = self._get_init_conc(
                measurement, SpeciesTypes.INHIBITOR)
            species = self._get_species(measurement, species_type)
            if inhib_conc == old_inhib_conc:
                replicates = []
                for replicate in species.data:
                    replicates.append(replicate.values)
                data.append(replicates)
            inhibitors.append(data)
            old_inhib_conc = inhib_conc

        return inhibitors

    def visualize_model_overview(
        self, visualized_species: SpeciesTypes = None, save_path: str = None
    ):
        fig = go.Figure()
        colors = matplotlib.colors.to_rgba_array(px.colors.qualitative.T10_r)

        # Select which species to plot
        if visualized_species is None:
            visualized_species = self._measured_species

        # Get measured data for visualization
        measurement_data, measurement_time = self._get_timecourse_data(
            visualized_species)

        init_substrates = []
        init_enzymes = []
        init_inhibitors = []
        for measurement in self.data.measurements:
            init_substrates.append(self._get_init_conc(
                measurement, SpeciesTypes.SUBSTRATE))
            init_enzymes.append(self._get_init_conc(
                measurement, SpeciesTypes.ENZYME))
            init_inhibitors.append(self._get_init_conc(
                measurement, SpeciesTypes.INHIBITOR))

        # plot replicates

        for i, (data, time, init_sub) in enumerate(
            zip(measurement_data, measurement_time, init_substrates)
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
        if any(len(meas) > 1 for meas in measurement_data):

            for color, replicate_data, init_sub, time in zip(
                colors, measurement_data, init_substrates, measurement_time
            ):
                mean = np.mean(replicate_data, axis=0)
                std = np.std(replicate_data, axis=0)
                fig.add_trace(
                    go.Scatter(
                        x=time[0],
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
                        x=time[0],
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
                        x=time[0],
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

        ydata = [ys["y"] for ys in fig.__dict__["_data_objs"]]
        xdata = [xs["x"] for xs in fig.__dict__["_data_objs"]]

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
                integration_data = model.integrate(
                    model._fit_result.params,
                    self.subset_time,
                    mean_y0s,
                )

                for data, time, color, init_sub in zip(
                    integration_data[:, :, 0],
                    measurement_time,
                    colors,
                    init_substrates,
                ):
                    color[-1] = 1
                    fig.add_trace(
                        go.Scatter(
                            x=time[0],
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
                param_name_map = dict(
                    k_cat="<b><i>k</i><sub>cat</sub>:</b>",
                    Km="<b><i>K</i><sub>M</sub>:</b>",
                    K_ie="<b><i>K</i><sub>ie</sub>:</b>",
                    K_ic="<b><i>K</i><sub>ic</sub>:</b>",
                    K_iu="<b><i>K</i><sub>iu</sub>:</b>",
                )
                params = ""
                for parameter in model.result.parameters:
                    print(parameter)
                    params = (
                        params
                        + f"{param_name_map[parameter.name]} "
                        + f"{parameter.value:.3f} "
                        + f"{self._format_unit(parameter.unit)} \n\n"
                    )
                annotations.append(
                    go.layout.Annotation(
                        font=dict(color="black", size=10),
                        x=-0.05 * np.nanmax(xdata),
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
                    dict(title=f"Measured data",
                         annotations=[empty_annotation]),
                ],
                label="",
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
            yaxis_range=[-0.05 * np.nanmax(ydata), 1.05 * np.nanmax(ydata)],
            xaxis_range=[-0.05 * np.nanmax(xdata), 1.05 * np.nanmax(xdata)],
        )

        # Add title, legend...
        # species = self._get_species(visualized_species)
        # substrate = self._get_species("substrate")
        # fig.update_layout(
        #     showlegend=True,
        #     title="Measured data",
        #     yaxis_title=f"{species.name} ({self._format_unit(species.conc_unit)})",
        #     xaxis_title=f"time ({self._format_unit(self._time_unit)})",
        #     hovermode="closest",
        #     legend_title_text=f"Initial {substrate.name} <br>({self._format_unit(self._conc_unit)})</br>",
        #     hoverlabel_namelength=-1,
        # )

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
            enzyme = measurement.getProtein(protein_id)
            measured_species = measurement.getReactant(measured_species_id)

            substrate_species = Species(
                id=substrate_id,
                name=enzmldoc.getReactant(substrate_id).name,
                initial_conc=substrate.init_conc,
                conc_unit=get_unit(substrate.unit),
                species_type=SpeciesTypes.SUBSTRATE,
            )

            replicates = [
                Series(
                    values=rep.data,
                    values_unit=get_unit(substrate.unit),
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
                    species_type=SpeciesTypes.PRODUCT,
                    data=replicates,
                )

            enzyme_species = Species(
                id=protein_id,
                name=enzmldoc.getProtein(protein_id).name,
                initial_conc=enzyme.init_conc,
                species_type=SpeciesTypes.ENZYME,
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
