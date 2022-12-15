from typing import List, Dict, Optional

from pyenzyme import EnzymeMLDocument
from EnzymePynetics.core.enzymekineticsexperiment import EnzymeKineticsExperiment
from EnzymePynetics.core.stoichiometrytypes import StoichiometryTypes
from EnzymePynetics.core.series import Series
from EnzymePynetics.core.measurement import Measurement
from EnzymePynetics.tools.kineticmodel import *

import numpy as np
from pandas import DataFrame
from scipy.integrate import odeint
from lmfit import minimize, report_fit
from IPython.display import display
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt


class ParameterEstimator():

    def __init__(self, data: EnzymeKineticsExperiment):
        self.data = data
        self.models: Dict[str, KineticModel] = None
        self._initialize_measurement_data()
        self._check_negative_concentrations()
        self.initial_kcat = self._calculate_kcat()
        self.initial_Km = self._calculate_Km()
        # TODO shapcheck function to check for consistent array lengths

    def fit_models(
        self,
        initial_substrate_concs: list = None,
        start_time_index: int = None,
        stop_time_index: int = None,
        enzyme_inactivation: bool = False,
        only_irrev_MM: bool = False,
        display_output: bool = True
        ) -> None:
        """Fits the measurement data to a set of kinetic models.

        Args:
            initial_substrate_concs (list, optional): Enables to subset the measurement data by choosing one ore multiple initial substrate concentrations. Defaults to None.
            start_time_index (int, optional): Subset the data by choosing the index of thee first measurement point which should be condidered. Defaults to None.
            stop_time_index (int, optional): Choose last measurement point which should be considered. Defaults to None.
            enzyme_inactivation (bool, optional): _description_. Defaults to False.

        Prints DataFrame of all kinetic parameters of all fitted models sorted by Akaike information criterion.
        """

        self.enzyme_inactivation = enzyme_inactivation

        # Subset data if one or multiple attributes are passed to the function
        if np.any([initial_substrate_concs, start_time_index, stop_time_index]):
            self.subset_substrate, self.subset_product, self.subset_enzyme, self.subset_initial_substrate, self.subset_time, self.subset_inhibitor = self._subset_data(
                initial_substrates=initial_substrate_concs,
                start_time_index=start_time_index,
                stop_time_index=stop_time_index)            
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
            only_irrev_MM=only_irrev_MM)

        self._run_minimization(display_output)

        self.result_dict = self._result_overview()
        if display_output:
            display(self.result_dict)

    def visualize(
        self,
        model_name: Optional[str] = None,
        path: Optional[str] = None,
        title: Optional[str] = None,
        visualize_species: Optional[str] = None,
        plot_means: bool = True,
        ax: plt.Axes = None,
        **plt_kwargs) -> None:
        """Visualizes the measurement data as well as the fitted model. By default the fest fitting model is choosen for visualization 
        (based on Akaice criterion)

        Args:
            model_name (Optional[str], optional): Specify the name of the model which should be visualized. Defaults to None.
            path (Optional[str], optional): Provide path, were to save the plot.. Defaults to None.
            title (Optional[str], optional): Choose alternative title for the plot. Defaults to None.
            visualize_species (Optional[str], optional): Choose whether 'substrate' or 'product' should be visualized. By default the measured species is visualized. Defaults to None.
            plot_means (bool, optional): Decide whether all measured data should be plotted or mean values with their respective standard deviation. Defaults to True.
        """
        if ax is None:
            ax_provided = False
            ax = plt.gca()
        else:
            ax_provided = True


        # Select which model to visualize
        best_model = self.result_dict.index[0]
        if model_name is None:
            model_name = best_model
        model = self.models[model_name]

        # Visualization modes
        plot_modes = {
            "substrate": [self.subset_substrate,0, self.data.reactant_name], # Substrate
            "product": [self.subset_product, 2, self.data.reactant_name], # TODO Product
        }

        if visualize_species is None:
            experimental_data, reactant, name = plot_modes[self.data.stoichiometry]
        else:
            experimental_data, reactant, name = plot_modes[visualize_species]

        def g(t, w0, params):

            '''
            Solution to the ODE w'(t)=f(t,w,p) with initial condition w(0)= w0 = cS
            '''

            w = odeint(model.model, w0, t, args=(params, self.enzyme_inactivation,))
            return w

        if plot_means:
            unique_substrates = np.unique(self.subset_initial_substrate)
            unique_inhibitors = np.unique(self.subset_inhibitor)

            cS, cE, cP, cI = [self._mean_w0(data, self.subset_initial_substrate) for data in model.w0.values()]

            # Markers
            unique_inhibitors = np.unique(self.inhibitor)
            markers = ["o", "x", "D", "X", "d"]
            marker_mapping = dict(zip(unique_inhibitors, markers[:len(unique_inhibitors)]))
            marker_vector = [marker_mapping[item] for item in self.inhibitor[:,0]]

            unique_concs = np.unique(self.subset_initial_substrate)
            colors = get_cmap("tab10").colors
            color_mapping = dict(zip(unique_concs, colors[:len(unique_concs)]))

            for inhibitor, marker in zip(unique_inhibitors, markers):
                # get substrates
                init_inhibitor = self.subset_inhibitor[:,0]

                inhibitor_mask = np.where(init_inhibitor == inhibitor)[0]

                for substrate, color in zip(unique_substrates, colors):
                    idx = inhibitor_mask[np.where(self.subset_initial_substrate[inhibitor_mask] == substrate)[0]]
                    mean = np.mean(experimental_data[idx,:], axis=0)
                    std = np.std(experimental_data[idx,:], axis=0)

                    ax.errorbar(self.subset_time, mean, std, label=substrate, fmt=marker, color=color, **plt_kwargs)

                    cS, cE, cP, cI = [np.mean(data[idx], axis=0) for data in model.w0.values()]
                    w0 = (cS[0], cE[0], cP[0], cI[0])

                    data_fitted = g(t=self.subset_time, w0=w0, params=model.result.params)
                    ax.plot(self.subset_time, data_fitted[:,reactant], color = color)

        else:
            cS, cE, cP, cI = model.w0.values()
            color_vector = [color_mapping[item] for item in self.subset_initial_substrate]

            for i, data in enumerate(experimental_data):

                w0 = (cS[i,0], cE[i,0], cP[i,0], cI[i,0])

                ax.scatter(x=self.subset_time, y=data, label=self.subset_initial_substrate[i], marker=marker_vector[i], color=color_vector[i], **plt_kwargs)

                data_fitted = g(t=self.subset_time, w0=w0, params=model.result.params)

                # Plot model
                ax.plot(self.subset_time, data_fitted[:,reactant], color = color_vector[i])

        if title is None:
            ax.set_title(self.data.title)
        else:
            ax.set_title(title)



        # Legend
        if ax_provided == False:

            ax.set_ylabel(f"{name} [{self.data.data_conc_unit}]")
            ax.set_xlabel(f"time [{self.data.time_unit}]")

            handles, labels = ax.get_legend_handles_labels()

            new_handles, new_labels = [[],[]]
            for handle, label in zip(handles, labels):
                if len(new_labels) == 0:
                    new_labels.append(label)
                    new_handles.append(handle)
                else:
                    if label == new_labels[-1]:
                        pass
                    else:
                        new_labels.append(label)
                        new_handles.append(handle)
            ax.legend(title = f"initial substrate [{self.data.data_conc_unit}]", handles=new_handles, labels=new_labels, loc='center left', bbox_to_anchor=(1, 0.5))
        
        if path != None:
            ax.savefig(path, format="svg")

        if ax_provided == False:
            report_title = f"Fit report for {model.name} model"
            print(f"{report_title}")
            report_fit(model.result)

    def get_model_results(self, model: str = None):
        if model == None:
            model = self.result_dict.index[0]

        return self.models[model].result.params

    def _initialize_measurement_data(self):
        """
        Extracts data from data objects and reshapes it for fitting.
        """

        measurement_data = []
        initial_substrate = []
        enzyme = []
        inhibitor = []

        for measurement in self.data.measurements:
            for replica in measurement.data:
                measurement_data.append(replica.values)
                initial_substrate.append(measurement.initial_substrate_conc)
                enzyme.append(measurement.enzyme_conc)
                if measurement.inhibitor_conc != None:
                    inhibitor.append(measurement.inhibitor_conc)
                else:
                    inhibitor.append(0)

        measurement_shape = np.array(measurement_data).shape

        self.time = np.array(self.data.time)
        self.initial_substrate = np.array(initial_substrate)
        self.enzyme = np.array(enzyme)
        self.inhibitor = np.repeat(np.array(inhibitor),measurement_shape[1]).reshape(measurement_shape)
        self.enzyme = np.repeat(np.array(enzyme),measurement_shape[1]).reshape(measurement_shape)

        if self.data.stoichiometry == StoichiometryTypes.SUBSTRATE.value:
            self.substrate = np.array(measurement_data)
            self.product = np.array(self._calculate_missing_species(self.substrate))
        elif self.data.stoichiometry == StoichiometryTypes.PRODUCT.value: 
            self.product = np.array(measurement_data)
            self.substrate = np.array(self._calculate_missing_species(self.product))
        else:
            raise AttributeError("Please define whether measured data is substrate or product data.")

    def _calculate_missing_species(self, existing_species: np.ndarray) -> List[list]:
        """
        Calcualtes data or the species that was not measured, based on the assumption of mass conservation.
        E.g.: If the substrate depletion was measured, the corresponding product concentration for each 
        measurement point is calculated by:
        P(t) = S0 - S(t) with
        P(t): product concentration at time-point 't'
        S0: initial substrate concentration of the measurement
        S(t): measured substrate concentration at time-point 't'

        Args:
            existing_species (np.ndarray): Existing measurement data.

        Returns:
            List[list]: Calculated species
        """
        species = []
        for measurement, initial_substrate in zip(existing_species, self.initial_substrate):
            species.append(
                [initial_substrate - value for value in measurement])
        return species

    def _calculate_rates(self):
        """
        Calculates the change per time-unit between all measurement points.
        """
        concentration_intervals = np.diff(self.substrate)
        time_intervals = np.diff(self.data.time)
        rates = abs(concentration_intervals / time_intervals)
        return rates

    def _calculate_kcat(self) -> float:
        rates = self._calculate_rates()
        #initial_enzyme_tile = np.repeat(self.enzyme, rates.shape[1]).reshape(rates.shape)
        kcat = np.nanmax(rates / self.enzyme[:,:-1])
        return kcat

    def _calculate_Km(self):
        return np.nanmax(self._calculate_rates()) / 2

    def _check_negative_concentrations(self):
        """
        Checks if negative concentrations are entered or calculated by the "_calculate_missing_species" function.
        If negatice concentrations are calculated, the data might not be blanced correctly, or the 'initial_substrate_concentration'
        is lower than the measurement data.
        """
        self.deactivate_substrate_inhibition = False
        self.deactivate_product_inhibition = False

        try:
            assert np.any(self.substrate<0) == False
        except AssertionError:
            #print("Provided product concentration is higher than specified initial substrate concentration. Calculated substrate concentration results in negative values. Therefore, substrate inhibition models are excluded.")
            self.deactivate_substrate_inhibition = True
        try:
            assert np.any(self.product<0) == False
        except AssertionError:
            #print("Provided substrate concentration is higher than specified initial substrate concentration. Calculated product concentration results in negative values. Therefore, product inhibition models are excluded.")
            self.deactivate_product_inhibition = True

    def _subset_data(self, initial_substrates: list = None, start_time_index: int = None, stop_time_index: int = None) -> tuple:
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
                        f"{concentration} not found in initial substrate concentrations. \nInitial substrate concentrations are {list(np.unique(self.initial_substrate))}")
                else:
                    idx = np.append(idx, np.where(self.initial_substrate == concentration)[0])
        idx = idx.astype(int)

        new_substrate = self.substrate[idx,start_time_index:stop_time_index]
        new_product = self.product[idx,start_time_index:stop_time_index]
        new_enzyme = self.enzyme[idx,start_time_index:stop_time_index]
        new_initial_substrate = self.initial_substrate[idx]
        new_time = self.time[start_time_index:stop_time_index]
        new_inhibitor = self.inhibitor[idx]

        return (new_substrate, new_product, new_enzyme, new_initial_substrate, new_time, new_inhibitor)

    def get_parameter_dict(self):
        return self.models[self.result_dict.index[0]].result.params

    def _initialize_models(self, substrate, product, enzyme, inhibitor, only_irrev_MM: bool) -> Dict[str, KineticModel]:
        """Initializer for all kinetic models. If an inhibitor is provided, inhibition models are initialized. If no inhibitor is specified,
        inhibitory models for substrate and product inhibition are inizialized additionally to the irreversible Michaelis Menten model.
        """

        irreversible_Michaelis_Menten = KineticModel(
            name="irreversible Michaelis Menten",
            params=[],
            w0={"cS": substrate, "cE": enzyme,
                "cP": product, "cI": inhibitor},            
            kcat_initial=self.initial_kcat,
            Km_initial=self.initial_Km,
            model=irreversible_model,
            enzyme_inactivation=self.enzyme_inactivation
        )

        if only_irrev_MM:
            return {irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten}

        if np.all(self.inhibitor == 0):

            competitive_product_inhibition = KineticModel(
                name="competitive product inhibition",
                params=["K_ic"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": product},            
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_product_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )
            uncompetitive_product_inhibition = KineticModel(
                name="uncompetitive product inhibition",
                params=["K_iu"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": product},            
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_product_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )
            noncompetitive_product_inhibition = KineticModel(
                name="non-competitive product inhibition",
                params=["K_iu", "K_ic"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": product},            
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_product_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )
            substrate_inhibition = KineticModel(
                name="substrate inhibition",
                params=["K_iu"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": substrate},            
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=substrate_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )

            model_dict = {irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten}
            if not self.deactivate_product_inhibition:
                model_dict[competitive_product_inhibition.name] = competitive_product_inhibition
                model_dict[uncompetitive_product_inhibition.name] = uncompetitive_product_inhibition
                model_dict[noncompetitive_product_inhibition.name] = noncompetitive_product_inhibition
            if not self.deactivate_substrate_inhibition:
                model_dict[substrate_inhibition.name] = substrate_inhibition

            return model_dict
               
        else:
            competitive_inhibition = KineticModel(
                name="competitive inhibition",
                params=["K_ic"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": inhibitor},            
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )
            uncompetitive_inhibition = KineticModel(
                name="uncompetitive inhibition",
                params=["K_iu"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": inhibitor},            
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )
            noncompetitive_inhibition = KineticModel(
                name="non-competitive inhibition",
                params=["K_iu", "K_ic"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": inhibitor},
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )
            partially_competitive_inhibition = KineticModel(
                name="partially competitive inhibition",
                params=["K_ic", "K_iu"],
                w0={"cS": substrate, "cE": enzyme,
                    "cP": product, "cI": inhibitor},            
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=partially_competitive_inhibition_model,
                enzyme_inactivation=self.enzyme_inactivation
            )
            return {
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                competitive_inhibition.name: competitive_inhibition,
                uncompetitive_inhibition.name: uncompetitive_inhibition,
                noncompetitive_inhibition.name: noncompetitive_inhibition,
                partially_competitive_inhibition.name: partially_competitive_inhibition,
            }

    def _run_minimization(self, display_output: bool) -> DataFrame:
        """Performs non-linear least-squared minimization to fit the data to the kinetic 
        models by adjusting the kinetic parameters of the models.

        Returns:
            DataFrame: Overfiew of the kinetic parameters of all kinetic models.
        """

        if display_output:
            print("Fitting data to:")
        for kineticmodel in self.models.values():
            if display_output:
                print(f" - {kineticmodel.name} model")

            def g(time: np.ndarray, w0: tuple, params):
                '''
                Solution to the ODE w'(t)=f(t,w,p) with initial condition w(0)= w0 (= [S0])
                '''
                w = odeint(kineticmodel.model, w0, time, args=(params, self.enzyme_inactivation,))
                return w

            def residual(params, time: np.ndarray, substrate: np.ndarray):
                """
                Calculated the distance between measured and modeled data.
                """
                residuals = 0.0 * substrate
                for i, measurement in enumerate(substrate):

                # Calculate residual for each measurement
                    cS, cE, cP, cI = kineticmodel.w0.values()
                    w0 = (cS[i, 0], cE[i, 0], cP[i, 0], cI[i, 0])

                    model = g(time, w0, params)  # solve the ODE with sfb.

                    # get modeled substrate
                    model = model[:, 0]

                    # compute distance to measured data
                    residuals[i] = measurement-model

                return residuals.flatten()

            kineticmodel.result = minimize(residual, kineticmodel.parameters, args=(
                self.subset_time, self.subset_substrate), method='leastsq', nan_policy='omit')
            

    def _result_overview(self) -> DataFrame:
        """
        Prettifies the results of all kinetic models and organized them in an pandas DataFrame.
        """

        if np.all(self.inhibitor == 0):
            inhibitor_unit = self.data.data_conc_unit
        else:
            inhibitor_unit = self.data.measurements[0].inhibitor_conc_unit

        parameter_mapper = {
            "k_cat": f"kcat [1/{self.data.time_unit}]",
            "Km": f"Km [{self.data.data_conc_unit}]",
            "K_ic": f"Ki competitive [{inhibitor_unit}]",
            "K_iu": f"Ki uncompetitive [{inhibitor_unit}]",
            "K_ie": f"ki time-dep enzyme-inactiv. [1/{self.data.time_unit}]",

        }

        result_dict = {}
        for model in self.models.values():
            name = model.name
            aic = round(model.result.aic)

            parameter_dict = {}
            for parameter in model.result.params.values():
                name = parameter_mapper[parameter.name]
                value = parameter.value
                stderr = parameter.stderr

                try:
                    percentual_stderr = stderr / value * 100
                except TypeError:
                    percentual_stderr = float("nan")

                if name.startswith("Ki time-dep"):
                    parameter_dict[name] = f"{value:.5f} +/- {percentual_stderr:.2f}%"
                else:
                    parameter_dict[name] = f"{value:.3f} +/- {percentual_stderr:.2f}%"

                if parameter.name == "k_cat":
                    kcat = parameter.value
                    kcat_stderr = parameter.stderr
                if parameter.name == "Km":
                    Km = parameter.value
                    Km_stderr = parameter.stderr

            if Km_stderr is None or kcat_stderr is None:
                kcat_Km_stderr = float("nan")
                kcat_Km = float("nan")
                percentual_kcat_Km_stderr = float("nan")
            else:
                kcat_Km = kcat / Km
                kcat_Km_stderr =((kcat_stderr / kcat)**2+(Km_stderr / Km)**2)**0.5 * kcat_Km
                percentual_kcat_Km_stderr = kcat_Km_stderr / kcat_Km * 100

            parameter_dict[f"kcat / Km [1/{self.data.time_unit} * 1/{self.data.data_conc_unit}]"] = f"{kcat_Km:.3f} +/- {percentual_kcat_Km_stderr:.2f}%"


            result_dict[model.name] = {"AIC": aic, **parameter_dict}

        df = DataFrame.from_dict(result_dict).T.sort_values("AIC", ascending=True)
        df.fillna('-', inplace=True)

        return df

    def _calculate_mean_std(self, data: np.ndarray):
        """Calculated mean and stddev of data for visualization.
        """
        mean_data = np.array([])
        std_data = np.array([])
        unique_initial_substrates = np.unique(self.subset_initial_substrate)
        if np.any(self.inhibitor > 0):
            unique_inhibitors = np.unique(self.inhibitor)
        for concentration in unique_initial_substrates:
            idx = np.where(self.subset_initial_substrate == concentration)
            mean_data = np.append(mean_data, np.mean(data[idx], axis=0))
            std_data = np.append(std_data, np.std(data[idx], axis=0))
            self.inhibitor

        mean_data = mean_data.reshape(len(unique_initial_substrates), int(len(mean_data)/len(unique_initial_substrates)))
        std_data = std_data.reshape(len(unique_initial_substrates), int(len(std_data)/len(unique_initial_substrates)))
        
        return mean_data, std_data

    def _mean_w0(self, measurement_data: np.ndarray, init_substrate):
        """Calculates the mean values of substrate product and inhibitor to initialize the fitter, if "plot_means" is set in the 'visualize' method.
        """
        unique = np.unique(init_substrate)
        mean_array = np.array([])
        for concentration in unique:
            idx = np.where(init_substrate == concentration)
            mean = np.mean(measurement_data[idx], axis=0)
            if mean_array.size == 0:
                mean_array = np.append(mean_array, mean)
            else:
                mean_array = np.vstack([mean_array, mean])
        return mean_array

    @classmethod
    def from_EnzymeML(
        cls,
        enzmldoc: EnzymeMLDocument,
        reactant_id: str,
        measured_species: StoichiometryTypes,
        substrate_id: str = "s0",
        inhibitor_id: str = None,
        protein_id: str = "p0"
        ):

        pH = enzmldoc.getReaction("r0").ph
        temperature = enzmldoc.getReaction("r0").temperature
        temperature_unit = enzmldoc.getReaction("r0").temperature_unit

        measurements = []
        for measurement in enzmldoc.measurement_dict.values():
            substrate = measurement.getReactant(substrate_id)
            reactant = measurement.getReactant(reactant_id)
            protein = measurement.getReactant(protein_id)
            if inhibitor_id != None:
                inhibitor = measurement.getReactant(inhibitor_id)
                inhibitor_conc = inhibitor.init_conc
                inhibitor_conc_unit = inhibitor.unit
            else:
                inhibitor_conc = 0
                inhibitor_conc_unit = None

            reps = [Series(values=reps.data) for reps in reactant.replicates]

            measurements.append(Measurement(
                initial_substrate_conc=substrate.init_conc,
                enzyme_conc=protein.init_conc,
                inhibitor_conc=inhibitor_conc,
                inhibitor_conc_unit=inhibitor_conc_unit,
                data=reps))

        experimental_data = EnzymeKineticsExperiment(
            data_conc_unit=reactant.replicates[0].data_unit,
            time_unit=reactant.replicates[0].time_unit,
            title=enzmldoc.name,
            temperature=temperature,
            temperature_unit=temperature_unit,
            pH=pH,
            reactant_name=enzmldoc.getReactant(reactant_id).name,
            measurements=measurements,
            stoichiometry=measured_species,
            time = reactant.replicates[0].time
        )

        return cls(experimental_data)