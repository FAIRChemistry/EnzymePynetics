
from typing import List, Dict
from abc import ABC, abstractmethod

from lmfit import Parameters, minimize, report_fit
from scipy.integrate import odeint
import numpy as np


from EnzymePynetics.tools.kineticmodel import KineticModel
from EnzymePynetics.tools.kineticmodel import integrated_MM_model, irreversible_model

class Fitter(ABC):

    @abstractmethod
    def _initialize_models():
        "Defines how model equations are initialized."

    @abstractmethod
    def _residuals():
        "Defines how model equations are initialized."

class RateMM(Fitter):

    def __init__(self,
                 Km_initial: float,
                 k_cat_initial: float,
                 substrate,
                 product,
                 enzyme,
                 inhibitor,
                 substrate_initial,
                 time,
                 enzyme_inactivation: bool):

        self.Km_initial = Km_initial
        self.k_cat_initial = k_cat_initial
        self.substrate = substrate
        self.product = product
        self.enzyme = enzyme
        self.inhibitor = inhibitor
        self.substrate_initial = substrate_initial
        self.time = time
        self.enzyme_inactivation = enzyme_inactivation

        self.w0 = {"cS": substrate,
                   "cE": enzyme,
                   "cP": product,
                   "cI": inhibitor,
                   "cS0": substrate_initial}

        self.models = self._initialize_models()
        self._fit_models()

        print(self.w0)

    def _initialize_models(self) -> Dict[str, KineticModel]:
        irreversible_Michaelis_Menten = KineticModel(
            name="irreversible Michaelis Menten",
            params=[],
            w0=self.w0,
            kcat_initial=self.k_cat_initial,
            Km_initial=self.Km_initial,
            model=irreversible_model,
            enzyme_inactivation=self.enzyme_inactivation
        )

        return {irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten}
    
    def _residuals(self, params: Parameters, time, w0, kinetic_model: KineticModel):


        def integrate(time: np.ndarray, w0: tuple, params: Parameters):
            '''
            Solution to the ODE w'(t)=f(t,w,p) with initial condition w(0)= w0 (= [S0])
            '''
            w = odeint(kinetic_model.model, w0, time, args=(params, False))
            return w
        
        substrate, enzyme, product, inhibitor, init_substrate = w0.values()
        residuals = 0.0 * substrate

        for i, (cS, cE, cP, cI, cS0) in enumerate(zip(substrate, enzyme, product, inhibitor, init_substrate)):
            w0 = (cS0, cE, cP, cI, cS0)

            model = integrate(time, w0, params)

            # get modeled substrate
            model = model[:, 0]

            # compute distance to measured data
            residuals[i] = cS-model
            print("lol")

        return residuals.flatten()
        

    def _fit_models(self):
        for kinetic_model in self.models.values():
            result = minimize(self._residuals, kinetic_model.parameters, args=(
            self.time, self.w0, kinetic_model), method='leastsq', nan_policy='omit')
        
        return result


class IntegratedMM(Fitter):

    def __init__(self,
                Km_initial: float,
                k_cat_initial: float,
                substrate,
                product,
                enzyme,
                inhibitor,
                substrate_initial,
                time,
                enzyme_inactivation: bool):

        self.Km_initial = Km_initial
        self.k_cat_initial = k_cat_initial
        self.substrate = substrate
        self.product = product
        self.enzyme = enzyme
        self.inhibitor = inhibitor
        self.substrate_initial = substrate_initial
        self.time = time
        self.enzyme_inactivation = enzyme_inactivation

        self.w0 = {"cS": substrate,
                   "cE": enzyme,
                   "cP": product,
                   "cI": inhibitor,
                   "cS0": substrate_initial}

        self.models = self._initialize_models()
        self._fit_models()

    def _initialize_models(self) -> Dict[str, KineticModel]:
        irreversible_Michaelis_Menten = KineticModel(
            name="irreversible Michaelis Menten (with time-offset)",
            params=["t_0", "k_inact"],
            w0=self.w0,
            kcat_initial=self.k_cat_initial,
            Km_initial=self.Km_initial,
            model=integrated_MM_inactivation,
            enzyme_inactivation=self.enzyme_inactivation
        )

        return {irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten}

    def _residuals(self, params: Parameters, time, kinetic_model: KineticModel):

        substrate, enzyme, product, inhibitor, init_substrate = self.w0.values()
        residuals = substrate * 0.0

        for i, (cS, cE, cP, cI, cS0) in enumerate(zip(substrate, enzyme, product, inhibitor, init_substrate)):
            w0 = (cS, cE, cP, cI, cS0)
        
            model_time = integrated_MM_inactivation(w0=w0, params=params)

            residuals[i] = model_time - time

        return residuals.flatten()
    
    def _fit_models(self):
        for kinetic_model in self.models.values():
            result = minimize(self._residuals, kinetic_model.parameters, args=(
            self.time, kinetic_model), method='leastsq', nan_policy='omit')
        
        return result

## Kinetic models ##

def integrated_MM(w0, params: Parameters) -> float:

    # parameters
    params = params.valuesdict()
    K_m = params["Km"]
    k_cat = params["k_cat"]
    t_0 = params["t_0"]

    # concentrations
    cS, cE, cP, cI, cS0 = w0

    return -1/(k_cat*cE)*(K_m* np.log(cS/cS0) + (cS-cS0)) + t_0

def integrated_MM_inactivation(w0, params: Parameters) -> float:

    # parameters
    params = params.valuesdict()
    K_m = params["Km"]
    k_cat = params["k_cat"]
    t_0 = params["t_0"]
    k_inact = params["k_inact"]

    # concentrations
    cS, cE, cP, cI, cS0 = w0

    #n1 = (K_m * np.log(cS/cS0))/(cE*k_cat)
    #n2 = cS0 / cE*k_cat
    #n3 = cS / cE*k_cat

    #result = -k_inact * np.log(-(-n1 + n2 -n3)/t_0)

    result = k_inact * np.log(-(-(K_m * np.log(cS/cS0))/(cE * k_cat) + cS0/(cE * k_cat) - cS/(cE * k_cat))/t_0)

    return result

if __name__ == "__main__":
    from sdRDM import DataModel
    import json
    import matplotlib.pyplot as plt


    path = "/Users/max/code/EnzymePynetics/irrev_MM_enzyme_inactivation_with_offset.json"
    with open(path, "r") as f:
        data = json.loads(f.read())

    enzyme = []
    initial_substrate = []
    substrate = []    
    enzyme = []    

    for measurement in data["measurements"]:
        time = measurement["time"]
        enzyme.append(measurement["enzyme_conc"])

        for species in measurement["species"]:
            initial_substrate.append(species["initial_conc"])

            for data in species["data"]:
                substrate.append(data["values"])

    print(initial_substrate)
    print(enzyme)
    print(substrate)

    substrate=np.array(substrate)
    product=np.array(enzyme)
    enzyme=np.array(enzyme)
    inhibitor=np.array(enzyme)
    substrate_initial=np.array(initial_substrate)
    time=np.array(time)



    fit = IntegratedMM(Km_initial=2,
                 k_cat_initial=10,
                 substrate=substrate,
                 product=enzyme,
                 enzyme=enzyme,
                 inhibitor=enzyme,
                 substrate_initial=initial_substrate,
                 time=time,
                 enzyme_inactivation=False)
    
    result = fit._fit_models()

    print(report_fit(result))


    # Plot
    plot = False
    if plot:
        def g(t, w0, params, model):

                '''
                Solution to the ODE w'(t)=f(t,w,p) with initial condition w(0)= w0 = cS
                '''

                w = odeint(model, w0, t, args=(params, False,))
                return w
        

        for i, (cS, cE, cP, cI, cS0) in enumerate(zip(substrate, enzyme, product, inhibitor, substrate_initial)):
            substrate = cS
            w0 = (cS0, cE, cP, cI, cS0)



            data_fitted = g(t=time, w0=w0, params=result.params, model=fit.models[next(iter(fit.models))].model)

            #plt.plot(time, data_fitted.T[0])
            plt.scatter(time, cS, label=f"{cS0}")
            plt.legend(title="initial substrate [mM]")

        plt.ylabel("substrate [mM]")
        plt.xlabel("time [min]")

        plt.savefig("test.svg", transparent=True)

        plt.show()












