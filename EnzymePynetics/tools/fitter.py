
from typing import List
from abc import ABC, abstractmethod

from lmfit import Parameters, minimize, report_fit
import numpy as np

from EnzymePynetics.tools.parameterestimator import ParameterEstimator
from EnzymePynetics.tools.kineticmodel import KineticModel

class Fitter(ABC):

    @abstractmethod
    def residuals():
        "Defines how model equations are initialized."

class RateMM(Fitter):

    def __init__(self, kinetic_model: List[KineticModel]):
        self.kinetic_model = kinetic_model

    def residuals(self, params: Parameters, time, substrate, s0, enzyme):
        pass


class IntegratedMM(Fitter):

    def __init__(self, kinetic_model: List[KineticModel]):
        self.kinetic_model = kinetic_model


class Fitter():
    def __init__(self, solver, measurements) -> None:
        self.solver = solver
        self.measurements = measurements

    def residuals()

    def fit_models(self):
        pass




def irrev_MM(k_cat: float,
             K_m: float,
             enzyme: float,
             substrate: float,
             init_substrate: float,
             t_0: float,
             k_inactivation: float,
             enzyme_inactivation: bool = False,
             ) -> float:
    
    if enzyme_inactivation:
        enzyme = np.exp(-k_inactivation * enzyme)
    
    return -1/(k_cat*enzyme)*(K_m* np.log(substrate/init_substrate) + (substrate-init_substrate)) + t_0



if __name__ == "__main__":
    import matplotlib.pyplot as plt










