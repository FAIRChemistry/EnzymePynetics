
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

    def __init__(self, kinetic_model: List[KineticModel], **kwargs):
        self.kinetic_model = kinetic_model

    def residuals(self):
        pass


class IntegratedMM(Fitter):

    def __init__(self, kinetic_model: List[KineticModel]):
        self.kinetic_model = kinetic_model


class Fitter():
    def __init__(self, solver: Solver, measurements) -> None:
        self.solver = solver
        self.measurements = measurements

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
        e = np.exp(-k_inactivation * e)
    
    return -1/(k_cat*enzyme)*(K_m* np.log(substrate/init_substrate) + (substrate-init_substrate)) + t_0



if __name__ == "__main__":
    import matplotlib.pyplot as plt


    concentration = np.linspace(10,1)
    s0 = 10
    km = 5
    kcat = 34
    t0 = 55
    e = 0.04

    def irrev_mm(st: float):
        return -1/(kcat*e)*(km* np.log(st/s0) + (st-s0)) + t0


    time = [irrev_mm(x) for x in concentration]
    time = np.random.normal(0,0.3, size=concentration.size) +time










