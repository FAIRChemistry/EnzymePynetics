

from EnzymePynetics.core.enzymekinetics import EnzymeKinetics
from EnzymePynetics.tools.fitter import Fitter

class PE():
    def __init__(self, kinetics: EnzymeKinetics):
        self.kinetics = kinetics

    def fit_models(self, fitter: Fitter):
        pass
