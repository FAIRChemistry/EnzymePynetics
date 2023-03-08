from lmfit import Parameters
from lmfit.minimizer import MinimizerResult
from typing import Dict, Callable, Tuple
from numpy import ndarray

class Model():
    def __init__(self,
                 name: str,
                 model: Callable,
                 params: Parameters,
                 kcat_initial: float,
                 Km_initial: float,
                 w0: Dict[str, ndarray],
                 enzyme_inactivation: bool = False
                 ) -> None:
        
        self.name = name
        self.model = model
        self.params = Parameters()
        self.enzyme_inactivation = enzyme_inactivation
        self.w0 = w0
        self.kcat_initial = kcat_initial
        self.Km_initial = Km_initial
        self.parameters = self._set_parameters(params)
        self.result: MinimizerResult = None