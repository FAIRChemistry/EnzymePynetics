from EnzymeKinetics.core.enzymekinetics import EnzymeKinetics
from EnzymeKinetics.core.stoichiometrytypes import StoichiometryTypes

import numpy as np


class ParameterEstimator():
    def __init__(self, data: 'EnzymeKinetics'):
        self.data = data
        self._initialize_measurement_data()
        self._check_negative_concentrations()
        self.initial_kcat = self._calculate_kcat()
        self.initial_Km = self._calculate_Km()

    def _initialize_measurement_data(self):

        measurement_data = []
        initial_substrate = []
        enzyme = []

        for measurement in self.data.measurements:
            for replica in measurement.data:
                measurement_data.append(replica.values)
                initial_substrate.append(measurement.initial_substrate_conc)
                enzyme.append(measurement.enzyme_conc)

        self.initial_substrate = np.array(initial_substrate)
        self.enzyme = np.array(enzyme)

        if self.data.stoichiometry == StoichiometryTypes.SUBSTRATE.value:
            self.substrate = np.array(measurement_data)
            self.product = np.array(self._calculate_product())
        elif self.data.stoichiometry == StoichiometryTypes.PRODUCT.value: 
            self.product = np.array(measurement_data)
            self.substrate = np.array(self._calculate_substrate())
        else:
            raise AttributeError("Please define whether measured data is substrate or product data.")

    def _calculate_substrate(self):
        substrate = []
        for product, initial_substrate in zip(self.product, self.initial_substrate):
            substrate.append(
                [initial_substrate - value for value in product])
        return substrate

    def _calculate_product(self):
        product = []
        for substrate, initial_substrate in zip(self.substrate, self.initial_substrate):
            product.append(
                [initial_substrate - value for value in substrate])
        return product

    def _calculate_rates(self):
        concentration_intervals = np.diff(self.substrate)
        time_intervals = np.diff(self.data.time)
        rates = abs(concentration_intervals / time_intervals)
        return rates

    def _calculate_kcat(self) -> float:
        rates = self._calculate_rates()
        initial_enzyme_tile = np.repeat(self.enzyme, rates.shape[1]).reshape(rates.shape)
        kcat = np.nanmax(rates / initial_enzyme_tile)
        return kcat

    def _calculate_Km(self):
        return np.nanmax(self._calculate_rates()) / 2

    def _check_negative_concentrations(self):
        if np.any(self.substrate<0):
            raise ValueError(
                "Substrate data contains negative concentrations. Check data.")        
        if np.any(self.product<0):
            raise ValueError(
                "Product data contains negative concentrations. Check data.")

    def fit_models(
        self,
        initial_substrate_concs: list = None,
        start_time_index: int = None,
        stop_time_index: int = None,
        flag_enzyme_inactivation: bool = False,
        ):

        # Subset data if any attribute is passed to the function
        if initial_substrate_concs != None or np.any([start_time_index, stop_time_index]):
            substrate, product, enzyme, initial_substrate = self._subset_data(
                initial_substrates=initial_substrate_concs,
                start_time_index=start_time_index,
                stop_time_index=stop_time_index)            
        else:
            substrate = self.substrate
            product = self.product
            enzyme = self.enzyme
            initial_substrate = self.initial_substrate

        return substrate
        

    def _subset_data(self, initial_substrates: list = None, start_time_index: int = None, stop_time_index: int = None) -> tuple:
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
        new_enzyme = self.enzyme[idx]
        new_initial_substrate = self.initial_substrate[idx]

        return (new_substrate, new_product, new_enzyme, new_initial_substrate)

if __name__ == "__main__":
    from EnzymeKinetics.core.measurement import Measurement
    import matplotlib.pyplot as plt
    import numpy as np



    m1 = Measurement(
        initial_substrate_conc=100,
        enzyme_conc=0.05,
    )
    m1.add_to_data([0.0,1,2,3,4,5])
    m1.add_to_data([0.3,1.3,2.3,3.3,4.3,5.3])
    m1.add_to_data([0.6,1.6,2.6,3.6,4.6,5.6])

    m2 = Measurement(
    initial_substrate_conc=200,
    enzyme_conc=0.05,
    )
    m2.add_to_data([1,2,3,4,5,6])
    m2.add_to_data([1.3,2.3,3.3,4.3,5.3,6.3])
    m2.add_to_data([1.6,2.6,3.6,4.6,5.6,6.6])

    m3 = Measurement(
    initial_substrate_conc=300,
    enzyme_conc=0.11,
    )
    m3.add_to_data([10,11,21,31,41,51])
    m3.add_to_data([10.3,11.3,21.3,31.3,41.3,51.3])
    m3.add_to_data([10.6,11.6,float("nan"),31.6,41.6,71.6])

    testdata = EnzymeKinetics(
        data_conc_unit="mmole / l",
        stoichiometry="product",
        data_conc="mmole / l",
        time_unit="min",
        title="test data",
        reactant_name="test product",
        measurements=[m1,m2,m3],
        time=[0,2,4,6,8,9]
    )


    test = ParameterEstimator(testdata)
    print(test.fit_models(initial_substrate_concs=[],start_time_index=4))




