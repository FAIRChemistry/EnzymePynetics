from EnzymeKinetics.core.enzymekinetics import EnzymeKinetics
from EnzymeKinetics.core.stoichiometrytypes import StoichiometryTypes

import numpy as np


class Kinetics():
    def __init__(self, data: 'EnzymeKinetics'):
        self.data = data
        self._initialize_measurement_data()

    def _initialize_measurement_data(self):
        print(self.data.stoichiometry)
        print(StoichiometryTypes.SUBSTRATE.value)

        measurement_data = []
        initial_substrate = []
        enzyme = []

        for measurement in self.data.measurements:
            for replica in measurement.data:
                measurement_data.append(replica.values)
                initial_substrate.append(measurement.initial_substrate_conc)
                enzyme.append(measurement.enzyme_conc)

        self.initial_substrate = initial_substrate
        self.enzyme = enzyme

        if self.data.stoichiometry == StoichiometryTypes.SUBSTRATE.value:
            self.substrate = measurement_data
            self.product = self._calculate_product()
        elif self.data.stoichiometry == "product": #TODO Fix missing Stoichiometry type
            self.product = measurement_data
            self.substrate = self._calculate_substrate()
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


        

if __name__ == "__main__":
    from EnzymeKinetics.core.measurement import Measurement
    import matplotlib.pyplot as plt
    import numpy as np



    m1 = Measurement(
        initial_substrate_conc=1,
        enzyme_conc=0.05,
    )
    m1.add_to_data([0.0,1,2,3,4,5])
    m1.add_to_data([0.3,1.3,2.3,3.3,4.3,5.3])
    m1.add_to_data([0.6,1.6,2.6,3.6,4.6,5.6])

    m2 = Measurement(
    initial_substrate_conc=2,
    enzyme_conc=0.05,
    )
    m2.add_to_data([1,2,3,4,5,6])
    m2.add_to_data([1.3,2.3,3.3,4.3,5.3,6.3])
    m2.add_to_data([1.6,2.6,3.6,4.6,5.6,6.6])

    m3 = Measurement(
    initial_substrate_conc=2,
    enzyme_conc=0.11,
    )
    m3.add_to_data([10,11,21,31,41,51])
    m3.add_to_data([10.3,11.3,21.3,31.3,41.3,51.3])
    m3.add_to_data([10.6,11.6,21.6,31.6,41.6,51.6])

    testdata = EnzymeKinetics(
        data_conc_unit="mmole / l",
        stoichiometry="substrate",
        data_conc="mmole / l",
        time_unit="min",
        title="test data",
        reactant_name="test product",
        measurements=[m1,m2,m3],
        time=[0,1,2,3,4,5]
    )

    test = Kinetics(testdata)
    print(test.product)
    print(test.substrate)


    for s in test.substrate:
        plt.scatter(np.arange(len(s)), s)
    for p in test.product:
        plt.scatter(np.arange(len(p)), p)
    plt.show()

