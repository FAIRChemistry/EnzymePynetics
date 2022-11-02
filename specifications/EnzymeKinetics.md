# EnzymeKineticsDatamodel

Data-model of pyEnzymeKinetics.

### EnzymeKinetics

- __title__
  - Type: string
  - Description: Title of the kinetic experiment.
- __reactant_name__
  - Type: string
  - Description: Name of the measured reactant.
- __measurements*__
  - Type: Measurement
  - Description: Measurement data for a given initial substrate concentration.
  - Multiple: True
- __stoichiometry__
  - Type: StoichiometryTypes
  - Description: Define whether "substrate" or "product" concentration was measured.
- __data_conc_unit*__
  - Type: ConcentrationTypes
  - Description: Molar concentration unit of the measured data.
- __time*__
  - Type: float
  - Description: Time array corresponding to time-course data.
  - Multiple: True
- __time_unit*__
  - Type: TimeTypes
  - Description: Time data unit.

### Measurement

A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling.

- __initial_substrate_conc*__
  - Type: float
  - Description: Initial substrate concentration of the measurement.
- __enzyme_conc*__
  - Type: float
  - Multiple: True
  - Description: Enzyme concentration in the measurement.
- __inhibitor_conc__
  - Type: float
  - Description: inhibitor concentration for the measurement, if inhibitor was present.
- __inhibitor_conc_unit__
  - Type: str
  - Description: Concentration unit of the inhibitior.
- __data*__
  - Type: Series
  - Description: One or multiple time-course concentration data arrays
  - Multiple: True

### Series
- __values__
  - Type: float
  - Description: Time-course data of an individual reaction.
  - Multiple: True

#### StoichiometryTypes

```python
SUBSTRATE = "substrate"
PRODUCT = "product"
```

#### TimeTypes

```python
S = "s"
MIN = "min"
H = "h"
```

#### ConcentrationTypes

```python
MOLAR = "mole / l"
MILLIMOLAR = "mmole / l"
MICROMOLAR = "umole / l"
NANAMOLAR = "nmole / l"
```
