# EnzymeKineticsDatamodel

### EnzymeKinetics

Base class, dealing with measurement data of an enzyme kinetics assay.

- __title__
  - Type: string
  - Description: Title of the kinetic experiment
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

- __initial_conc__
  - Type: float
  - Description: Initial substrate concentration of the measurement.

#### StoichiometryTypes

Description.

```python
SUBSTRATE = "substrate"
PRODUCT = "product"
```

#### ConcentrationTypes

```python
MOLAR = "mole / l"
MILLIMOLAR = "mmole / l"
MICROMOLAR = "umole / l"
NANAMOLAR = "nmole / l"
```

#### TimeTypes

Description.

```python
S = "s"
MIN = "min"
H = "h"
```
