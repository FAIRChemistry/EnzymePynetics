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

Description.

- __test__
  - Type: string
  - Description: Test field.

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
