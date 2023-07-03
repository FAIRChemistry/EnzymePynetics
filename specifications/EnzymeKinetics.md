# Enzyme kinetics data model

The following data model describes (meta-)data from an enzyme kinetics experiment in a structured way. It constists out of multiple ```measurements```, describing one or multiple measurements at different initial substrate and/or enzyme concentrations. If and inhibitor was applied to the reaction, its concentration can be speciefied as well.

### EnzymeKineticsExperiment

Base class, dealing with measurement data of an enzyme kinetics assay.

- __title__
  - Type: string
  - Description: Title of the kinetic experiment
- __reactant_name__
  - Type: string
  - Description: Name of the measured reactant.
- __temperature__
  - Type: float
  - Description: Temperature of the reaction.
- __temperature_unit__
  - Type: string
  - Description: Temperature unit.
- __pH__
  - Type: float
  - Description: pH of the reaction.
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
- __enzyme_conc__
  - Type: float
  - Description: Enzyme concentration in the reaction.
- __data*__
  - Type: Series
  - Description: One or multiple time-course concentration data arrays.
  - Multiple: True
- __inhibitor_conc__
  - Type: float
  - Description: Inhibitor concentration, if applied to the reaction.
- __inhibitor_conc_unit__
  - Type: ConcentrationTypes
  - Description: Inhibitor concentration in the reaction, if applied.

### Series

Time-course data of an individual reaction.

- __values__
  - Type: float
  - Description: Time-course data of an individual reaction.
  - Multiple: True
- __test__
  - Type: string
  - Description: Test field

#### StoichiometryTypes

Measurement data can eighter be substrate or product.

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
GRAMLITER = "g / l"
MILLIGRAMLITER = "mg / l"
MICROGRAMLITER = "ug / l"
NANGRAMLITER = "ng / l"
```

#### TimeTypes

Allowed time types.

```python
S = "s"
MIN = "min"
H = "h"
```
