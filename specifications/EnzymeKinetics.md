# EnzymePynetics data model

The following data model describes (meta-)data from an enzyme kinetics experiment in a structured way. It constists out of multiple ```measurements```, describing one or multiple measurements at different initial substrate and/or enzyme concentrations. If and inhibitor was applied to the reaction, its concentration can be speciefied as well.

## Objects

### EnzymeKinetics

Base class, dealing with measurement data of an enzyme kinetics assay.

- title
  - Type: string
  - Description: Title of the kinetic experiment.
- model_results
  - Type: ModelResult
  - Description: Fitted kinetic models which were used for parameter estimation.
  - Multiple: True
- measurements
  - Type: Measurement
  - Description: Measurement data for a given initial substrate concentration.
  - Multiple: True

### Species

- name
  - Type: string
  - Description: name of the reactant.
- conc_unit
  - Type: ConcentrationTypes
  - Description: Concentration unit of the measurement data.
- time_unit
  - Type: TimeTypes
  - Description: Time data unit.
- initial_conc
  - Type: float
  - Description: Initial concentration of the reactant.
- species_type
  - Type: SpeciesTypes
  - Description: Define the role of the species in the reaction.
- data
  - Type: Series
  - Description: One or multiple time-course measurement data arrays.
  - Multiple: True

### Measurement

A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling.

- species
  - Type: Species
  - Description: Reactants of the reaction.
  - Multiple: True
- temperature
  - Type: float
  - Description: Temperature of the reaction.
- temperature_unit
  - Type: string
  - Description: Temperature unit.
- pH
  - Type: float
  - Description: pH of the reaction

### ModelResult

Description of a kinetic model

- name
  - Type: string
  - Description: Name of the kinetic model.
- equations
  - Type: string
  - Description: Equation of the kinetic model.
  - Multiple: True
- parameters
  - Type: Parameter
  - Description: Kinetic parameters of the model.
  - Multiple: True
- fit_success
  - Type: bool
  - Description: Whether or not model fitting was possible.
- AIC
  - Type: float
  - Description: Akaike information criterion.
- BIC
  - Type: float
  - Description: Bayesian information criterion.
- RMSD
  - Type: float
  - Description: Root mean square deviation between model and measurement data.

### Parameter

Defines a kinetic parameter.

- name
  - Type: string
  - Description: Name of the kinetic parameter.
- value
  - Type: float
  - Description: Value of the kinetic parameter.
- unit
  - Type: string
  - Description: Unit of the parameter.
- standard_deviation
  - Type: float
  - Description: 1 sigma standard deviation of the kinetic parameter.
- correlations
  - Type: Correlation
  - Multiple: True
  - Descritpion: Correlation of parameter to other parameters of a model.
- upper_limit
  - Type: float
  - Description: Upper limit for parameter value.
- lower_limit
  - Type: float
  - Description: lower limit for parameter value.

### Correlation

- parameter
  - Type: string
  - Description: Name of the parameter.
- value
  - Type: float
  - Description: Correlation value between -1 and 1.

### Series

Time-course data of an individual reaction.

- values
  - Type: float
  - Description: Time-course data of an individual reaction.
  - Multiple: True
- time
  - Type: float
  - Description: Time array corresponding to time-course data.
  - Multiple: True

## Enumerations

### SpeciesTypes

Possible roles of a species in a reaction.

```python
SUBSTRATE = "substrate"
PRODUCT = "product"
INHIBITOR = "inhibitor"
ENZYME = "enzyme"
```

### ConcentrationTypes

```python
MOLAR = "mol / l"
MILLIMOLAR = "mmol / l"
MICROMOLAR = "umol / l"
NANAMOLAR = "nmol / l"
GRAMLITER = "g / l"
MILLIGRAMLITER = "mg / l"
MICROGRAMLITER = "ug / l"
NANOGRAMLITER = "ng / l"
```

### TimeTypes

```python
S = "s"
MIN = "min"
H = "h"
```
