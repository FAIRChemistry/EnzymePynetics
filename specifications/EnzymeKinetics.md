# EnzymePanetics data model

The following data model describes (meta-)data from an enzyme kinetics experiment in a structured way. It constists out of multiple ```measurements```, describing one or multiple measurements at different initial substrate and/or enzyme concentrations. If and inhibitor was applied to the reaction, its concentration can be speciefied as well.

### EnzymeKinetics

Base class, dealing with measurement data of an enzyme kinetics assay.

- __title__
  - Type: string
  - Description: Title of the kinetic experiment.
- __kinetic_models__
  - Type: KineticModel
  - Description: Kinetic moodels which were used for parameter estimation.
  - Multiple: True
- __measurements*__
  - Type: Measurement
  - Description: Measurement data for a given initial substrate concentration.
  - Multiple: True

### Species

- __name*__
  - Type: string
  - Description: name of the reactant. 
- __conc_unit*__
  - Type: ConcentrationTypes
  - Description: Concentration unit of the measurement data.
- __initial_conc*__
  - Type: float
  - Description: Initial concentration of the reactant.
- __reactant_type__
  - Type: ReactantTypes
  - Description: Define whether "substrate" or "product" concentration was measured.
- __data__
  - Type: Series
  - Description: One or multiple time-course measurement data arrays.
  - Multiple: True

### Measurement

A Measurement object contains information about the applied enzyme concentration and one or multiple time-course concentration measurements. Additionally, the initial substrate concentration should be specified. This is neccessary to derive the substrate concentration for the modeling process. If an inhibitor was applied to the measurement, its concentration and the respective conetration unit can be specified to account for inhibition in kinetic modeling.

- __species__
  - Type: Species
  - Description: Reactants of the reaction.
  - Multiple: True
- __enzyme_conc__
  - Type: float
  - Description: Enzyme concentration in the reaction.
- __temperature__
  - Type: float
  - Description: Temperature of the reaction.
- __temperature_unit__
  - Type: string
  - Description: Temperature unit.
- __pH__
  - Type: float
  - Description: pH of the reaction
- __time*__
  - Type: float
  - Description: Time array corresponding to time-course data.
  - Multiple: True
- __time_unit*__
  - Type: TimeTypes
  - Description: Time data unit.

### KineticModel

Description of a kinetic model

- __name__
  - Type: string
  - Description: Name of the kinetic model.
- __equation__
  - Type: string
  - Description: Equation of the kinetic model.
- __parameters__
  - Type: Parameter
  - Description: Kinetic parameters of the model.
  - Multiple: True
- __AIC__
  - Type: float
  - Description: Akaike information criterion.
- __BIC__
  - Type: float
  - Description: Bayesian information criterion.
- __RMSD__
  - Type: float
  - Description: Root mean square deviation between model and measurement data.

### Parameter

Defines a kinetic parameter.

- __name__
  - Type: string
  - Description: Name of the kinetic parameter
- __value__
  - Type: float
  - Description: Value of the kinetic parameter.
- __standard_deviation__
  - Type: float
  - Description: Standard deviation of the kinetic parameter.

### Series

Time-course data of an individual reaction.

- __values__
  - Type: float
  - Description: Time-course data of an individual reaction.
  - Multiple: True

#### ReactantTypes

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
NANOGRAMLITER = "ng / l"
```

#### TimeTypes

```python
S = "s"
MIN = "min"
H = "h"
```
