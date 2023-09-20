# EnzymePynetics

## Objects

### Estimator

- name
    - Type: string
    - Description: Title of the kinetic experiment
- measured_reactant: 
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@Reactant
    - Description: Reactant that is measured in the experiment
- reaction_systems
    - Type: ReactionSystem
    - Description: Reactions of multiple species
    - Multiple: True
- species
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@AbstractSpecies, https://github.com/EnzymeML/enzymeml-specifications.git@Protein, https://github.com/EnzymeML/enzymeml-specifications.git@Reactant
    - Description: Reactants, Inhibitor, Activators and Catalysts of the reaction
    - Multiple: True
- reaction
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@Reaction
    - Description: Reaction proceeding in measurements
- models
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@KineticModel
    - Description: Kinetic model options used to fit to measurement data.
    - Multiple: True
- measurements
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@Measurement
    - Description: Measurement data for a given initial substrate concentration.
    - Multiple: True

### ReactionSystem

- name
    - Type: string
    - Description: Name of the reaction system
- reactions
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@Reaction
    - Description: Reactions of the reaction system
    - Multiple: True
- result
    - Type: ModelResult
    - Description: Result of the kinetic model fitting.


### ModelResult

Description of a kinetic model

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
- correlations
  - Type: Correlation
  - Multiple: True
  - Description: Correlation of parameter to other parameters of a model.

### Correlation

- parameter_name
  - Type: string
  - Description: Name of the parameter.
- value
  - Type: float
  - Description: Correlation value between -1 and 1.