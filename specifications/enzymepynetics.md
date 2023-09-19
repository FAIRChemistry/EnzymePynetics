# EnzymePynetics

## Objects

### Estimator

- name
    - Type: string
    - Description: Title of the kinetic experiment
- measured_reactant: 
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@Reactant
    - Description: Reactant that is measured in the experiment
- species
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@AbstractSpecies, https://github.com/EnzymeML/enzymeml-specifications.git@Protein, https://github.com/EnzymeML/enzymeml-specifications.git@Reactant
    - Description: Reactants, Inhibitor, Activators and Catalysts of the reaction
    - Multiple: True
- reactions
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@Reaction
    - Description: Reaction proceeding in measurements
    - Multiple: True
- models
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@KineticModel
    - Description: Kinetic model options used to fit to measurement data.
    - Multiple: True
- measurements
    - Type: https://github.com/EnzymeML/enzymeml-specifications.git@Measurement
    - Description: Measurement data for a given initial substrate concentration.
    - Multiple: True