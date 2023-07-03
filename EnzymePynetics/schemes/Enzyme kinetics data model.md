```mermaid
classDiagram
<<<<<<< HEAD
    EnzymeKineticsExperiment *-- Measurement
    EnzymeKineticsExperiment *-- StoichiometryTypes
    EnzymeKineticsExperiment *-- ConcentrationTypes
    EnzymeKineticsExperiment *-- TimeTypes
    Measurement *-- Series
    Measurement *-- ConcentrationTypes
    
    class EnzymeKineticsExperiment {
        +string title
        +string reactant_name
        +float temperature
        +string temperature_unit
        +float pH
        +Measurement[0..*] measurements*
        +StoichiometryTypes stoichiometry
        +ConcentrationTypes data_conc_unit*
=======
    EnzymeKinetics *-- KineticModel
    EnzymeKinetics *-- Measurement
    Species *-- ConcentrationTypes
    Species *-- ReactantTypes
    Species *-- Series
    Measurement *-- Species
    Measurement *-- TimeTypes
    KineticModel *-- Parameter
    
    class EnzymeKinetics {
        +string title
        +KineticModel[0..*] kinetic_models
        +Measurement[0..*] measurements*
    }
    
    class Species {
        +string name*
        +ConcentrationTypes conc_unit*
        +float initial_conc*
        +ReactantTypes reactant_type
        +Series[0..*] data
    }
    
    class Measurement {
        +Species[0..*] species
        +float enzyme_conc
        +float temperature
        +string temperature_unit
        +float pH
>>>>>>> dev
        +float[0..*] time*
        +TimeTypes time_unit*
    }
    
<<<<<<< HEAD
    class Measurement {
        +float initial_substrate_conc*
        +float enzyme_conc
        +Series[0..*] data*
        +float inhibitor_conc
        +ConcentrationTypes inhibitor_conc_unit
=======
    class KineticModel {
        +string name
        +string equation
        +Parameter[0..*] parameters
        +float AIC
        +float BIC
        +float RMSD
    }
    
    class Parameter {
        +string name
        +float value
        +float standard_deviation
>>>>>>> dev
    }
    
    class Series {
        +float[0..*] values
<<<<<<< HEAD
        +string test
    }
    
    class StoichiometryTypes {
=======
    }
    
    class ReactantTypes {
>>>>>>> dev
        << Enumeration >>
        +SUBSTRATE = "substrate"
        +PRODUCT = "product"
    }
    
<<<<<<< HEAD
    class StoichiometryTypes {
=======
    class ReactantTypes {
>>>>>>> dev
        << Enumeration >>
        +SUBSTRATE = "substrate"
        +PRODUCT = "product"
    }
    
    class ConcentrationTypes {
        << Enumeration >>
        +MOLAR = "mole / l"
        +MILLIMOLAR = "mmole / l"
        +MICROMOLAR = "umole / l"
        +NANAMOLAR = "nmole / l"
        +GRAMLITER = "g / l"
        +MILLIGRAMLITER = "mg / l"
        +MICROGRAMLITER = "ug / l"
        +NANGRAMLITER = "ng / l"
    }
    
    class ConcentrationTypes {
        << Enumeration >>
        +MOLAR = "mole / l"
        +MILLIMOLAR = "mmole / l"
        +MICROMOLAR = "umole / l"
        +NANAMOLAR = "nmole / l"
        +GRAMLITER = "g / l"
        +MILLIGRAMLITER = "mg / l"
        +MICROGRAMLITER = "ug / l"
        +NANGRAMLITER = "ng / l"
    }
    
    class TimeTypes {
        << Enumeration >>
        +S = "s"
        +MIN = "min"
        +H = "h"
    }
    
```