```mermaid
classDiagram
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
        +float[0..*] time
        +TimeTypes time_unit
    }
    
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
    }
    
    class Series {
        +float[0..*] values
    }
    
    class ReactantTypes {
        << Enumeration >>
        +SUBSTRATE = "substrate"
        +PRODUCT = "product"
    }
    
    class ReactantTypes {
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
        +NANOGRAMLITER = "ng / l"
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
        +NANOGRAMLITER = "ng / l"
    }
    
    class TimeTypes {
        << Enumeration >>
        +S = "s"
        +MIN = "min"
        +H = "h"
    }
    
```