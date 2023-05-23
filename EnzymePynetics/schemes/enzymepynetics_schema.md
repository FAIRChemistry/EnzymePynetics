```mermaid
classDiagram
    EnzymeKinetics *-- Measurement
    EnzymeKinetics *-- ModelResult
    Species *-- SpeciesTypes
    Species *-- ConcentrationTypes
    Species *-- Series
    Measurement *-- TimeTypes
    Measurement *-- Species
    ModelResult *-- Parameter
    
    class EnzymeKinetics {
        +string title
        +ModelResult[0..*] model_results
        +Measurement[0..*] measurements
    }
    
    class Species {
        +string name
        +ConcentrationTypes conc_unit
        +float initial_conc
        +SpeciesTypes species_type
        +Series[0..*] data
    }
    
    class Measurement {
        +Species[0..*] species
        +float temperature
        +string temperature_unit
        +float pH
        +float[0..*] time
        +TimeTypes time_unit
    }
    
    class ModelResult {
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
        +float upper_limit
        +float lower_limit
    }
    
    class Series {
        +float[0..*] values
    }
    
    class SpeciesTypes {
        << Enumeration >>
        +SUBSTRATE
        +PRODUCT
        +INHIBITOR
        +ENZYME
    }
    
    class ConcentrationTypes {
        << Enumeration >>
        +MOLAR
        +MILLIMOLAR
        +MICROMOLAR
        +NANAMOLAR
        +GRAMLITER
        +MILLIGRAMLITER
        +MICROGRAMLITER
        +NANOGRAMLITER
    }
    
    class TimeTypes {
        << Enumeration >>
        +S
        +MIN
        +H
    }
    
```