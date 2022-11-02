```mermaid
classDiagram
    EnzymeKinetics *-- Measurement
    EnzymeKinetics *-- StoichiometryTypes
    EnzymeKinetics *-- ConcentrationTypes
    EnzymeKinetics *-- TimeTypes
    
    class EnzymeKinetics {
        +string title
        +string reactant_name
        +Measurement[0..*] measurements*
        +StoichiometryTypes stoichiometry
        +ConcentrationTypes data_conc_unit*
        +float[0..*] time*
        +TimeTypes time_unit*
    }
    
    class Measurement {
        +float initial_conc
    }
    
    class StoichiometryTypes {
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
    }
    
    class TimeTypes {
        << Enumeration >>
        +S = "s"
        +MIN = "min"
        +H = "h"
    }
    
```