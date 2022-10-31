```mermaid
classDiagram
    EnzymeKinetics *-- Measurement
    EnzymeKinetics *-- StoichiometryTypes
    EnzymeKinetics *-- ConcentrationTypes
    EnzymeKinetics *-- TimeTypes
    Measurement *-- Series
    
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
        +float initial_substrate_conc*
        +float[0..*] enzyme_conc*
        +float inhibitor_conc
        +str inhibitor_conc_unit
        +Series[0..*] data*
        +float[0..*] values
    }
    
    class Series {
        << Enumeration >>
    }
    
    class StoichiometryTypes {
        << Enumeration >>
        +SUBSTRATE = "substrate"
        +PRODUCT = "product"
    }
    
    class TimeTypes {
        << Enumeration >>
        +S = "s"
        +MIN = "min"
        +H = "h"
    }
    
    class ConcentrationTypes {
        << Enumeration >>
        +M = "mole / l"
        +mM = "mmole / l"
        +uM = "umole / l"
        +nM = "nmole / l"
    }
    
```