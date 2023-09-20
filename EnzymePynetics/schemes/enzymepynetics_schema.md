```mermaid
classDiagram
    AbstractSpecies <-- Protein
    AbstractSpecies <-- Complex
    AbstractSpecies <-- Reactant
    AbstractSpecies <-- Protein
    AbstractSpecies <-- Reactant
    EnzymeMLDocument *-- Creator
    EnzymeMLDocument *-- Vessel
    EnzymeMLDocument *-- Protein
    EnzymeMLDocument *-- Complex
    EnzymeMLDocument *-- Reactant
    EnzymeMLDocument *-- Reaction
    EnzymeMLDocument *-- KineticParameter
    EnzymeMLDocument *-- Measurement
    EnzymeMLDocument *-- File
    AbstractSpecies *-- Vessel
    Protein *-- SBOTerm
    Complex *-- SBOTerm
    Reactant *-- SBOTerm
    Reaction *-- SBOTerm
    Reaction *-- ReactionElement
    Reaction *-- KineticModel
    ReactionElement *-- SBOTerm
    ReactionElement *-- AbstractSpecies
    KineticModel *-- SBOTerm
    KineticModel *-- KineticParameter
    KineticParameter *-- SBOTerm
    Measurement *-- MeasurementData
    MeasurementData *-- AbstractSpecies
    MeasurementData *-- Replicate
    Replicate *-- DataTypes
    Replicate *-- AbstractSpecies
    Estimator *-- ReactionSystem
    Estimator *-- AbstractSpecies
    Estimator *-- Protein
    Estimator *-- Reactant
    Estimator *-- Reaction
    Estimator *-- KineticModel
    Estimator *-- Measurement
    ReactionSystem *-- ModelResult
    ReactionSystem *-- Reaction
    ModelResult *-- Parameter
    Parameter *-- Correlation
    AbstractSpecies *-- Vessel
    Protein *-- SBOTerm
    Reactant *-- SBOTerm
    Reaction *-- SBOTerm
    Reaction *-- ReactionElement
    Reaction *-- KineticModel
    ReactionElement *-- SBOTerm
    ReactionElement *-- AbstractSpecies
    KineticModel *-- SBOTerm
    KineticModel *-- KineticParameter
    KineticParameter *-- SBOTerm
    Measurement *-- MeasurementData
    MeasurementData *-- AbstractSpecies
    MeasurementData *-- Replicate
    Replicate *-- DataTypes
    Replicate *-- AbstractSpecies
    
    class Estimator {
        +string name
        +Reactant measured_reactant:
        +ReactionSystem[0..*] reaction_systems
        +AbstractSpecies, Protein, Reactant[0..*] species
        +Reaction reaction
        +KineticModel[0..*] models
        +Measurement[0..*] measurements
    }
    
    class ReactionSystem {
        +string name
        +Reaction[0..*] reactions
        +ModelResult result
    }
    
    class ModelResult {
        +string[0..*] equations
        +Parameter[0..*] parameters
        +bool fit_success
        +float AIC
        +float BIC
        +float RMSD
    }
    
    class Parameter {
        +string name
        +Correlation[0..*] correlations
    }
    
    class Correlation {
        +string parameter_name
        +float value
    }
    
    class Vessel {
        +string name*
        +posfloat volume*
        +string unit*
        +StrictBool constant*
        +string uri
        +string creator_id
    }
    
    class AbstractSpecies {
        +string name*
        +Vessel vessel_id*
        +float init_conc
        +StrictBool constant*
        +string unit
        +string uri
        +string creator_id
    }
    
    class Protein {
        +string sequence*
        +string ecnumber
        +string organism
        +string organism_tax_id
        +string uniprotid
        +SBOTerm ontology*
    }
    
    class Reactant {
        +string smiles
        +string inchi
        +string chebi_id
        +SBOTerm ontology*
    }
    
    class Reaction {
        +string name*
        +bool reversible*
        +float temperature
        +string temperature_unit
        +float ph
        +SBOTerm ontology*
        +string uri
        +string creator_id
        +KineticModel model
        +ReactionElement[0..*] educts
        +ReactionElement[0..*] products
        +ReactionElement[0..*] modifiers
    }
    
    class ReactionElement {
        +AbstractSpecies species_id*
        +posfloat stoichiometry*
        +bool constant*
        +SBOTerm ontology
    }
    
    class KineticModel {
        +string name*
        +string equation*
        +KineticParameter[0..*] parameters
        +SBOTerm ontology
    }
    
    class KineticParameter {
        +string name*
        +float value*
        +string unit*
        +float initial_value
        +float upper
        +float lower
        +bool is_global*
        +float stdev
        +bool constant*
        +SBOTerm ontology
    }
    
    class Measurement {
        +string name*
        +float temperature*
        +string temperature_unit*
        +float ph*
        +MeasurementData[0..*] species
        +float[0..*] global_time*
        +string global_time_unit*
        +string uri
        +string creator_id
    }
    
    class MeasurementData {
        +float init_conc*
        +string unit*
        +string measurement_id*
        +AbstractSpecies species_id
        +Replicate[0..*] replicates
    }
    
    class Replicate {
        +AbstractSpecies species_id*
        +string measurement_id*
        +DataTypes data_type*
        +string data_unit*
        +string time_unit*
        +float[0..*] time*
        +float[0..*] data*
        +bool is_calculated*
        +string uri
        +string creator_id
    }
    
    class SBOTerm {
        << Enumeration >>
        +BIOCHEMICAL_REACTION
        +ACID_BASE_REACTION
        +CONFORMATIONAL_TRANSITION
        +CONVERSION
        +DEGRADATION
        +DISSOCIATION
        +IONISATION
        +ISOMERISATION
        +NON_COVALENT_BINDING
        +REDOX_REACTION
        +SPONTANEOUS_REACTION
        +PROTEIN
        +GENE
        +SMALL_MOLECULE
        +ION
        +RADICAL
        +INTERACTOR
        +SUBSTRATE
        +PRODUCT
        +CATALYST
        +INHIBITOR
        +ESSENTIAL_ACTIVATOR
        +NON_ESSENTIAL_ACTIVATOR
        +POTENTIATOR
        +MACROMOLECULAR_COMPLEX
        +PROTEIN_COMPLEX
        +DIMER
        +MICHAELIS_MENTEN
        +K_CAT
        +K_M
        +V_MAX
    }
    
    class DataTypes {
        << Enumeration >>
        +CONCENTRATION
        +ABSORPTION
        +FEED
        +BIOMASS
        +CONVERSION
        +PEAK_AREA
    }
    
    class https://github.com/EnzymeML/enzymeml-specifications.git {
        << External Object >>
        +Repository <sdRDM.markdown.markdownparser.MarkdownParser object at 0x140a70b10>
    }
    
```