```mermaid
classDiagram
    AbstractSpecies <-- Protein
    AbstractSpecies <-- Complex
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
        +Reactant measured_reactant
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
        +boolean constant*
        +string uri
        +string creator_id
    }
    
    class AbstractSpecies {
        +string name*
        +Vessel vessel_id*
        +float init_conc
        +boolean constant*
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
    
    class ParamType {
        << Enumeration >>
        +K_M
        +K_CAT
        +K_IC
        +K_IU
        +K_IE
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
        +Repository objects=[{'name': 'Vessel', 'docstring': 'This object describes vessels in which the experiment has been carried out. These can include any type of vessel used in biocatalytic experiments.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'Name of the used vessel.', 'template_alias': 'Name'}, {'name': 'volume', 'required': True, 'type': ['posfloat'], 'description': 'Volumetric value of the vessel.', 'template_alias': 'Volume value'}, {'name': 'unit', 'required': True, 'type': ['string'], 'description': 'Volumetric unit of the vessel.', 'template_alias': 'Volume unit'}, {'name': 'constant', 'required': True, 'type': ['boolean'], 'description': 'Whether the volume of the vessel is constant or not.', 'default': 'True'}, {'name': 'uri', 'required': False, 'default': None, 'type': ['string'], 'description': 'URI of the vessel.'}, {'name': 'creator_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unique identifier of the author.'}], 'type': 'object', 'subtypes': [], 'module': 'species'}, {'name': 'AbstractSpecies', 'docstring': 'This object is used to inherit basic attributes common to all species used in the data model.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'None'}, {'name': 'vessel_id', 'required': True, 'reference': 'Vessel.id', 'type': ['Vessel'], 'description': 'None'}, {'name': 'init_conc', 'required': False, 'default': None, 'type': ['float'], 'description': 'None'}, {'name': 'constant', 'required': True, 'type': ['boolean'], 'description': 'None'}, {'name': 'unit', 'required': False, 'default': None, 'type': ['string'], 'description': 'None'}, {'name': 'uri', 'required': False, 'default': None, 'type': ['string'], 'description': 'None'}, {'name': 'creator_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'None'}], 'type': 'object', 'subtypes': [], 'module': 'species'}, {'name': 'Protein', 'docstring': 'This objects describes the proteins that were used or produced in the course of the experiment.', 'attributes': [{'name': 'sequence', 'required': True, 'type': ['string'], 'description': 'Amino acid sequence of the protein', 'template_alias': 'Sequence'}, {'name': 'ecnumber', 'required': False, 'default': None, 'type': ['string'], 'description': 'EC number of the protein.', 'pattern': '(\\d+.)(\\d+.)(\\d+.)(\\d+)', 'template_alias': 'EC Number'}, {'name': 'organism', 'required': False, 'default': None, 'type': ['string'], 'description': 'Organism the protein was expressed in.', 'template_alias': 'Source organism'}, {'name': 'organism_tax_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'Taxonomy identifier of the expression host.'}, {'name': 'uniprotid', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unique identifier referencing a protein entry at UniProt. Use this identifier to initialize the object from the UniProt database.', 'template_alias': 'UniProt ID'}, {'name': 'ontology', 'required': True, 'type': ['SBOTerm'], 'description': 'None', 'default': 'SBOTerm.CATALYST'}], 'type': 'object', 'subtypes': [], 'parent': 'AbstractSpecies', 'module': 'species'}, {'name': 'Reactant', 'docstring': 'This objects describes the reactants that were used or produced in the course of the experiment.', 'attributes': [{'name': 'smiles', 'required': False, 'default': None, 'type': ['string'], 'description': 'Simplified Molecular Input Line Entry System (SMILES) encoding of the reactant.', 'template_alias': 'SMILES'}, {'name': 'inchi', 'required': False, 'default': None, 'type': ['string'], 'description': 'International Chemical Identifier (InChI) encoding of the reactant.', 'template_alias': 'InCHI'}, {'name': 'chebi_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unique identifier of the CHEBI database. Use this identifier to initialize the object from the CHEBI database.'}, {'name': 'ontology', 'required': True, 'type': ['SBOTerm'], 'description': 'None', 'default': 'SBOTerm.SMALL_MOLECULE'}], 'type': 'object', 'subtypes': [], 'parent': 'AbstractSpecies', 'module': 'species'}, {'name': 'Reaction', 'docstring': 'This object describes a chemical or enzymatic reaction that was investigated in the course of the experiment. All species used within this object need to be part of the data model.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'Name of the reaction.', 'template_alias': 'Name'}, {'name': 'reversible', 'required': True, 'type': ['bool'], 'description': 'Whether the reaction is reversible or irreversible', 'default': 'False', 'template_alias': 'Reversible'}, {'name': 'temperature', 'required': False, 'default': None, 'type': ['float'], 'description': 'Numeric value of the temperature of the reaction.', 'template_alias': 'Temperature value'}, {'name': 'temperature_unit', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unit of the temperature of the reaction.', 'pattern': 'kelvin|Kelvin|k|K|celsius|Celsius|C|c', 'template_alias': 'Temperature unit'}, {'name': 'ph', 'required': False, 'default': None, 'type': ['float'], 'description': 'PH value of the reaction.', 'template_alias': 'pH value', 'inclusiveminimum': '0', 'inclusivemaximum': '14'}, {'name': 'ontology', 'required': True, 'type': ['SBOTerm'], 'default': 'SBOTerm.BIOCHEMICAL_REACTION', 'description': 'Ontology defining the role of the given species.'}, {'name': 'uri', 'required': False, 'default': None, 'type': ['string'], 'description': 'URI of the reaction.'}, {'name': 'creator_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unique identifier of the author.'}, {'name': 'model', 'required': False, 'default': None, 'type': ['KineticModel'], 'description': 'Kinetic model decribing the reaction.'}, {'name': 'educts', 'required': False, 'type': ['ReactionElement'], 'default_factory': 'ListPlus()', 'multiple': 'True', 'description': 'List of educts containing ReactionElement objects.', 'template_alias': 'Educts'}, {'name': 'products', 'required': False, 'type': ['ReactionElement'], 'default_factory': 'ListPlus()', 'multiple': 'True', 'description': 'List of products containing ReactionElement objects.', 'template_alias': 'Products'}, {'name': 'modifiers', 'required': False, 'type': ['ReactionElement'], 'default_factory': 'ListPlus()', 'multiple': 'True', 'description': 'List of modifiers (Proteins, snhibitors, stimulators) containing ReactionElement objects.', 'template_alias': 'Modifiers'}], 'type': 'object', 'subtypes': [], 'module': 'enzyme_reaction'}, {'name': 'ReactionElement', 'docstring': 'This object is part of the Reaction object and describes either an educt, product or modifier. The latter includes buffers, counter-ions as well as proteins/enzymes.', 'attributes': [{'name': 'species_id', 'required': True, 'reference': 'AbstractSpecies.id', 'type': ['AbstractSpecies'], 'description': 'Internal identifier to either a protein or reactant defined in the EnzymeMLDocument.', 'references': 'EnzymeMLDocument.reactants.id'}, {'name': 'stoichiometry', 'required': True, 'type': ['posfloat'], 'description': 'Positive float number representing the associated stoichiometry.', 'default': '1.0'}, {'name': 'constant', 'required': True, 'type': ['bool'], 'description': 'Whether or not the concentration of this species remains constant.', 'default': 'False'}, {'name': 'ontology', 'required': False, 'default': None, 'type': ['SBOTerm'], 'description': 'Ontology defining the role of the given species.'}], 'type': 'object', 'subtypes': [], 'module': 'enzyme_reaction'}, {'name': 'KineticModel', 'docstring': 'This object describes a kinetic model that was derived from the experiment.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'Name of the kinetic law.'}, {'name': 'equation', 'required': True, 'type': ['string'], 'description': 'Equation for the kinetic law.'}, {'name': 'parameters', 'required': False, 'type': ['KineticParameter'], 'default_factory': 'ListPlus()', 'multiple': 'True', 'description': 'List of estimated parameters.'}, {'name': 'ontology', 'required': False, 'default': None, 'type': ['SBOTerm'], 'description': 'Type of the estimated parameter.'}], 'type': 'object', 'subtypes': [], 'module': 'modelling'}, {'name': 'KineticParameter', 'docstring': 'This object describes the parameters of the kinetic model and can include all estimated values.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'Name of the estimated parameter.'}, {'name': 'value', 'required': True, 'type': ['float'], 'description': 'Numerical value of the estimated parameter.'}, {'name': 'unit', 'required': True, 'type': ['string'], 'description': 'Unit of the estimated parameter.'}, {'name': 'initial_value', 'required': False, 'default': None, 'type': ['float'], 'description': 'Initial value that was used for the parameter estimation.'}, {'name': 'upper', 'required': False, 'default': None, 'type': ['float'], 'description': 'Upper bound of the estimated parameter.'}, {'name': 'lower', 'required': False, 'default': None, 'type': ['float'], 'description': 'Lower bound of the estimated parameter.'}, {'name': 'is_global', 'required': True, 'type': ['bool'], 'description': 'Specifies if this parameter is a global parameter.', 'default': 'False'}, {'name': 'stdev', 'required': False, 'default': None, 'type': ['float'], 'description': 'Standard deviation of the estimated parameter.'}, {'name': 'constant', 'required': True, 'type': ['bool'], 'description': 'Specifies if this parameter is constant', 'default': 'False'}, {'name': 'ontology', 'required': False, 'default': None, 'type': ['SBOTerm'], 'description': 'Type of the estimated parameter.'}], 'type': 'object', 'subtypes': [], 'module': 'modelling'}, {'name': 'Measurement', 'docstring': 'This object describes the result of a measurement, which includes time course data of any type defined in DataTypes. It includes initial concentrations of all species used in a single measurement.', 'attributes': [{'name': 'name', 'required': True, 'type': ['string'], 'description': 'Name of the measurement'}, {'name': 'temperature', 'required': True, 'type': ['float'], 'description': 'Numeric value of the temperature of the reaction.', 'template_alias': 'Temperature value'}, {'name': 'temperature_unit', 'required': True, 'type': ['string'], 'description': 'Unit of the temperature of the reaction.', 'pattern': 'kelvin|Kelvin|k|K|celsius|Celsius|C|c'}, {'name': 'ph', 'required': True, 'type': ['float'], 'description': 'PH value of the reaction.', 'inclusiveminimum': '0', 'inclusivemaximum': '14'}, {'name': 'species', 'required': False, 'type': ['MeasurementData'], 'default_factory': 'ListPlus()', 'multiple': 'True', 'description': 'Species of the measurement.'}, {'name': 'global_time', 'required': True, 'type': ['float'], 'multiple': 'True', 'description': 'Global time of the measurement all replicates agree on.'}, {'name': 'global_time_unit', 'required': True, 'type': ['string'], 'description': 'Unit of the global time.'}, {'name': 'uri', 'required': False, 'default': None, 'type': ['string'], 'description': 'URI of the reaction.'}, {'name': 'creator_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unique identifier of the author.'}], 'type': 'object', 'subtypes': [], 'module': 'time course data handling'}, {'name': 'MeasurementData', 'docstring': 'This object describes a single entity of a measurement, which corresponds to one species. It also holds replicates which contain time course data.', 'attributes': [{'name': 'init_conc', 'required': True, 'type': ['float'], 'description': 'Initial concentration of the measurement data.'}, {'name': 'unit', 'required': True, 'type': ['string'], 'description': 'The unit of the measurement data.'}, {'name': 'measurement_id', 'required': True, 'type': ['string'], 'description': 'Unique measurement identifier this dataset belongs to.'}, {'name': 'species_id', 'required': False, 'default': None, 'reference': 'AbstractSpecies.id', 'type': ['AbstractSpecies'], 'description': 'The identifier for the described reactant.'}, {'name': 'replicates', 'required': False, 'type': ['Replicate'], 'default_factory': 'ListPlus()', 'multiple': 'True', 'description': 'A list of replicate objects holding raw data of the measurement.'}], 'type': 'object', 'subtypes': [], 'module': 'time course data handling'}, {'name': 'Replicate', 'docstring': 'This object contains the measured time course data as well as metadata to the replicate itself.', 'attributes': [{'name': 'species_id', 'required': True, 'reference': 'AbstractSpecies.id', 'type': ['AbstractSpecies'], 'description': 'Unique identifier of the species that has been measured.'}, {'name': 'measurement_id', 'required': True, 'type': ['string'], 'description': 'Unique identifier of the measurement that the replicate is part of.'}, {'name': 'data_type', 'required': True, 'type': ['DataTypes'], 'description': 'Type of data that was measured (e.g. concentration)', 'default': 'DataTypes.CONCENTRATION'}, {'name': 'data_unit', 'required': True, 'type': ['string'], 'description': 'SI unit of the data that was measured.'}, {'name': 'time_unit', 'required': True, 'type': ['string'], 'description': 'Time unit of the replicate.'}, {'name': 'time', 'required': True, 'type': ['float'], 'multiple': 'True', 'description': 'Time steps of the replicate.'}, {'name': 'data', 'required': True, 'type': ['float'], 'multiple': 'True', 'description': 'Data that was measured.'}, {'name': 'is_calculated', 'required': True, 'type': ['bool'], 'description': 'Whether or not the data has been generated by simulation.', 'default': 'False'}, {'name': 'uri', 'required': False, 'default': None, 'type': ['string'], 'description': 'URI of the protein.'}, {'name': 'creator_id', 'required': False, 'default': None, 'type': ['string'], 'description': 'Unique identifier of the author.'}], 'type': 'object', 'subtypes': [], 'module': 'time course data handling'}] enums=[{'name': 'SBOTerm', 'docstring': 'These are a small fraction of the SBOTerms defined for the SBML markup language.', 'mappings': [{'key': 'BIOCHEMICAL_REACTION', 'value': '"SBO:0000176"'}, {'key': 'ACID_BASE_REACTION', 'value': '"SBO:0000208"'}, {'key': 'CONFORMATIONAL_TRANSITION', 'value': '"SBO:0000181"'}, {'key': 'CONVERSION', 'value': '"SBO:0000182"'}, {'key': 'DEGRADATION', 'value': '"SBO:0000179"'}, {'key': 'DISSOCIATION', 'value': '"SBO:0000180"'}, {'key': 'IONISATION', 'value': '"SBO:0000209"'}, {'key': 'ISOMERISATION', 'value': '"SBO:0000377"'}, {'key': 'NON_COVALENT_BINDING', 'value': '"SBO:0000177"'}, {'key': 'REDOX_REACTION', 'value': '"SBO:0000200"'}, {'key': 'SPONTANEOUS_REACTION', 'value': '"SBO:0000672"'}, {'key': 'PROTEIN', 'value': '"SBO:0000252"'}, {'key': 'GENE', 'value': '"SBO:0000251"'}, {'key': 'SMALL_MOLECULE', 'value': '"SBO:0000247"'}, {'key': 'ION', 'value': '"SBO:0000327"'}, {'key': 'RADICAL', 'value': '"SBO:0000328"'}, {'key': 'INTERACTOR', 'value': '"SBO:0000336"'}, {'key': 'SUBSTRATE', 'value': '"SBO:0000015"'}, {'key': 'PRODUCT', 'value': '"SBO:0000011"'}, {'key': 'CATALYST', 'value': '"SBO:0000013"'}, {'key': 'INHIBITOR', 'value': '"SBO:0000020"'}, {'key': 'ESSENTIAL_ACTIVATOR', 'value': '"SBO:0000461"'}, {'key': 'NON_ESSENTIAL_ACTIVATOR', 'value': '"SBO:0000462"'}, {'key': 'POTENTIATOR', 'value': '"SBO:0000021"'}, {'key': 'MACROMOLECULAR_COMPLEX', 'value': '"SBO:0000296"'}, {'key': 'PROTEIN_COMPLEX', 'value': '"SBO:0000297"'}, {'key': 'DIMER', 'value': '"SBO:0000607"'}, {'key': 'MICHAELIS_MENTEN', 'value': '"SBO:0000028"'}, {'key': 'K_CAT', 'value': '"SBO:0000025"'}, {'key': 'K_M', 'value': '"SBO:0000027"'}, {'key': 'V_MAX', 'value': '"SBO:0000186"'}], 'type': 'enum'}, {'name': 'DataTypes', 'docstring': 'These values are used to determine the type of time course data.', 'mappings': [{'key': 'CONCENTRATION', 'value': '"conc"'}, {'key': 'ABSORPTION', 'value': '"abs"'}, {'key': 'FEED', 'value': '"feed"'}, {'key': 'BIOMASS', 'value': '"biomass"'}, {'key': 'CONVERSION', 'value': '"conversion"'}, {'key': 'PEAK_AREA', 'value': '"peak-area"'}], 'type': 'enum'}] inherits=[{'parent': 'AbstractSpecies', 'child': 'Protein'}, {'parent': 'AbstractSpecies', 'child': 'Complex'}, {'parent': 'AbstractSpecies', 'child': 'Reactant'}] compositions=[{'container': 'EnzymeMLDocument', 'module': 'Creator'}, {'container': 'EnzymeMLDocument', 'module': 'Vessel'}, {'container': 'EnzymeMLDocument', 'module': 'Protein'}, {'container': 'EnzymeMLDocument', 'module': 'Complex'}, {'container': 'EnzymeMLDocument', 'module': 'Reactant'}, {'container': 'EnzymeMLDocument', 'module': 'Reaction'}, {'container': 'EnzymeMLDocument', 'module': 'KineticParameter'}, {'container': 'EnzymeMLDocument', 'module': 'Measurement'}, {'container': 'EnzymeMLDocument', 'module': 'File'}, {'container': 'AbstractSpecies', 'module': 'Vessel'}, {'container': 'Protein', 'module': 'SBOTerm'}, {'container': 'Complex', 'module': 'SBOTerm'}, {'container': 'Reactant', 'module': 'SBOTerm'}, {'container': 'Reaction', 'module': 'SBOTerm'}, {'container': 'Reaction', 'module': 'ReactionElement'}, {'container': 'Reaction', 'module': 'KineticModel'}, {'container': 'ReactionElement', 'module': 'SBOTerm'}, {'container': 'ReactionElement', 'module': 'AbstractSpecies'}, {'container': 'KineticModel', 'module': 'SBOTerm'}, {'container': 'KineticModel', 'module': 'KineticParameter'}, {'container': 'KineticParameter', 'module': 'SBOTerm'}, {'container': 'Measurement', 'module': 'MeasurementData'}, {'container': 'MeasurementData', 'module': 'AbstractSpecies'}, {'container': 'MeasurementData', 'module': 'Replicate'}, {'container': 'Replicate', 'module': 'DataTypes'}, {'container': 'Replicate', 'module': 'AbstractSpecies'}] external_objects={}
    }
    
```