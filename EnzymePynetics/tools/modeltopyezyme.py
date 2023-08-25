from typing import Dict, List
import pyenzyme as pe
from pyenzyme.enzymeml.core.ontology import SBOTerm
from pyenzyme.enzymeml.core.enzymereaction import EnzymeReaction
from pyenzyme.enzymeml.models import KineticModel as KM
from pyenzyme.enzymeml.models import KineticParameter

from EnzymePynetics.tools.kineticmodel import KineticModel


def _add_results(estimator: "ParameterEstimator", kinetic_model: KineticModel, enzymeml: pe.EnzymeMLDocument):

    # Check if EnzymeML has reaction
    if not enzymeml.getReaction("r0"):
        reaction = _create_substrate_reaction(kinetic_model, enzymeml, e)
        enzymeml.addReaction(reaction)
    else:
        substrate_reaction = enzymeml.reaction_dict["r0"]

    # Add kinetic model and parameters
    species_mapping = _get_ontology_id_dict(estimator, enzymeml)

    _add_model_results(
        kinetic_model=kinetic_model,
        reaction=substrate_reaction,
        species_mapping=species_mapping,
        equation=kinetic_model.substrate_rate_law)

    if kinetic_model.enzyme_rate_law:
        _create_enzyme_reaction(
            kinetic_model=kinetic_model,
            enzymeml=enzymeml,
            species_mapping=species_mapping
        )

    return enzymeml


def _create_enzyme_reaction(kinetic_model: KineticModel, enzymeml: pe.EnzymeMLDocument, species_mapping: dict):
    inactive_protein_id = _add_inactive_protein_species(
        enzymeml=enzymeml,
        active_protein_id=species_mapping["enzyme"],
    )

    enzyme_inactivation = _create_inactivation_reaction(
        enzymeml=enzymeml,
        active_protein_id=species_mapping["enzyme"],
        inactive_protein_id=inactive_protein_id,
    )

    _add_model_results(
        kinetic_model=kinetic_model,
        reaction=enzyme_inactivation,
        species_mapping=species_mapping,
        equation=kinetic_model.enzyme_rate_law
    )


def _add_inactive_protein_species(enzymeml: pe.EnzymeMLDocument, active_protein_id: str) -> str:

    active_protein = enzymeml.getProtein(active_protein_id)
    inactive_protein = pe.Protein(**active_protein.dict())

    # Modify protein
    inactive_protein.name += " inactive"
    inactive_protein.init_conc = None
    inactive_protein.constant = False
    inactive_protein.id = "p1"

    return enzymeml.addProtein(inactive_protein)


def _create_inactivation_reaction(enzymeml: pe.EnzymeMLDocument, active_protein_id: str, inactive_protein_id: str) -> EnzymeReaction:

    # Create new reaction
    reaction_params = enzymeml.getReaction("r0")

    enzyme_inactivation = pe.EnzymeReaction(
        name="Time-dependent enzyme inactivation",
        reversible=False,
        temperature=reaction_params.temperature,
        temperature_unit=reaction_params.temperature_unit,
        ph=reaction_params.ph
    )

    # Add reaction species
    enzyme_inactivation.addEduct(active_protein_id, 1.0, enzymeml)
    enzyme_inactivation.addProduct("p1", 1.0, enzymeml)
    enzyme_inactivation_id = enzymeml.addReaction(enzyme_inactivation)

    return enzymeml.getReaction(enzyme_inactivation_id)


def _create_substrate_reaction(estimator: "ParameterEstimator", kinetic_model: KineticModel, enzymeml: pe.EnzymeMLDocument):

    # Get reaction conditions
    measurement = next(iter(enzymeml.measurement_dict.values()))
    ph = measurement.ph
    temperature = measurement.temperature
    temperature_unit = measurement.temperature_unit

    substate = estimator.substrate_name
    substrate_id = _get_species_id(enzymeml.reactant_dict, substate)
    product = estimator.product_name
    product_id = _get_species_id(enzymeml.reactant_dict, product)
    enzyme = estimator.enzyme_name
    enzyme_id = _get_species_id(enzymeml.protein_dict, enzyme)
    inhibitor = estimator.inhibitor_name
    inhibitor_id = _get_species_id(enzymeml.reactant_dict, inhibitor)

    if kinetic_model.enzyme_rate_law:
        enzyme_is_constant = False
    else:
        enzyme_is_constant = True

    model = KM(
        name=kinetic_model.name,
        eqautions=kinetic_model.substrate_rate_law,
    )
    for param in kinetic_model.result.parameters:
        model.addParameter(
            name=param.name,
            value=param.value,
            unit=_fix_mole(param.unit),
            upper=param.upper_limit,
            lower=param.lower_limit,
            stdev=param.standard_deviation,
        )

    # Create reaction
    reaction = pe.EnzymeReaction(
        name="substrate reaction",
        temperature=temperature,
        temperature_unit=temperature_unit,
        ph=ph,
        reversible=False,
        ontology=SBOTerm.BIOCHEMICAL_REACTION,
        model=model,
    )

    # Add reaction elements
    reaction.addEduct(
        species_id=substrate_id,
        stoichiometry=1,
        constant=False,
        ontology=SBOTerm.SUBSTRATE,
    )

    reaction.addProduct(
        species_id=product_id,
        stoichiometry=1,
        constant=False,
        ontology=SBOTerm.PRODUCT,
    )

    reaction.addModifier(
        species_id=enzyme_id,
        stoichiometry=1,
        constant=enzyme_is_constant,
        ontology=SBOTerm.CATALYST,
    )

    if inhibitor_id:
        reaction.addModifier(
            species_id=inhibitor_id,
            stoichiometry=1,
            constant=False,
            ontology=SBOTerm.INHIBITOR,
        )


def _add_model_results(
        kinetic_model: KineticModel,
        reaction: EnzymeReaction,
        species_mapping: dict,
        equation: str
) -> EnzymeReaction:

    equation = _map_ids_to_equation(equation, species_mapping)

    params = [
        KineticParameter(
            name=param.name,
            value=param.value,
            unit=_fix_mole(param.unit),
            upper=param.upper_limit,
            lower=param.lower_limit,
            stdev=param.standard_deviation,
            constant=True
        )
        for param in kinetic_model.result.parameters
        if param.name in equation
    ]

    reaction.model = KM(
        name=kinetic_model.name,
        parameters=params,
        equation=equation
    )


def _get_species_id(species_dict: dict, species_name: str):

    for species_id, species in species_dict.items():
        if species.name == species_name:
            return species_id

    else:
        return ""


def _get_ontology_id_dict(estimator: "ParameterEstimator", enymeml: pe.EnzymeMLDocument) -> Dict[str, str]:

    reactant_dict = enymeml.reactant_dict
    protein_dict = enymeml.protein_dict

    return dict(
        substrate=_get_species_id(reactant_dict, estimator.substrate_name),
        product=_get_species_id(reactant_dict, estimator.product_name),
        enzyme=_get_species_id(protein_dict, estimator.enzyme_name),
        inhibitor=_get_species_id(reactant_dict, estimator.inhibitor_name),
    )


def _map_ids_to_equation(equation: str, ontology_id_dict: Dict[str, str]) -> str:

    for key, value in ontology_id_dict.items():
        equation = equation.replace(key, value)

    return equation


def _fix_mole(string: str) -> str:
    if "mol" in string and "mole" not in string:
        return string.replace("mol", "mole")
    else:
        return string
