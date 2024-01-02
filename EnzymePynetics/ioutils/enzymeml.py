import copy
import pyenzyme as pe
from typing import Tuple
from sdRDM import DataModel
from EnzymePynetics.core.reactionsystem import ReactionSystem
from EnzymePynetics.core.protein import Protein
from EnzymePynetics.core.reactant import Reactant

# Specify EnzymeML version
URL = "https://github.com/EnzymeML/enzymeml-specifications.git"
COMMIT = "5e5f05b9dc76134305b8f9cef65271e35563ac76"

EnzymeML = DataModel.from_git(URL, COMMIT)
SBOTerm = EnzymeML.enums.SBOTerm
DataTypes = EnzymeML.enums.DataTypes


def parse_enzymeml(
    cls: "Estimator", enzymeml: "EnzymeML.EnzymeMLDocument", measured_reactant: Reactant
) -> Tuple["Estimator", "EnzymeML.EnzymeMLDocument"]:
    if isinstance(enzymeml, str) and isinstance(measured_reactant, str):
        enzymeml, _ = DataModel.parse(enzymeml)
        measured_reactant = get_measured_species(enzymeml, measured_reactant)

    species = []
    for reactant in enzymeml.reactants:
        species.append(Reactant(**reactant.to_dict()))

    for protein in enzymeml.proteins:
        species.append(Protein(**protein.to_dict()))

    return (
        cls(
            name=enzymeml.name,
            measured_reactant=measured_reactant,
            measurements=enzymeml.measurements,
            species=species,
        ),
        enzymeml,
    )


def get_measured_species(
    enzymeml: "EnzymeML.EnzymeMLDocument", measured_reactant: str
) -> Reactant:
    """Checks if 'measured_reactant' is a valid reactant.name in
    the 'EnzymeML.EnzymeMLDocument'.

    Args:
        enzymeml (EnzymeML.EnzymeMLDocument): EnzymeML document
        measured_reactant (str): Reactant measured in the experiment

    Raises:
        ValueError: If Reactant not found in EnzymeML document

    Returns:
        Reactant: Reactant object
    """

    for species in enzymeml.reactants:
        if species.name.lower().strip() == measured_reactant.lower().strip():
            return species

    raise ValueError(
        f"'{measured_reactant}' not found in EnzymeML document.",
        f"Available reactants: {[reactant.name for reactant in enzymeml.reactants]}",
    )


def _to_enzymeml(
    enzymeml: "EnzymeML.EnzymeMLDocument",
    reaction_system: ReactionSystem,
    out_path: str,
) -> "EnzymeML.EnzymeMLDocument":
    if isinstance(enzymeml, str):
        out_path = enzymeml
        enzymeml, _ = DataModel.parse(enzymeml)

    for reaction in reaction_system.reactions:
        enzymeml.reactions.append(reaction)

    if out_path:
        with open(out_path, "w") as f:
            f.write(enzymeml.json())

    return enzymeml


def map_to_pyenzyme(enzymeml: "EnzymeML.EnzymeMLDocument") -> pe.EnzymeMLDocument:
    doc = pe.EnzymeMLDocument(name=enzymeml.name)

    for vessel in enzymeml.vessels:
        doc.addVessel(pe.Vessel(**vessel.to_dict()))

    for creator in enzymeml.creators:
        doc.addCreator(pe.Creator(**creator.to_dict()))

    for reactant in enzymeml.reactants:
        doc.addReactant(pe.Reactant(**reactant.to_dict()))

    for protein in enzymeml.proteins:
        doc.addProtein(pe.Protein(**protein.to_dict()))

    for measurement in enzymeml.measurements:
        measurement_dict = copy.deepcopy(measurement.to_dict())
        pyenz_measurement = pe.Measurement(**measurement_dict)
        for species in measurement.species:
            species_dict = copy.deepcopy(species.to_dict())
            species_dict.pop("id")
            species_dict.pop("measurement_id")
            species_dict.pop("__source__")
            species_id = species_dict.pop("species_id")
            change_mol_to_mole(species_dict)

            if species_id in doc.reactant_dict.keys():
                species_dict["reactant_id"] = species_id
            elif species_id in doc.protein_dict.keys():
                species_dict["protein_id"] = species_id
            else:
                raise ValueError(
                    f"Species with id '{species_id}' not found in EnzymeML document."
                )
            pyenz_measurement.addData(**species_dict)

        doc.addMeasurement(pyenz_measurement)

    for reaction in enzymeml.reactions:
        reaction = copy.deepcopy(reaction.to_dict())
        change_mol_to_mole(reaction)
        doc.addReaction(pe.EnzymeReaction(**reaction))

    return doc


def change_mol_to_mole(data):
    if isinstance(data, dict):
        for key, value in data.items():
            if key == "unit" or key == "data_unit":
                data[key] = value.replace("mol", "mole")
            elif isinstance(value, list):
                for item in value:
                    change_mol_to_mole(item)
            elif isinstance(value, dict):
                change_mol_to_mole(value)
            elif isinstance(data, list):
                for item in data:
                    change_mol_to_mole(item)


def _to_omex(
    enzymeml: "EnzymeML.EnzymeMLDocument",
    reaction_system: ReactionSystem,
    out_path: str,
) -> pe.EnzymeMLDocument:
    enzml = _to_enzymeml(enzymeml, reaction_system, out_path=None)
    doc = map_to_pyenzyme(enzml)
    doc.toFile(
        out_path, name=f"{enzymeml.name} {reaction_system.reactions[0].temperature}C"
    )
    return doc


if __name__ == "__main__":
    path = "/Users/max/Documents/GitHub/kinetic_modeling_workflow/out.json"
    doc, _ = DataModel.parse(path)
    map_to_pyenzyme(doc)
