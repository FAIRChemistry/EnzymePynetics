from typing import Tuple
from sdRDM import DataModel
from EnzymePynetics.modified.reactionsystem import ReactionSystem
from EnzymePynetics.modified.protein import Protein
from EnzymePynetics.modified.reactant import Reactant

# Specify EnzymeML version
URL = "https://github.com/EnzymeML/enzymeml-specifications.git"
COMMIT = "72c3d8be4a094983667a7aa62fb599fbc9f7351c"

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
        with open(path, "w") as f:
            f.write(enzymeml.json())

    return enzymeml
