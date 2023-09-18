from sdRDM import DataModel
from EnzymePynetics.modified.protein import Protein
from EnzymePynetics.modified.reactant import Reactant

# Specify EnzymeML version
URL = "https://github.com/EnzymeML/enzymeml-specifications.git"
COMMIT = "72c3d8be4a094983667a7aa62fb599fbc9f7351c"

EnzymeML = DataModel.from_git(URL, COMMIT)
SBOTerm = EnzymeML.enums.SBOTerm
DataTypes = EnzymeML.enums.DataTypes


def parse_enzymeml(
    cls: "Estimator",
    enzymeml: "EnzymeML"
) -> "Estimator":

    if isinstance(enzymeml, str):
        enzymeml, _ = DataModel.parse(enzymeml)

    species = []
    for reactant in enzymeml.reactants:
        species.append(Reactant(**reactant.to_dict()))

    for protein in enzymeml.proteins:
        species.append(Protein(**protein.to_dict()))

    return cls(
        name=enzymeml.name,
        measurements=enzymeml.measurements,
        species=species,
    )
