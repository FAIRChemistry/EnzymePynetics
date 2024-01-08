# EnzymePynetics - Python Tool for Kinetic Model Fitting

[![Generate API](https://github.com/haeussma/EnzymePynetics/actions/workflows/generate_api.yaml/badge.svg)](https://github.com/haeussma/EnzymePynetics/actions/workflows/generate_api.yaml)

## 🛤 What is EnzymePynetics?

EnzymePynetics is a Python-based tool designed for fitting time-course data of enzyme-catalyzed reactions to various kinetic models.

### Key Features of EnzymePynetics:

- **Estimator Initialization**: Utilizes an `Estimator` object for each dataset, derived from an [EnzymeMLDocument](https://github.com/EnzymeML/enzymeml-specifications). This document is created by MTPHandler and includes vital information on measurement data and reaction components.
- **Reaction Equation Definition**: Allows for the specification of educts, products, and enzymes in the reaction equation. Users can define different substrate and enzyme rate laws.
- **Kinetic Parameter Initialization**: Includes parameters like turnover number ($k_{cat}$), Michaelis constant ($K_M$), competitive inhibition constant ($K_{ic}$), uncompetitive inhibition constant ($K_{iu}$), and time-dependent enzyme inactivation rate ($k_{ie}$). These parameters are initialized with specific values and bounds based on the dataset.
- **Parameter Estimation**: Features a robust fitting process using substrate concentration data, integrating the rate laws of a `ReactionSystem`. Utilizes the Lmfit implementation of the Levenberg-Marquardt algorithm for fitting.
- **Comprehensive Analysis and Visualization**: After fitting, it provides an overview of parameters, their standard errors, and the Akaike Information Criterion (AIC) for each system. Additionally, a correlation matrix and interactive visualizations of fitted systems are available for evaluating the fit quality.
- **Integration with EnzymeMLDocument**: Selected kinetic models, along with estimated parameters and uncertainties, are added to the EnzymeMLDocument. This document is then serialized as an SBML-compliant .omex archive for each experimental condition.

## ⚡️ Quick Start

To begin using EnzymePynetics, clone the repository:

```Bash
git clone https://github.com/FAIRChemistry/EnzymePynetics/
