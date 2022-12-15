import setuptools
from setuptools import setup

setup(
    name="EnzymePynetics",
    version="0.0.2",
    description="Parameter estimation for enzyme kinetics.",
    url="https://github.com/haeussma/EnzymePynetics",
    author="Haeussler, Max",
    author_email="st171427@stud.uni-stuttgart.de",
    packages=setuptools.find_packages(),
    
    install_requires=[
        "lmfit",
        "matplotlib",
        "numpy",
        "pandas",
        "pydantic",
        "scipy",
        "setuptools"]
    )
