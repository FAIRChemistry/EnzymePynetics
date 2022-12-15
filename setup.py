import setuptools
from setuptools import setup

setup(
    name="EnzymePynetics",
    version="0.0.1",
    description="Parameter estimation for enzyme kinetics.",
    url="https://github.com/haeussma/EnzymePynetics",
    author="Haeussler, Max",
    author_email="st171427@stud.uni-stuttgart.de",
    packages=setuptools.find_packages(),
    
    install_requires=[
        "ipython==8.6.0",
        "lmfit==1.0.3",
        "matplotlib==3.6.2",
        "numpy==1.22.3",
        "pandas==1.4.2",
        "pydantic==1.9.0",
        "scipy==1.8.1",
        "setuptools==61.2.0"]
    )
