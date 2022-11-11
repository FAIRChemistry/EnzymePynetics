from kineticmodel import KineticModel, auto_inhibition_model_dict

class Fitter:
    def __init__(self, name: str):
        self.name = name

if __name__ == "__main__":
    test = Fitter("max")
    print(test)