<h2 align="center">
  Markdown Example for sdRDM
</h2>

<p align="center"> 
This is an example of how to set up a data model using the Software-Driven Research Data Management (sdRDM) library which is based on abstract object models. Furthermore, the sdRDM library supports the conversion of data models defined in the Markdown format.</p>

### ✍️ Syntax

Data models defined in the Markdown format follow these conventions:

- **Modules** are denoted by a heading level 1 ```#```
- **Objects** are started with a heading level 3 ```###``` 
- Each object contains **fields** in bold as a list &rarr; ```- __name__```
- **Required fields** are denoted with an asterix &rarr; ```- __name*__```
- Each field has **options** as a list of name to value mapping &rarr; ```- Type: string```

### ⚙️ Field options

Each field in an object can hold options relevant for mapping to another data model (e.g. a standardized format) and general information such as its type and description. In the following is a collection of all native and required fields:

- **Type** - Required option to denote the data type. Please note, this can also contain other objects defined in this document.
- **Multiple** - Whether or not this field can contain multiple values. Setting to ```True```will result in a ```List[dtype]``` annoatation in the software.
- **Description** - Required option to describe the field. This should be a brief description that explains what the attribute is about.

### 🧬 Inheritance

In order to inherit attributes to another object, the object definition additionally includes the name of the parent object in italic wrapped with brackets:

&rarr; ```## Child [_Parent_]```

In the following an [example](https://github.com/JR-1991/sdrdm-template/tree/main/specifications) data model is defined using above rules. Feel free to use this example also as a template for your own application.

### 👁 How can I use it by myself?

You can experiment and use this [example](https://github.com/JR-1991/sdrdm-template/tree/main/specifications) repository right away to get familiar with teh concept. This repository includes an [action](https://github.com/JR-1991/sdrdm-template/blob/main/.github/workflows/generate_api.yaml) that is triggered whenever changes are pushed. Thus, when you introduce changes to the markdown document, these will directly be reflected onto the generated software. Follow these steps to start out:

1. Fork this repository into your own profile. This will create an exact copy, but you have all rights to modify it without affecting the original.

![](https://www.earthdatascience.org/images/earth-analytics/git-version-control/githubguides-bootcamp-fork.png)

2. Open the ```specifications/Example.md``` file and edit it according to the syntax. You can also press ```Preview``` to inspect the rendered Markdown.
   
![](https://docs.github.com/assets/cb-118903/images/help/repository/edit-file-edit-dropdown.png)

3. Commit changes to the ```main``` branch or create a new one from it. By creating a new branch you can safely work without affecting the original. Once your modifications are done, you can merge these into the ```main``` branch.

![](https://docs.github.com/assets/cb-32137/images/help/repository/choose-commit-branch.png)

4. Watch your changes being reflected onto the API. You can also directly fetch this model using the [sdRDM](https://github.com/JR-1991/software-driven-rdm) library. For this, you can use the following example code that should run as is. 

```python
from sdRDM import DataModel

Root, Author, Parameter = DataModel.from_git(
    url="https://github.com/JR-1991/sdrdm-template.git",
    import_modules = ["Author", "Parameter"]
)

# Visualize the data model
Root.visualize_tree()

# Enter your data
dataset = Root(title="Some Title", description="Some Description")
dataset.add_to_authors(name="Jan Range", affiliation="SimTech")
dataset.add_to_parameters(key="Param", value=10.0)

# Inspect your dataset
print(dataset.yaml())

# Option: Link your dataset to an option --> Dataverse
dataset.to_dataverse()

# Option: Export your dataset to another format
with open("my_dataset.json", "w") as f:
    f.write(dataset.json())

# Re-opening your dataset using sdRDM will cause the library
# to re-build the software state in which the dataset was created

```

*(Images were taken from GitHub's ["Editing Files" tutorial](https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files))*
