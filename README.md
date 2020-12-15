# CoVigator (Corona Navigator)

CoVigator is a fully automatized Corona Analysis Pipeline integrating the following major steps:

* Generation/Download of Reference DBs
* Starting of different workflows depending of input data type (e.g. RNA-Seq, Assembly)

# Configuration

`config.py` contains all the essential configuration paths

# Usage

```
python processing.py -i <path_to_input_fastqs> -w <working_dir>
```

whereas <input_paths> should be a list of FASTQs or a path where the FASTQs are stored and
<working_dir> describes the path to store the results in.

## Developer guide

### Setup your development environment

1. Create a virtual environment `virtualenv venv` making sure you are using a Python >= 3.6 interpreter.
2. Activate your virtual environment `source venv/bin/activate`
3. Install all dependencies `pip install -r requirements.txt`

### Run unit tests

Unit tests can be run from an IDE like PyCharm or otherwise from the commmand line `python -m unittests discover covigator.tests`.

### Install the application

1. Build the binary: `python setup.py bdist_wheel`
2. Install covigator: `pip install dist/covigator-x.y.z-py3-none-any.whl`

After installation there will be two endpoints available in the path: `covigator-download` and `covigator-pipeline`.