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

### Install the Python package

1. Build the binary: `python setup.py bdist_wheel`
2. Install covigator: `pip install dist/covigator-x.y.z-py3-none-any.whl`

After installation there will be two endpoints available in the path: `covigator-download` and `covigator-pipeline`.

### Setup the database

Given a working Postgres database:

1. Connect to the default database: `psql -W postgres`
2. Create a user for covigator `CREATE USER covigator WITH PASSWORD 'covigator';`
3. Create a database: CREATE DATABASE covigator OWNER covigator;

### Configure access to the database

The application expects some environment variables to be configured with the credentials to the database:
```
COVIGATOR_DB_HOST=0.0.0.0
COVIGATOR_DB_NAME=covigator
COVIGATOR_DB_USER=covigator
COVIGATOR_DB_PASSWORD=covigator
COVIGATOR_DB_PORT=5432
```

If these are not provided it will use the values shown above as default values.