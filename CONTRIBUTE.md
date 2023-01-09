# CoVigator developer guide

## Setup your development environment

 1. Create a virtual environment `virtualenv venv` making sure you are using a Python >= 3.6 interpreter.
 2. Activate your virtual environment `source venv/bin/activate`
 3. Install all dependencies `pip install -r requirements.txt`

## Tests

There are two type of tests:
- **Unit tests**. This need to have two attributes: being fast (in the order of milliseconds) and not depending on any external resource. The objective is to run these constantly during development and in a continuous integration environment and to provide a fast feedback loop. These are under `covigator.tests.unit_tests`
- **Integration tests**. This tests are normally more complex, involve multiple components of the application and/or external resources, and they can be slow. These tests are not intended for automation.  These are under `covigator.tests.integration_tests`

Tests can be run from an IDE like PyCharm or otherwise from the commmand line.

Run all tests as follows:
`python -m unittests discover covigator.tests.unit_tests`

Run a specific test as follows:
`python -m unittest covigator.tests.unit_tests.test_vcf_loader`

**NOTE**: unit tests can make use of the database by initialising it as `Database(test=True)` will start an empty in memory SQLite database.

## Install the Python package

 1. Build the binary: `python setup.py bdist_wheel`
 2. Install covigator: `pip install dist/covigator-x.y.z-py3-none-any.whl`

After installation there will be two endpoints available in the path: `covigator-download` and `covigator-pipeline`.

## Setup the database

Given a working Postgres database:

 1. Connect to the default database: `psql -W postgres`
 2. Create a user for covigator `CREATE USER covigator WITH PASSWORD 'covigator';`
 3. Create a database: CREATE DATABASE covigator OWNER covigator;

## Configure access to the database

The application expects some environment variables to be configured with the credentials to the database:
```
COVIGATOR_DB_HOST=0.0.0.0
COVIGATOR_DB_NAME=covigator
COVIGATOR_DB_USER=covigator
COVIGATOR_DB_PASSWORD=covigator
COVIGATOR_DB_PORT=5432
```

If these are not provided it will use the values shown above as default values.

## Deploy to uWGSI

 1. Install covigator and uwsgi using pip
 2. Configure CoVigator using the necessary environment variables
 3. Go to the folder where the file `covigator/wsgi.py` is located (this can be copied elsewhere) and start the uWSGI server as follows:
```
uwsgi --socket 0.0.0.0:8080 --protocol=http -w wsgi:server
```