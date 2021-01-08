# CoVigator (Corona Navigator)

CoVigator is a fully automatized SARS-CoV-2 analysis pipeline integrating the following major steps:

* Generation/Download of Reference DBs
* Starting of different workflows depending of input data type (e.g. RNA-Seq, Assembly)

The system architecture has the following components:
- Accessor: queries external systems for new data and creates an entry in the local database with all necessary metadata
- Processor: reads from the database samples to be processed and triggers a workflow with multiple steps
    - Download: downloads the necessary FASTQ files
    - Pipeline: triggers a variant calling pipeline which outputs VCF files
    - Delete: deletes the FASTQ files
    - Load: loads the variants from VCF files into the database
- Application Programming Interface (API): entry point to fetch data from the database
- Front end: a web application reading data from the API and presenting the data through tables and visualizations

While the accessor and the processor are backend processes that are intended to run asynchronously and periodically, the API and the front end are accessible by end users.

Although the initial use case for Covigator is SARS-CoV-2 data, it is intended to be usable with other infectious organisms.

![Covigator system design](docs/resources/system_design.png "Covigator system design")

## Accessor

The accessor queries external systems, checks in the database which samples are new and creates the required entries in the database.
The accessor is intended to run periodically.
The accessor currently uses the European Nucleotide Archive (ENA) to fetch NGS raw data and it is intended to use GSAID to fetch assemblies in the future.
The accessor does not download large raw files, but only the sample metadata, the URLs and MD5 check sums required to download the data.
The accessor implements a retry mechanism with an exponential backoff to manage temporary failures of 
service of any external data provider.

### Input data

- The organism taxonomic identifier (eg: for SARS-CoV-2 the taxonomic identifier is 2697049)
- The host organism taxonomic identifier (eg: for Homo sapiens the taxonomic identifier is 9606)

The taxonomic identifiers for the different organisms is available through EMBL-EBI as described here https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/taxon-api.html or through NCBI here https://www.ncbi.nlm.nih.gov/taxonomy.

### Configuration

The only configuration required by the accessor is the database configuration. The configuration is done through environment variables.

- `COVIGATOR_DB_HOST`: the host of the database (default value 0.0.0.0)
- `COVIGATOR_DB_NAME`: the database name inside Postgres (default value: covigator)
- `COVIGATOR_DB_USER`: the database user (default value: covigator)
- `COVIGATOR_DB_PASSWORD`: the database password (default value: covigator)
- `COVIGATOR_DB_PORT`: the database port (default value: 5432)

### Usage

`covigator-ena-accessor --tax-id 2697049 --host-tax-id 9606`

## Processor

The processor is in charge of orchestrating the load of samples into the database. 
The processor is intended to run periodically.
Whenever the accessor finds a new sample not present in the database it stores all the required sample metadata and it creates a new job. 
This job is the starting point for the processor which orchestrates the flow of the job through its life cycle.

The happy path of a job is the following:
- `PENDING`: newly created job by the accessor
- `QUEUED`: the job has already been read by the processor and the subsequent actions are scheduled
- `DOWNLOADED`: intermediate state between the downloader and the pipeline
- `PROCESSED`: intermediate state between the pipeline and the clean up
- `LOADED`: final state once the pipeline results have been loaded into the database
  
The failure states are the following and are self descriptive
- `FAILED_DOWNLOAD`
- `FAILED_PROCESSING`
- `FAILED_LOAD`

There is a cleanup task that happens after processing in parallel to the load process, 
but this is not a state of the job and failures in the cleanup although logged in the database do not stop a job from 
reaching the final state `LOADED`.

A timestamp is stored for every change of state. The failure states also stores an error message.

![Covigator system design](docs/resources/processor_state_diagram.png "Covigator system design")

The above workflow is orchestrated using dask library for parallelization. 
The different processes have different priorities to ensure that the graph of tasks is processed in a depth first order.
This is relevant when processing a large amount of jobs to avoid that all jobs are first downloaded, thus requiring that all raw files are stored in the file system simultaneously.
The tasks managed in this workflow are not computationally intensive and each uses a single CPU and a low amount of memory.
The computationally intensive tasks are send by the pipeline to a cluster.

### Input data

- The number of available CPUs for the processor. This corresponds to the number of samples that will be processed simultaneously.

### Configuration

The processor requires the database configuration. See above section in the accessor configuration for details. 

Additionally, the downloader requires the path where downloaded files will be stored.
The configuration is done through environment variables.

- `COVIGATOR_STORAGE_FOLDER`: the folder where files will be stored by the downloaded (default value `./data/covigator`)

### Usage

`covigator-processor --num-cpus 32`

### Processes within the workflow

#### Downloader

The downloader takes an ENA run, downloads all of its FASTQs, 
stores them in a folder structures such as `$COVIGATOR_STORAGE_FOLDER/${run_accession}` 
(eg: `/covigator/data/ERR12345/ERR12345.fastq.gz`) and then checks the MD5 checksum of the 
downloaded files.

The downloader streams large files into disk without holding them into memory, thus enabling the download of large files.
The downloader implements a retry mechanism with an exponential backoff to manage temporary failures of 
service of any external data provider.

#### Pipeline

*Configuration*

`config.py` contains all the essential configuration paths

*Usage*

```
python processing.py -i <path_to_input_fastqs> -w <working_dir>
```

whereas <input_paths> should be a list of FASTQs or a path where the FASTQs are stored and
<working_dir> describes the path to store the results in.

#### Clean up

Deletes the FASTQ files and eventually intermediate files left behind by the pipeline.

#### Loader

**TODO**: teads the outcome of the pipeline (in VCF files?) and loads it into the database.



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