## Requirements

Provided by Martin Löwer on 01/12/2020

1) Download as much sequence data as possible, this will include NGS data (from SRA) and assembled sequences (e.g. from GISAID or NCBI)
2) Process data adequately for each data type, e.g. for NGS data we can have RNA or DNA, amplicon seq ....

3) Data processing should result in a list of mutations, this might also be "somatic" variants relative to the original Wuhan sample, otherwise we these are SNP like variants
(special ephasis on subclonal variants, i.e. intra-host variability)
4) Özlem is trying to retrain her model for using such amplicon seq NGS data for variant detection, having a deep learning model here will be a nice feature

5) variant annotation, in the future we can add NeoFox for epitope/innumogenicity prediction
6) we need to store these results and redo the process every week (or day), analyzing only the new data (of course)
7) after this is running, we need to present the data to the world, can be as simple as a tsv on github in the beginning, but the next step should be a website with fancy user interface
8) and everything should be virus-agnostic, so we can make a fluvigator etc afterwards


## Suggestions

- From Thomas Bukur an example on how to integrate Jupyter notebooks into a web site.
https://covidtracker.fr/covidtracker-france/
Source: https://github.com/rozierguillaume/covid-19


## Design

Layers:
- Download
- Database
- Pipeline
- Presentation

Generally use Python if possible.

(System design)[system_design.png]

### Database

- Postgres
- SQLalchemy as an Object Relational Mapper in Python
- Jobs table (status: DONE | FAILED | DOWNLOADING | PIPELINE_PROCESSING) (NOTE: if processing is lenghty we may need to split this processing step in multiple steps)
- Define other tables to hold the data (ie: variants, samples, etc.)

### Download

- GSAID ???? Manually as no API available. Confirm if this is the case.
- NCBI Virus (NCBI protein), fetch nucleotides instead of aminoacids. Available through Entrez API (XML) in Biopython.
- SRA. SRA fastq-dump from SRA toolkit for download, Entrez API through Biopython to check for the existence of new samples.
- cronjob to schedule download process to get started (probably move to Django eventually)
- Job management using the database

### Installation

- Fetch references using Entrez through Biopython, this may be a scripted but manual task.

### Pipeline

- Use slurm as an HPC queue system
- Clarify approximate execution time per sample
- Prototyping with dask for sending jobs to HPC, clarify if dask allows to determine dependencies between jobs.
- Think on alternatives to dask: queue.py from easyfuse, Nextflow

### API

- TODO

### Presentation

- To get started an export from the database into a TSV will be sufficient
- Plotly for interactive plots
- Bokeh also for interactive plots
- Django?
- Sphinx as a simple alternative?