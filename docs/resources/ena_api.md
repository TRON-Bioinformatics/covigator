# European Nucleotide Archive (ENA) API

The ENA API is available here https://www.ebi.ac.uk/ena/portal/api/
Documentation for the API is available here https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html

## Search 

The most relevant endpoint from the API is the search endpoint that can be used to retrieve all results following certain criteria.

'''
/search?result=read_run&query=tax_eq({tax_id})&limit={page_size}&offset={offset}&fields={fields}&format=json
'''

where:
- `{tax_id}`: the taxomic identifier of the organism of interest. Sars-CoV-2 has the id `2697049`. See Fetching taxomic identifiers
- `{page_size}`: the number of results to be returned
- `{offset}`: the offset in the results, ie: this is used to paginate over the results
- `{fields}`: this is a comma separated lists of metadata fields to be returned, eg: `scientific_name,study_accession,experiment_accession,first_created`. The available metadata fields are described in the API too here https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run&format=json

Example fetching 3 entries with an offset of 10:
https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq(2697049)&limit=3&offset=10&fields=scientific_name,study_accession,experiment_accession,first_created,collection_date,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,sample_collection,sequencing_method,center_name,fastq_ftp,fastq_md5,host_tax_id,host_sex,host_body_site,host_gravidity,host_phenotype,host_genotype,lat,lon,country&format=json

### Metadata

scientific_name,study_accession,experiment_accession,first_created,collection_date,instrument_platform,instrument_model,sample_collection,sequencing_method,center_name,fastq_ftp,fastq_md5,
# data on host
host_tax_id,
host_sex,
host_body_site,
host_gravidity,
host_phenotype,
host_genotype,
# geographical data
lat,
lon,
country


## Fetching taxonomic identifiers

Taxonomic identifiers can be fetched searching by the scientific name of the organism of interest

https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/Severe%20acute%20respiratory%20syndrome%20coronavirus%202

or

https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/Homo%20sapiens