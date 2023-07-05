![CoVigator logo](covigator/dashboard/assets/CoVigator_logo_txt_reg_no_bg.png "CoVigator logo")

-----------------

# CoVigator: monitoring SARS-CoV-2 mutations

[![PyPI version](https://badge.fury.io/py/covigator.svg)](https://badge.fury.io/py/covigator)
[![Run unit tests](https://github.com/TRON-Bioinformatics/covigator/actions/workflows/unit_tests.yml/badge.svg?branch=main)](https://github.com/TRON-Bioinformatics/covigator/actions/workflows/unit_tests.yml)
[![codecov](https://codecov.io/gh/TRON-Bioinformatics/covigator/branch/main/graph/badge.svg?token=J5Q8UV65PD)](https://codecov.io/gh/TRON-Bioinformatics/covigator)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/d6735b902b7b42e0a7cd423cebff69d2)](https://www.codacy.com/gh/TRON-Bioinformatics/covigator/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=TRON-Bioinformatics/covigator&amp;utm_campaign=Badge_Grade)
[![Powered by Dash](https://img.shields.io/badge/powered%20by-Dash-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://dash.plotly.com/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![Documentation Status](https://readthedocs.org/projects/covigator/badge/?version=latest)](https://covigator.readthedocs.io/en/latest/?badge=latest)

![CoVigator art](covigator/dashboard/assets/wordcloud.png)

**CoVigator dashboard**: [https://covigator.tron-mainz.de](https://covigator.tron-mainz.de)

**CoVigator documentation**: [https://covigator.readthedocs.io/](https://covigator.readthedocs.io)

Human infections with SARS-CoV-2 are spreading globally since the beginning of 2020, necessitating preventive or 
therapeutic strategies and first steps towards an end to this pandemic were done with the approval of the first mRNA 
vaccines against SARS-CoV-2. 
The accumulation of virus samples that have been sequenced in a short time frame is unprecedented (see Figure 1).
This is the first pandemic recorded at a molecular level with such level of detail giving us the opportunity to develop
new tools for the monitoring of its evolution.

We want to provide an up-to-date interactive view on SARS-CoV-2 mutations to support global efforts in preventing or 
treating infections. 
Monitoring the appearance of relevant new mutations is key to enable a fast reaction to new strains and for that 
purpose we enable the exploration of these mutations and their annotations (see Figure 2). 
Thus, we envision to help guiding global vaccine design efforts to overcome the threats of this pandemic.

CoVigator is a monitoring system for SARS-CoV-2 which integrates a full variant calling pipeline, 
a database that stores all relevant information about mutations in SARS-CoV-2 and finally a dashboard to enable 
visual analytics.

![CoVigator sample accumulation](docs/source/_static/figures/screencast_01_samples_by_country_tab.gif)

<p align = "center">
<b>Figure 1: Sample accumulation by country and dN/dS evolution</b>
</p>

![CoVigator gene S view](docs/source/_static/figures/screencast_03_recurrent_mutations_tab.gif)

<p align = "center">
<b>Figure 2: Recurrent mutations</b>
</p>

CoVigator loads publicly available SARS-CoV-2 DNA sequences from the 
[European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena) providing raw reads in FASTQ format.

ENA enables a high resolution into the SARS-CoV-2 mutation details through the individual 
reads. This allows us to annotate mutations with a Variant Allele Frequency (VAF) and explore intrahost 
mutations. Importantly, we only process the Illumina samples from ENA. 
This means excluding all of the Oxford Nanopore samples and hence having a partial view of all the available data.

The dashboard is implemented in the visualization framework [Dash](https://dash.plotly.com/). 
The computation is distributed through our cluster with a library of similar name [Dask](https://dask.org/).
The analysis pipeline is implemented in the [Nextflow](https://www.nextflow.io/) framework.

![CoVigator system](docs/source/_static/figures/system_design_manuscript.png)

<p align = "center">
<b>Figure 3: System design</b>
</p>

The CoVigator project was developed at the Biomarker Development Center at 
[TRON (Translational Oncology at the University Medical Center of the Johannes Gutenberg University gGmbH)](https://tron-mainz.de/). 
The project was kindly supported by 
[Intel´s Pandemic Response Technology Initiative](https://newsroom.intel.com/tag/pandemic-response-technology-initiative).

## How to cite

*   Bukur T., Riesgo-Ferreiro P., Sorn P, Gudimella R., Hausmann J., Rösler T., Löwer M., Schrörs B., & Sahin U. 
CoVigator — A Knowledge Base for Navigating SARS-CoV-2 Genomic Variants. Viruses. 2023; 15(6):1391. [10.3390/v15061391](https://doi.org/10.3390/v15061391)

*   Schrörs, B., Riesgo-Ferreiro, P., Sorn, P., Gudimella, R., Bukur, T., Rösler, T., Löwer, M., & Sahin, U. (2021). 
Large-scale analysis of SARS-CoV-2 spike-glycoprotein mutants demonstrates the need for continuous screening of virus 
isolates. PLOS ONE, 16(9), e0249254. [10.1371/journal.pone.0249254](https://doi.org/10.1371/journal.pone.0249254)

### Acknowledgements

The lineage annotation in the covigator dashboard package uses "the description of constellations of mutations for the SARS-CoV-2 virus" by [cov-lineages](https://github.com/cov-lineages/constellations) provided and licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)
