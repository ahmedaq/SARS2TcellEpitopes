# Epitopes targeted by T cells in convalescent COVID-19 patients

## Table of Contents
*  [Overview](#overview)
*  [Details](#details)
*  [Requirements](#requirements)
*  [Usage](#usage)
*  [Troubleshooting](#troubleshooting)
*  [Citation](#citation)


## Overview
Summary of emerging data about experimentally-identified SARS-CoV-2 T cell epitopes compiled from 1,259 convalescent COVID-19 patients across 8 independent studies


## Details
#### Title of paper
Epitopes targeted by T cells in convalescent COVID-19 patients
#### Authors
Ahmed A. Quadeer, Syed Faraz Ahmed, and Matthew R. McKay

## Requirements
A PC with MATLAB (preferrably v2019b or later) installed on it with the following additional toolboxes:
    * Bioinformatics Toolbox
    * Statistics and Machine Learning Toolbox


## Usage
1.  Download the repository
2.  Open MATLAB
3.  Change directory to the downloaded repository 
4.  Run the script ProcessingRawData.m. 
    
    Inputs:
    1. "SARS-CoV-2-T-cell-epitopes.xlsx": This table will be periodically updated as new experimental studies reporting SARS-CoV-2 T cell epitopes in COVID-19 convalescent patients are published. 
    2. "AllSARSTepitopes.csv": This table includes details of all SARS-CoV epitopes downloaded from the VIPR database.
    
    Outputs:
    1. "SARS-CoV-2-T-cell-epitopes-RespFreq.xlsx": This table lists all reported unique SARS-CoV-2 epitopes along with their computed response frequencies and 95% confidence intervals.
    2. "SARS-CoV-2-T-cell-epitopes-Mapped-SARS.xlsx": This table lists the SARS-CoV-2 epitopes reported with some HLA allele information. Additional HLA alleles obtained from the genetically-identical SARS-CoV epitope data is also included.
    3. "SARS-CoV-2-T-cell-peptides-Mapped-SARS.xlsx": This table lists the SARS-CoV-2 epitopes reported with no HLA allele information. Precise epitopes and corresponding additional HLA alleles obtained from the genetically-identical SARS-CoV epitope data is also included.


## Troubleshooting
For any questions or comments, please email at eeaaquadeer(at)ust.hk or sfahmed(at)connect.ust.hk.


## Citation
#### Plain text
Ahmed A Quadeer, Syed Faraz Ahmed, Matthew R McKay, Epitopes targeted by T cells in convalescent COVID-19 patients, <i>bioRxiv</i>, 2020.


