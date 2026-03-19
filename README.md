# ICBResponse

`ICBResponse` is an R package that provides a Shiny application for exploring immunotherapy response data across PAN-T, CAR-T, and CRISPR analyses.

## Overview

The app supports interactive exploration of:

- PAN-T patient-level analysis
- PAN-T cell-level differential expression results
- CAR-T patient-level analysis
- CAR-T cell-level differential expression results
- CRISPR cell-level results

This repository contains the package source code and the Shiny app source.

## Repository structure

Typical contents of this repository include:

- `DESCRIPTION`  
  Package metadata

- `NAMESPACE`  
  Exported functions and imports

- `R/`  
  Package R functions

- `man/`  
  Package documentation

- `inst/app/`  
  Shiny application code

- `inst/app/data/`  
  External data used by the app for local/GitHub workflows

- `tests/`  
  Package tests

## Installation

You can install the package from a local source tarball:

```r
install.packages("ICBResponse_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Data setup

The CRAN package does not include the large Shiny app datasets.

Download the data folder: 
inst/app/data

Unzip it to a local folder, then run:

```r
ICBResponse::run_app(data_dir = "/path/to/data")
```