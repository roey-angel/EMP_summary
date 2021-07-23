Summarise the taxonomic distribution across different sample types in the Earth Microbiome Project
========

[![Twitter Follow](https://img.shields.io/twitter/follow/espadrine.svg?style=social&label=Follow)](https://twitter.com/RoeyAngel)   ![license](https://img.shields.io/github/license/mashape/apistatus.svg?style=flat-square)


The Earth Microbiome Project ([EMP][1]) is collaborative effort to curate and characterise microbial taxonomy and function on Earth.  
The aim of this script is to summarise the taxonomic distribution across different sample types, as defined by the EMP. The script summarises phyla by default but can also summarise the sequences at different taxonomy levels. The script uses the 3<sup>rd</sup> classification level of the [EMP Ontology][2] (`empo_3`), but this can also be altered.

Overview
--------
    ├── EMP_figures # Plots
    │   ├── plot bar-1.png
    │   ├── plot heatmanp-1.png
    │   ├── plot metacoder-1.png
    │   ├── plot violin-1.png
    ├── EMP_empo_3_Phylum.csv         # Taxonomic distribution by sample type table
    ├── EMP_summary.Rproj # R project setup file
    ├── LICENCE            # Copyright information
    ├── README.md          # Overview of the repo
    ├── Summarise_EMP_taxa.Rmd           # Executable R markdown script
    ├── Summarise_EMP_taxa.html          # HTML output of the RMD file
    └── references.bib          # Bibtex formatted refereces cited in the RMD file


Running the analysis
--------
The RMD file is best executed using [knitr](https://yihui.name/knitr/) on [RStudio](https://www.rstudio.com/). 


[1]: https://earthmicrobiome.org/
[2]: https://earthmicrobiome.org/protocols-and-standards/empo/
