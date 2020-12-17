How fast do biological rock crusts grow?
========

[![Twitter Follow](https://img.shields.io/twitter/follow/espadrine.svg?style=social&label=Follow)](https://twitter.com/RoeyAngel)   ![license](https://img.shields.io/github/license/mashape/apistatus.svg?style=flat-square)


Data, sequence analysis ploting for the microbiome figures included in the paper: [Nimrod Wieler, Tali Erickson Gini, Roey Angel and Osnat Gillor: Estimating bacterial rock crust growth rate using a well-dated millennia-old arid archaeological site]


Overview
--------
    ├── Data          # Primary data
    │   ├── Shivta_site_silva.nrv119.taxonomy
    │   ├── Shivta_site_otuTab2.txt
    │   ├── Shivta_metadata.csv
    │   └── Shivta_site_otuReps.filtered.align.treefile
    ├── LICENCE            # Copyright information
    ├── README.md          # Overview of the repo
    ├── references.bib          # Bibtex formatted refereces cited in the RMD file
    ├── Results         # Result files generated while running the RMD file
    ├── BRC_growth_rate_microbiome_analysis.Rproj          # R project file
    ├── BRC_growth_rate_microbiome.html          # HTML output of the RMD file
    ├── BRC_growth_rate_microbiome.md          # MD output of the RMD file
    ├── BRC_growth_rate_microbiome.RMD            # Executable R markdown script
    └── BRC_growth_rate_generate_OTUS.sh          # Shell excecutable for processing the sequence data

Running the analysis
--------
The RMD file is best executed using [knitr](https://yihui.name/knitr/) on [RStudio](https://www.rstudio.com/). 
