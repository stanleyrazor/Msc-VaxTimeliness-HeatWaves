
# PROJECT STRUCTURE

1. `data` folder contains all datasets: 

- cds contains CDS ERA5 Data
- dhs contains DHS data (git-ignored) for now, 
- processed contains any dataset generated during the analysis and 
- shp contains shapefiles data.

2. `output` folder contains any documents, slides/presentations, manuscripts and plot images generated from analysis
3. `R` contains r scripts for doing analysis:

- `01-allvax-extraction` contains codes for extracting antigens data from DHS.
- `01-cds-data-fetch` contains codes for extracting CDS ERA5 data using their API and R package
- `01-tx5x-analysis` contains codes for analyzing CDS ERA5 data and computing heatwaves
- `o1-dhs-covariates-extraction` contains codes for extracting covariate information from DHS
- `01-merge-N-model` contains codes for merging ERA5 data and DHS data and the model code itself
- `dump-file` for dumping code I might regret deleting later

3. `report` filder contains the MCS Quarto template for writing

4. `resources` folder contains several resources such as DHS guide statistics, Nigeria's Immunization schedule and so on
