
# PROJECT STRUCTURE

The main analysis was conducted in the root folder, and the thesis report was done in the `report` folder. Here is a description of the project folder:


1. `data` folder contains all datasets: 

- `cds` contains CDS ERA5 Data. The 2-metre air temperature is in the cds folder, while the 2-metre dewpoint temperature is in its own folder `2m_dewpoint_temperature`. All files are labelled using the year they represent.
- `dhs` contains DHS data (git-ignored) for now (as they are quite heavy).
- `GMST` folder contains data on the Global Mean Surface Temperature used as a covariate in the heatwave exposure model.
- `processed` contains any dataset generated during the analysis. These include the admin-2 and cluster processed heatwave data (we use the cluster processed one in our analysis), the DHS covariates file (`dhs-covariates.rds`) information, the heat index processed information, a master dataset (`master-survey-dataset.rds`) used for building the full sampling scheme, and the outcome information for all antigens (`vaxdata-components.rds`) 
- `shp` contains shapefiles data from GADM, for Nigeria.

2. `output` folder contains any documents, slides/presentations, manuscripts and plot images generated from analysis.

- `docs` contains powerpoint slides for this analysis.
- `draw.io` contains the heatwave exposure construction flowchart.
- `img` contains all images and plots produced in the analysis. images in the root folder are those that do not concern any specific antigen, antigen specific files, are in antigen specific folders, `heatmaps` contain maps and plots of the exposure variables
- `models` contains the rds files for antigen specific models. each file has a heat index model and a heatwave model.

3. `R` contains r scripts for doing analysis:

- `01-allvax-extraction` contains codes for extracting antigens data from DHS.
- `01-cds-data-fetch` contains codes for extracting CDS ERA5 data using their API and R package
- `01-tx5x-analysis` contains codes for analyzing CDS ERA5 data and computing heatwaves
- `01-heatindex-processing` contains code for the construction of the heat index exposure variable.
- `01-dhs-covariates-extraction` contains codes for extracting covariate information from DHS
- `01-merge-N-model` contains codes for merging ERA5 data and DHS data and the model code itself, and all plots concerned
- `autoReg-modifier.R` contains code for modifying functions such as `gaze, autoReg` from the `autoReg` package to work with models from `survey` package such as `svysurvreg`
- `dump-file` for dumping code I might regret deleting later

3. `report` folder contains the MCS Quarto template for writing

- `_book` contains the actual PDF report for the thesis
- `images` is a file created by quarto to host the images in the folder
- `img` contains a few miscellaneous images such as the logos
- `scripts` contains a few R scripts used to produce/reproduce some images, tables in the thesis itself
- `*.qmd` are the actual chapter files for my thesis

4. `resources` folder contains several resources such as DHS guide statistics, Nigeria's Immunization schedule and so on
