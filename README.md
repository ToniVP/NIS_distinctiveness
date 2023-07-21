# NIS_distinctiveness
This repository contains all the compiled data, code and further generated data from the work "**A trait-based approach to assess niche overlap and functional distinctiveness between non-indigenous and native species**" authored by Antoni Vivó-Pons, Mats Blomqvist, Anna Törnroos and Martin Lindegren.

## Data description
The folder **Raw_data** contains the original data files from which all the other databases included in the **Data** folder are derived. 
#### Raw_data
- _env_data.txt_ contains all the environmental data extracted from the ice-ocean NEMO nordic model.
- _sp_status.txt_ accounts for the name, phylum and origin (native or NIS) from all the included taxa.
- _sp_traits_raw.txt_ trait information for all species.
- _species_AFDW_2005_2020.txt_ species abundance measures (n individuals, wet weight, AFDW) and details for all the specific sampling events
- _species_site_AFDW_2005.txt_ site x species biomass (AFDW) matrix

#### Data


## Code description

[01_Descriptive_figures](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/01_Descriptive%20figures.R). This script contains the code related to the **descriptive figures of the data** used in the study.

[02_Distinctiveness_regional (Step I & II)](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/02_Distinctiveness_regional%20(Step%20I%20%26%20II).R). This script corresponds to the **1st** and **2nd steps** from the proposed framework, related to the computation of functional distinctiveness of NIS and natives at a **regional scale** (regional species pool) allowing us to obtain several outputs needed in further steps.

[03_Distinctiveness_local (Step III)](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/03_Distinctiveness_local%20(Step%20III).R). This script corresponds to the **3rd step** from the proposed framework, related to the computation of functional distinctiveness of NIS together with other community metrics at a **local scale** (local species pools). 

[04_Statistical_analysis (Step III)](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/04_Statistical_analysis%20(Step%20III).R). This script contains all the code related with the **multi- and single modeling approaches** included in the study to **detect potential drivers** of NIS distinctiveness at a **local scale**. 

[Env_data_extraction](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/Env_data_extraction.R). Example of extracting environmental data from several .NC files, downloaded from the [Copernicus Marine Service](https://data.marine.copernicus.eu/products).





