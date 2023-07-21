# NIS_distinctiveness
This repository contains all the compiled data, code and further generated data from the work "**A trait-based approach to assess niche overlap and functional distinctiveness between non-indigenous and native species**" authored by Antoni Vivó-Pons, Mats Blomqvist, Anna Törnroos and Martin Lindegren.

**_IMPORTANT_**: ALL the plots and figure panels shown in the manuscript were post-processed in Adobe Illustrator, therefore the "raw" versions of the plots obtained with the following code might have different variable names or not be gathered in the same way as they appear in the manuscript. Additionally, the NIS illustrations shown within some figures were also created using Adobe Illustrator. 
In any case, the provided code allows for obtaining identical, but raw versions of all the figures included in the manuscript. 

## Data description
The folder **Raw_data** contains the original data files from which all the other databases included in the **Data** folder are derived.

#### Raw_data
- _env_data.txt_; contains all the environmental data extracted from the ice-ocean NEMO nordic model.
- _sp_status.txt_; accounts for the name, phylum and origin (native or NIS) from all the included taxa.
- _sp_traits_raw.txt_; trait information for all species.
- _species_AFDW_2005_2020.txt_; species abundance measures (n individuals, wet weight, AFDW) and details for all the specific sampling events
- _species_site_AFDW_2005.txt_; site x species biomass (AFDW) matrix

#### Data
- _Di_metrics_station.txt_; several metrics obtained at local scale, in each sampling event, including NIS local distinctiveness, species richness and Shannon index, among others.
- _Trait_modalities_; list with all the trait and corresponding modalities.
- _dist_matrix_ovrll.txt_; overall pairwise matrix of functional distances between species.
- _sp_traits.txt_; cleaned species traits used in the analysis.
- _spe_index.txt_; distinctiveness values for all species in the regional pool.
- _traits.effects.txt_; effect of each trait on regional distinctiveness for each species. 

## Code description
[Libraries](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/Libraries.R). List of required libraries.

[Functions](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/Functions.R). List of needed functions.

[01_Descriptive_figures](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/01_Descriptive%20figures.R). This script contains the code related to the **descriptive figures of the data** used in the study.

[02_Distinctiveness_regional (Step I & II)](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/02_Distinctiveness_regional%20(Step%20I%20%26%20II).R). This script corresponds to the **1st** and **2nd steps** from the proposed framework, related to the computation of functional distinctiveness of NIS and natives at a **regional scale** (regional species pool).

[03_Distinctiveness_local (Step III)](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/03_Distinctiveness_local%20(Step%20III).R). This script corresponds to the **3rd step** from the proposed framework, related to the computation of functional distinctiveness of NIS together with other community metrics at a **local scale** (local species pools). 

[04_Statistical_analysis (Step III)](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/04_Statistical_analysis%20(Step%20III).R). This script contains all the code related with the **multi- and single modeling approaches** included in the study to **detect potential drivers** of NIS distinctiveness at a **local scale**. 

[Env_data_extraction](https://github.com/ToniVP/NIS_distinctiveness/blob/main/Code/Env_data_extraction.R). Example of extracting environmental data from several .NC files, downloaded from the [Copernicus Marine Service](https://data.marine.copernicus.eu/products).





