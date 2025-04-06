# FEED_farmer_disease_management_heterogeneity

This repository contains files for the analysis presented in the scientific paper "Incorporating heterogeneity in farmer disease control behaviour into a livestock disease transmission model" by Edward M. Hill, Naomi S. Prosser, Paul Brown, Eamonn Ferguson, Jasmeet Kaler, Martin J. Green, Matt J. Keeling and Michael J. Tildesley.

Zenodo DOI for the repository:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14229548.svg)](https://doi.org/10.5281/zenodo.14229548)

Publication details: EM Hill, NS Prosser, P Brown, E Ferguson, MJ Green, J Kaler, MJ Keeling, MJ Tildesley. (2023) Incorporating heterogeneity in farmer disease control behaviour into a livestock disease transmission model. *Preventive Veterinary Medicine*. 219: 106019. doi: 10.1016/j.prevetmed.2023.106019. URL: https://doi.org/10.1016/j.prevetmed.2023.106019.

## Corrections to the publication 

### 27 November 2024

Within the livestock disease model code there was a bug identified with the cell-to-cell distance calculation.
Rerunning the analysis with the amended cell-to-cell distance calculation, there was no change in the qualitative findings. There were marginal changes to the quantitative results. Updated versions of the main manuscript and supporting information are available. 
 - **Main manuscript**: https://github.com/EdMHill/FEED_farmer_disease_management_heterogeneity/blob/main/docs/2024-11-27-manuscript_corrections/2024-11-27-main_manuscript_corrections.pdf
 - **Supporting information**: https://github.com/EdMHill/FEED_farmer_disease_management_heterogeneity/blob/main/docs/2024-11-27-manuscript_corrections/2024-11-27-supporting_information_corrections.pdf

 A correction has also been published by *Preventive Veterinary Medicine*: https://doi.org/10.1016/j.prevetmed.2024.106408.

 ### 06 April 2025
Within the livestock disease model code there there was an implementation error in the conditional subsampling algorithm used to model the spatial transmission of infection. Rerunning the analysis with the amended code, our previous qualitative conclusions are unchanged. There were changes to the quantitative results, mainly a reduction in the outbreak sizes. Updated versions of the main manuscript and supporting information are available. 
 - **Main manuscript**: https://github.com/EdMHill/FEED_farmer_disease_management_heterogeneity/blob/main/docs/2025-04-06-manuscript_corrections/2025-04-06-main_manuscript_corrections.pdf
 - **Supporting information**: https://github.com/EdMHill/FEED_farmer_disease_management_heterogeneity/blob/main/docs/2025-04-06-manuscript_corrections/2025-04-06-supporting_information_corrections.pdf

## Livestock disease model 

The livestock disease model simulations are performed using the programming language Julia.
Julia makes use of environments, allowing bespoke package lists for separate projects. Documentation on working with environments and installing packages in the same state that is given by the project manifest: https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project

The livestock disease model analyses in the paper require holding location and demographic cattle data that are not publicly available. We instead provide synthetic cattle demography data for Great Britain, with the number of holdings matching the number used in our analyses, albeit randomly located in a 1000kmx250km region. This enables the user to run the main simulation script (run_GB_livestock_disease_spatial_model.jl). Note the results produced using the synthetic data do not match those from the paper.

We provide results files that allow the user to run the associated results scripts.

## Repository structure

Please find below an explainer of the directory structure within this repository.

### data

 - **datavis**  
Visualisations of interview data

 - **dummy_data**  
Directory containing synthetic data on livestock counts and premises location data.

 - **elicitation_interviews**  
The file *cluster_assignment_five_variable_model_with_herd_size.csv* has the herd sizes and behavioural cluster assignment for the five variable statistical model.

The file *elicitation_interview_data.csv* has the farmer-level behaviour and psychosocial score data from the elicitation interviews that was analysed for the associated paper.

 - **shapefiles**  
County and country boundary shapefiles for UK. Source: Office for National Statistics licensed under the Open Government Licence v.3.0. Contains OS data © Crown copyright and database right 2022. Counties and Unitary Authorities (December 2020) UK BGC dataset: https://www.arcgis.com/home/item.html?id=aff50e8d15364a7b82c62c14861eb240

### results

The file *elicitation_interview_analysis_script.R* contains the R code for analysis of the farmer elicitation interview data for the results presented in the paper and fed into the livestock disease model. 

#### livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs
Directories to store simulation outputs from the livestock disease model.

#### livestock_disease_model/generate_figures
The main plotting script is generate_plots_script.jl, which calls functions in generate_plots_supporting_fns.jl. 
The script compute_intervention_cost_stats.jl generates the threshold intervention costs.

 - **bar_plots**    
 Directory for bar plots

 - **JLD2_files**   
 Houses the cost related statistics file. Generated by compute_intervention_cost_stats.jl

 - **seed_region_spatial_maps**  
 Directory for spatial maps

 - **violin_plots**  
 Directory for violin plots

#### livestock_disease_model/table_summ_stats
The script generate_table_summ_stats_script.jl computes the statistics for the results tables associated with the analysis carried out using the livestock disease model.

### src

#### livestock_disease_model

 - **common_fns**  
Files that house commonly used functions across multiple spatial simulation projects

 - **GB_model_with_elicitation_data**  
Houses directories with mathematical model simulation files. Main simulation file is run_GB_livestock_disease_spatial_model.jl

    - **config_JLD2_files**     
    Simulation parameter information save file in JLD2 format

    - **config_log_files**      
    Simulation parameter information save file in txt format

    - **supporting_fns**    
    Collection of files that contain the functions called in the main outbreak simulation.
