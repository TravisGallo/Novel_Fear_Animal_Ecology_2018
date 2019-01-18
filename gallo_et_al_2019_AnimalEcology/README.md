# R script, JAGS models, and data used to assess how urbanization changes predator avoidance behaviors in Chicago, IL USA
## **Gallo, T., M. Fidino, E.W. Lehrer, and S. Magle. Urbanization alters predator avoidance behaviors. _Journal of Animal Ecology_ DOI: **


# **File Descriptions:**
**2019-01-15_gallo_et_al_JAE_plotting_results.R** – the only file that you need to open. Sources files below to load data sets, format data for JAGS model, run JAGS model, summarize posterior distributions, calculate derived parameters, and create Figure 2 from manuscript.

**2019-01-15_gallo_et_al_JAE_fit_3sp_JAGSmodel.R** - script that loads data sets and formats data to be fit in JAGS model. This script also calls the JAGS model and summarizes the output. Sourced in 2019-01-15_gallo_et_al_JAE_plotting_results.R.

**2019-01-15_gallo_et_al_JAE_3sp_interaction_JAGSmodel.R** - JAGS model used to estimate the co-occurrence vigilance in white-tailed deer and rabbits. Read in 2019-01-15_gallo_et_al_JAE_fit_3sp_JAGSmodel.R.

**2019-01-15_gallo_et_al_JAE_utility_script.R** - script that loads utility functions. Sourced in 2019-01-15_gallo_et_al_JAE_fit_3sp_JAGSmodel.R.

**2019-01-15_gallo_et_al_JAE_y_matrix_sp10_sp13.txt** - The number of days each species was detected and is supplied as data to the JAGS model so that each species detection probability can be calculated. This array is ordered by species by site by season. Species are in the same order as 2019-01-15_gallo_et_al_JAE_species_used_in_sp10_sp13_analysis.txt and sites are in the same order as 2019-01-15_gallo_et_al_JAE_sites_used_in_sp10_sp13_analysis.txt.

**2019-01-15_gallo_et_al_JAE_z_matrix_sp10_sp13.txt** - Data on whether or not each species were observed at the camera trapping sites each season. If they were detected the cell takes a value of 1, if they were not detected it takes a 0, and if the site was not sampled it takes an NA. This is a species by site by season array in the same orders as 2019-01-15_gallo_et_al_JAE_y_matrix_sp10_sp13.txt.

**2019-01-15_gallo_et_al_JAE_j_matrix_sp10_sp13.txt** - Data on the number of days a camera trap was active at each site and season. Used with the detection data to calculate detection probabilities, and is a site by season matrix. If a site was not sampled a zero is reported.

**2019-01-15_gallo_et_al_JAE_sites_used_in_sp10_sp13_analysis.txt** - A vector of the site names used in our analysis. 

**2019-01-15_gallo_et_al_JAE_species_used_in_sp10_sp13_analysis.txt** - A vector of the observed species whose data are contained in 2019-01-15_gallo_et_al_JAE_y_matrix_sp10_sp13.txt and 2019-01-15_gallo_et_al_JAE_z_matrix_sp10_sp13.txt.

**2019-01-15_gallo_et_al_JAE_urban_covs.csv** - Housing density, tree cover, and impervious cover data at each sampling site used to calculate the urbanization index used in our models.

**2019-01-15_gallo_et_al_JAE_activity.data.csv** - Data for each photograph of coyote, white-tailed deer, and Eastern cottontail used to estimate daily activity patterns.

**2019-01-15_gallo_et_al_JAE_deer_vigilance.csv** - Data indicating if a deer had their head up or down in a photograph.

**2019-01-15_gallo_et_al_JAE_rabbit_vigilance.csv** - Data indicating if a rabbit had their head up or down in a photograph.

**Note:** All of these files must be within your working directory for the analysis to work. Our analysis was done in parallel and used all but two of the cores in a computer. Therefore, if you have two or less cores on your computer you will need to adjust the settings annotated within 2019-01-15_gallo_et_al_JAE_fit_3sp_JAGSmodel.R.

![Eastern cottontail line drawing by mason](https://github.com/TravisGallo/Novel_Fear_Animal_Ecology_2018/blob/master/gallo_et_al_2019_AnimalEcology/coyote-deer-rabbit.png)
