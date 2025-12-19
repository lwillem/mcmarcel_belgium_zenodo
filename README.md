# README: RSV modelling project McMarcel (BEL)

`The objective of this RSV modelling project is to evaluate the impact and cost-effectiveness of potential maternal and neonatal RSV immunisation strategies in Belgium. The modelling structure is a static model that tracks monthly birth cohorts of children from birth to 5 years old. The acronym stands for **M**ulti-**C**ountry **M**odel **A**pplication for **R**SV **C**ost-**E**ffectiveness po**L**icy (McMarcel).

McMarcel has been conceived to study the application of maternal vaccination and monoclonal antibody prevention strategies for 73 Gavi-eligible countries. Later on, we adopted the code to perform studies (West-)European countries or regions. We are aware that our approach underestimates the benefits of vaccination (without additional adjustments to cope with potential herd immunity impacts), though we are keeping this strategy for feasibility reasons.

This software contains the model implementation and input to conduct a health economic analysis and is distributed under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 licence. If you use (parts of) this modelling project, please cite: **Getaneh et al. (2023) Cost-effectiveness of monoclonal antibody and maternal immunization against respiratory syncytial virus (RSV) in infants: Evaluation for six European countries. Vaccine. 41(9), 1623-1631.**

To facilitate public access to the model methods while preserving data confidentiality, we applied stochastic data perturbation to key epidemiological input data. Specifically, we added non-systematic noise to these files, which are marked with an "X" (e.g. BELX or RESCEUX). This approach enables the reproduction of output that is epidemiologically non-actionable while retaining internal consistency for methodological replication.

Note that the original framework includes automated output testing and parallel processing. To enhance transferability, we provide a serial implementation here. Please contact us for additional model features and post-processing code.`


## Where to start?

We organised our R-project upon **RSV_main.R**, which is the workbench to coordinate the model input, the number of runs, random seeds,... and all output. Please make sure that your [**Working Directory**](https://stat.ethz.ch/R-manual/R-devel/library/base/html/getwd.html) is specified as the location of this R file, since all links in de code are relative to the location of this script.

## Where are my results stored?

During the execution of the model (by using the RSV_main script), an **output** directory will be created. Each run will create a separate sub-directory to store all results.

## How to specify intervention options and countries?

The main script uses a csv file in the "config" directory, of which each row specifies one intervention program to be included in the cost-utility analysis for a specific country and scenario. Simulation results are grouped per "config_tag" to create figures and tables.

## Main directory

<p align="left">

| Directory or file name | Content                                                                                                                                                 |
|------------------|------------------------------------------------------|
| config                 | Directory with model configurations by country, scenario and mAb and maternal intervention characteristics (coverage, duration of protection, efficacy) |
| functions              | Directory with R help functions                                                                                                                         |
| input                  | Directory with input data to run the analysis                                                                                                           |
| README.pdf             | Readme file with the project introduction and reference                                                                                                 |
| RSV_main.R             | Main script to run the cost-effectiveness analysis                                                                                                      |
| LICENCE.txt            | GNU GENERAL PUBLIC LICENSE Version 3                                                                                                                    |

</p>

All code is tested with R Version 4.3.1 on MacOS 13.7


## Disclaimers

-   Discounting of future cost and effects is implemented with the assumption that the modelling time horizon aligns with the ageing cohort. This is OK when the intervention takes place at birth with an year-round or seasonal program. With a catch-up program, some cohorts will receive their immunization at e.g. age 6 months. As long as the catch-up program is assumed to prevent only burden below one year of age (e.g. catch-up at age 6 months with five months duration of protection), there is no issue with discounting. When the catch-up program is assumed to prevent burden beyond one year of age (e.g. catch-up at age 8 months with five months duration of protection), the discounted averted burden will be slightly underestimated.

-   In the config file, if the mean and stdev values specified for efficacy against an outpatient case are exactly the same as the values specified for efficacy against a hospitalized case, the uncertainty distribution around efficacy is sampled only once. If the mean or stdev values are different, two independent uncertainty distributions are sampled, one for efficacy against an outpatient case, and one for efficacy against a hospitalized case.


## Published work with this framework

-   Li, Roberfroid, Bilcke, Castanares-Zapatero, De Meester, Mao, Thiry, Willem, Beutels. Cost-effectiveness of Abrysvo® and Beyfortus® against RSV infections in Belgian infants. Health Technology Assessment (HTA) Brussels: Belgian Health Care Knowledge Centre (KCE). 2025. KCE Reports 402. D/2025/10.273/12. 

-   Li, Hodgson, Flaig, Kieffer, Herring, Beyhaghi, Willem, Jit, Bilcke, Beutels (2023). Cost-Effectiveness of Respiratory Syncytial Virus Preventive Interventions in Children: A Model Comparison Study. Value in Health. 26(4), 508--518.

-   Getaneh, Li, Mao, Johannesen, Barbieri, van Summeren, Wang, Tong, Baraldi, Phijffer, Rizzo, van Wijhe, Heikkinen, Bont, Willem, Jit, Beutels, Bilcke, RESCEU investigators (2023). Cost-effectiveness of monoclonal antibody and maternal immunization against respiratory syncytial virus (RSV) in infants: Evaluation for six European countries. Vaccine. 41(9), 1623-1631.

-   Li, Bilcke, Vázquez Fernández, Bont, Willem, Wisløff, Jit, Beutels, RESCEU investigators (2022). Cost-effectiveness of respiratory syncytial virus disease prevention strategies: Maternal vaccine versus seasonal or year-round monoclonal antibody program in Norwegian children. The Journal of Infectious Diseases, 226 (Supplement 1), S95--S101.

-   Li, Willem, Antillon, Bilcke, Jit, Beutels. (2020) Health and economic burden of Respiratory Syncytial Virus (RSV) disease and the cost-effectiveness of potential interventions against RSV among children under 5 years in 72 Gavi-eligible countries. BMC Medicine. 18, 1--16.

## Other references

-   Cromer et al. (2017) Burden of paediatric respiratory syncytial virus disease and potential effect of different immunisation strategies: a modelling and cost-effectiveness analysis for England, Lancet Public Health, 2017

-   Shi et al (2017) Global, regional, and national disease burden estimates of acute lower respiratory infections due to respiratory syncytial virus in young children in 2015: a systematic review and modelling study. The Lancet. 390(10098), 946-958.

------------------------------------------------------------------------

Copyright 2025, CHERMID, UNIVERSITY OF ANTWERP

Contact: Lander Willem or Xiao Li

Last update: 2025-12-19

