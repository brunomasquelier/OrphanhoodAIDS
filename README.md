# OrphanhoodAIDS
In countries without adequate death registration systems, adult mortality is often estimated using orphanhood-based methods. The HIV pandemic breaches several assumptions of these methods, for example, by increasing the correlation between maternal and child survival. Using microsimulations we generated 1,152 populations facing HIV epidemics and evaluated different orphanhood-based estimates against the underlying mortality rates. We regressed survivorship probabilities on proportions of respondents with surviving mothers, adjusting for trends in seroprevalence and coverage of antiretroviral therapy, to obtain new coefficients. We tested the different methods on survey and census data from 16 African countries with high HIV prevalence. The method is detailed in 

> Masquelier, B. and Timæus, I.M, Estimating adult mortality based on maternal orphanhood in populations with HIV/AIDS, Population Studies, _forthcoming_

A preprint is available here:
https://blogs.lshtm.ac.uk/iantimaeus/files/2024/05/MasquelierTimaeusOrphanhoodPreprint.pdf

This repository contains the following folders:
- R Code for microsimulations: R Code to recreate the set of microsimulations. Users interested in the microsimulations need a local installation of SOCSIM, which is available here: https://lab.demog.berkeley.edu/socsim/. It can also be revised to use the rsocsim package in R, available here: https://github.com/MPIDR/rsocsim.
- Data: Mortality and fertility standards needed to format the demographic rates to inform the microsimulations. 
- Results: CSV data files with estimates of nq25 from the different orphanhood approaches.
- Workbook: An Excel workbook to facilitate the calculations for practical applications based on tabulated data.
