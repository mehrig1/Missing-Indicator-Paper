# Component-wise censoring paper
Data simulation and analysis for assessing performance of Cox models for composite endpoints with component-wise censoring

Journal: Clinical Trials 

Paper link: https://doi.org/10.1177/17407745231177046

Title: Performance of Cox regression models for composite time-to-event endpoints with component-wise censoring in randomized trials 

Authors: Authors: Speiser, Jaime Lynn; Ambrosius, Walter T; Pajewski, Nicholas M.
Department of Biostatistics and Data Science, Wake Forest University School of Medicine, Winston-Salem, NC, USA.

Funding: National Institute on Aging, PRagmatic EValuation of evENTs And Benefits of Lipid-lowering in oldEr adults (PREVENTABLE) trial, Grant number: U19 AG065188. Wake Forest Older Americans Independence Center, Grant Number: P30 AG021332.

Trial website: www.preventabletrial.org

Abstract: Composite time-to-event endpoints are beneficial for assessing related outcomes jointly in clinical trials, but components of the endpoint may have different censoring mechanisms. For example, in the PRagmatic EValuation of evENTs And Benefits of Lipid-lowering in oldEr adults (PREVENTABLE) trial, one part of the composite outcome is right censored (all-cause mortality), whereas the other part is interval censored (dementia and persistent disability). Although Cox regression is an established method for time-to-event outcomes, it is unclear how models perform under a component-wise censoring scheme for large clinical trial data. The goal of this paper is to conduct a simulation study to investigate the performance of Cox models under different scenarios for composite endpoints with component-wise censoring. We simulate data by varying the strength and direction of the association between treatment and outcome for the two components, the proportion of events arising from the two components of the outcome (right censored and interval censored), and the method for including the interval censored component in the Cox model (upper value and midpoint of the interval). Under these scenarios, we compare the treatment effect estimate bias, confidence interval coverage, and power. Based on the simulation study, Cox models generally have adequate power to detect a significant treatment effect for composite outcomes with component-wise censoring. Performance was similar regardless of if the upper value or midpoint of the interval censored part of the composite outcome was used. 

## Code files
* data simulation midpt analysis: provides functions for simulating data that use the midpoint of the interval censored component in the composite outcome and aggregating results for different scenarios across the simulation runs
* data simulation upper value analysis: provides functions for simulating data that use the upper value of the interval censored component in the composite outcome and aggregating results for different scenarios across the simulation runs
* comparing midpt and upper results: provides code that reads in results files for the different scenarios and produces Bland Altman plots to compare the differences between use of midpoint of the interval versus use of the upper value of the interval in the calculation for the composite outcome
