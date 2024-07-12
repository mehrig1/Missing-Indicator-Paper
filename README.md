# Missing Indicator Method paper
Data simulation and analysis for assessing whether the missing indicator method is beneficial in longitudinal data modeling 

Title: Imputation and missing indicators for handling missing longitudinal data: A simulation study based on electronic health record data

Authors: Ehrig, M; Bullock, GS; Leng, X; Pajewski, NM; Speiser, JL
Department of Biostatistics and Data Science, Wake Forest University School of Medicine, Winston-Salem, NC, USA.

Funding: This project was supported in part by the NIH/NLM (R25 LM014214) in the Department of Biomedical Engineering and Center for Biomedical Informatics at Wake Forest University School of Medicine. Dr. Speiser is supported by the National Institute On Aging of the National Institutes of Health under Award Number K25AG068253. This study was supported in part by the Wake Forest Claude D. Pepper Older Americans Independence Center (P30 AG021332). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

Abstract: 
Introduction: Missing data in electronic health records (EHRs) is highly prevalent and results in analytical concerns such as heterogeneous sources of bias and loss of statistical power. One simple analytic method for addressing missing or unknown covariate values is to treat missing-ness for a particular variable as a category onto itself, which we refer to as the missing indicator method.  For cross-sectional analyses, recent work suggested that there was minimal benefit to the missing indicator method; however, it is unclear how this approach performs in the setting of longitudinal data, in which correlation among clustered repeated measures may be leveraged for potentially improved model performance.  

Methods: We conducted a simulation study aimed to evaluate whether the missing indicator method improved model performance and imputation accuracy for longitudinal data mimicking an application of developing a clinical prediction model for falls in older adults based on EHR data. We simulated a longitudinal binary outcome using mixed effects logistic regression that emulated a falls assessment at annual follow-up visits. Using multivariate imputation by chained equations, we simulated time-invariant predictors such as sex and medical history, as well as dynamic predictors such as physical function, body mass index, and medication use. We induced missing data in predictors under scenarios that had both random (MAR) and dependent missing-ness (MNAR). We evaluated aggregate performance using the area under the curve for models with and without missing indicators as predictors, as well as complete case analysis, across simulation replicates. We evaluated imputation quality using normalized root mean square error for continuous variables, and percent falsely classified for categorical variables.  

Results: Independent of the mechanism used to simulate missing data (MAR or MNAR), overall model performance via area under the curve was similar regardless of whether missing indicators were included in the model. The root mean square error and percent falsely classified measures were similar for models including missing indicators versus those without missing indicators. Model performance and imputation quality were similar regardless of whether the outcome was related to missingness. Imputation with or without missing indicators had similar mean values of area under the curve compared to complete case analysis, although complete case analysis had the largest range of values. 

Discussion: The results of this study suggest that the inclusion of missing indicators in longitudinal data modeling neither improve nor worsen overall performance or imputation accuracy.  Future research is needed to address whether the inclusion of missing indicators is useful in prediction modeling with longitudinal data in different settings, such as high dimensional data analysis.
 

## Code files
* MAR draft script v2 20 - MAR missing mechanism, 20% missing data, indicators not included when simulating the outcome; creates "MAR v2 20.csv"
* MAR draft script v2 50 - MAR missing mechanism, 50% missing data, indicators not included when simulating the outcome; creates "MAR v2 50.csv"
* MAR draft script v3 20 - MAR missing mechanism, 20% missing data, indicators included when simulating the outcome; creates "MAR v3 20.csv"
* MAR draft script v3 50 - MAR missing mechanism, 50% missing data, indicators included when simulating the outcome; creates "MAR v3 50.csv"
* MNAR draft script v2 20 - MNAR missing mechanism, 20% missing data, indicators not included when simulating the outcome; creates "MNAR v2 20.csv"
* MNAR draft script v2 50 - MNAR missing mechanism, 50% missing data, indicators not included when simulating the outcome; creates "MNAR v2 50.csv"
* MNAR draft script v3 20 - MNAR missing mechanism, 20% missing data, indicators included when simulating the outcome; creates "MNAR v3 20.csv"
* MNAR draft script v3 50 - MNAR missing mechanism, 50% missing data, indicators included when simulating the outcome; creates "MAR v3 50.csv"
* Results script - reads in the results csv files and creates 3 figures found in the paper
* Results tables - reads in the results csv files and creates supplementary tables
* Supplemental document updated - creates the supplemental document, which includes summary and descriptive statistics of one simulated dataset
