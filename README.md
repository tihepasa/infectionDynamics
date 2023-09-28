# infectionDynamics

Supplementary material for the paper *Spatio-temporal modeling of co-dynamics of smallpox, measles and pertussis in pre-healthcare Finland*.

## Files
The technical model specification in Stan is done in the file `infection_model.stan`, and the R code to fit the model with cmdstanr can be found from the file `infection_model.R`. The data used for the analysis is in the file `infectionpars.rds`. The data include the monthly, regional, dichotomous death occurrences from January 1820 to December 1850, and some additional parameters needed to fit the model. A more specific description of the contents of the parameter file can be found from the beginning of the file `infection_model.R`.

## Figures
![image](https://github.com/tihepasa/infectionDynamics/assets/65618755/8a1096cb-c8c1-417d-b990-468ed32e5aba)


![image](https://github.com/tihepasa/infectionDynamics/assets/65618755/c7d71380-b501-4d72-afb6-2fd8c17c84bd)

## Tables
![image](https://github.com/tihepasa/infectionDynamics/assets/65618755/914dce2f-929c-496b-9f8e-fe88a0546c17)
