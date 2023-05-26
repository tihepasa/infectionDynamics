# infectionDynamics

Supplementary material for the paper *Spatio-temporal co-dynamics of smallpox, measles and pertussis in pre-healthcare Finland*.

The folder `figures` includes maps illustrating the regional adjustment parameters.

The model is technically defined in the file `infection_model.stan`, and the file `infection_model.R` uses the definition to fit the model. The data for the model is in the file `infectionpars.rds` which is read into R in the file `infection_model.R`. The data include the monthly, regional, dichotomous death occurrences from January 1820 to December 1850, and some additional parameters needed to fit the model. A more specific description of the contents of the parameter file can be found in the file `infection_model.R`.
