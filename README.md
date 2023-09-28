# infectionDynamics

Supplementary material for the paper *Spatio-temporal modeling of co-dynamics of smallpox, measles and pertussis in pre-healthcare Finland*.

## Files
The technical model specification in Stan is done in the file `infection_model.stan`, and the R code to fit the model with cmdstanr can be found from the file `infection_model.R`. The data used for the analysis is in the file `infectionpars.rds`. The data include the monthly, regional, dichotomous death occurrences from January 1820 to December 1850, and some additional parameters needed to fit the model. A more specific description of the contents of the parameter file can be found from the beginning of the file `infection_model.R`.

## Figures
### Figure 1
![plot](./figures/figure1_b.png)

### Figure 2
![plot](./figures/figure2_c.png)

## Tables
### Table 1
![plot](./figures/table1.png)
