Replication of the simulation studies in "Inference for Cause-Specific Cox Model Absolute Risk Estimated from Case-Cohort Data", by Etiévant and Gail (2025). The methods described in the article are implemented in the file `help.functions.R`.

### Required packages 

```
dplyr, ggplot2, grid, gridExtra, gtable, nnet, parallel, survival, xtable.
```

### Scripts

* Script `simulations_exhaustive.R` replicates the simulations proposed by Etiévant and Gail (2024) in Section 5.1.1, 5.2.1, and Web Appendix E.4.1.

* Scripts `simulations_stratified.R` and `simulations_stratified_fraction.R` replicate the simulations in Section 5.1.2, 5.2.2, and Web Appendix E.4.2.

* Script `simulations_exhaustive_Weibull.R` replicates the simulations proposed by Etiévant and Gail (2024) in Web Appendix C.3.1.

* Scripts `simulations_stratified_Weibull.R` and `simulations_stratified_fraction_Weibull.R` replicate the simulations in Web Appendix C.3.2.

Each script relies on functions provided in `help.functions.R`.


### Instructions to run each script

* Save the chosen script(s) and file `help.functions.R` in the same directory.

* Open and run the whole script(s).

* The results of the simulations are saved in csv tables and Rdata files. For example, when running script `simulations_exhaustive.R`, file `details.logAR.y1_exhaustive.csv` will contain the simulation results displayed in Table 1 in Section 5, in Web Table 6 in Web Appendix C.3.1, and in Web Table 11 in Web Appendix E.4.1.


### Functions provided in `help.functions.R`

* **influences** - estimation of influences on the cause-specific log-relative hazard, cause-specific baseline hazard at each unique event time, as well as cause-specific cumulative baseline hazard until each event time of the event of primary interest in interval (t1,t2]. This function should be used with design weights under the case-cohort design with exhaustive sampling of cases, or under the case-cohort design with stratified sampling based on case status when using only the event times of the cases in the case-cohort (Breslow.weight). This function returns parameters estimates.

* **influences.generalized** - estimation of influences on the cause-specific log-relative hazard, cause-specific baseline hazard at each unique event time, as well as cause-specific cumulative baseline hazard until each event time of the event of primary interest in interval (t1,t2]. This function should be used with design weights under the case-cohort design with stratified sampling based on case status when using the event event times of all the cases in the cohort (Breslow.all). This function returns parameters estimates.

* **risk.estimation** - estimation of the absolute risk and pure risk for the primary event of interest.

* **risk.influences** - estimation of influences on the absolute risk and pure risk for the primary event of interest. This function should be used with design weights under the case-cohort design with exhaustive sampling of cases, or under the case-cohort design with stratified sampling based on case status when using only the event times of the cases in the case-cohort (Breslow.weight). This function returns parameters estimates.

* **risk.influences.generalized** - estimation of influences on the absolute risk and pure risk for the primary event of interest. This function should be used with design weights under the case-cohort design with stratified sampling based on case status when using the event event times of all the cases in the cohort (Breslow.all). This function returns parameters estimates.

* **influences.calib** - estimation of influences on the cause-specific log-relative hazard, cause-specific baseline hazard at each unique event time, as well as cause-specific cumulative baseline hazard until each event time of the event of primary interest in interval (t1,t2]. This function should be used with calibrated weights under the case-cohort design with exhaustive sampling of cases. This function returns parameters estimates.

* **influences.generalized.calib** - estimation of influences on the cause-specific log-relative hazard, cause-specific baseline hazard at each unique event time, as well as cause-specific cumulative baseline hazard until each event time of the event of primary interest in interval (t1,t2]. This function should be used with calibrated weights under the case-cohort design with stratified sampling  based on case status (Breslow.all or Breslow.weight). This function returns parameters estimates.

* **risk.influences.calib** - estimation of influences on the absolute risk and pure risk for the primary event of interest. This function should be used with calibrated weights. This function returns parameters estimates.

* **robustvariance** - estimation of the robust variance estimate, i.e., the sum of the squared influences.

* **variance** - estimation of the influence-based variance, that follows the complete variance decomposition. This function should be used with design weights under the case-cohort design with exhaustive sampling of cases.

* **variance.calib** - estimation of the influence-based variance, that follows the complete variance decomposition. This function should be used with calibrated weights under the case-cohort design with exhaustive sampling of cases.

* **variance.generalized** - estimation of the influence-based variance, that follows the complete variance decomposition. This function should be used with design weights or calibrated weights under the case-cohort design with stratified sampling based on case status.

* **conf.interval** - calibration of the design weights using the raking procedure. The Newton Raphson method is used for the optimization.

* **calibration** - calibration of the design weights using the raking procedure.

* **X.generation** - generation of a gaussian variable, with unit variance and given mean.


### Reference

Etiévant L, Gail MH (2025) Inference for Cause-Specific Cox Model Absolute Risk Estimated from Case-Cohort Data. preprint

### Additional scripts

* Script `PLCO analysis.R` replicates the data illustration proposed by Etiévant and Gail (2025) in Section 6.

* Scripts `PLCO supplementary analysis.R` replicates the data illustration in Web Appendix D.

Data underlying these analyses are maintained by the National Cancer Institute, Division of Cancer Epidemiology and Genetics and Division of Cancer Prevention, and are available to bona fide researchers upon submission and approval of a research proposal, and subsequent completion of a Data Transfer Agreement. Proposals can be submitted at https://cdas.cancer.gov/datasets/plco/2/.
