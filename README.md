# ESSD-DVQR

D-vine copula based postprocessing (DVQR) scripts for the ESSD benchmark. 

The code is provided as supplementary material with

* Demaeyer, J., Bhend, J., Lerch, S., Primo, C., Van Schaeybroeck, B., Atencia, A., Ben Bouallègue, Z., Chen, J., Dabernig, M., Evans, G., Faganeli Pucer, J., Hooper, B., Horat, N., Jobst, D., Merše, J., Mlakar, P., Möller, A., Mestre, O., Taillardat, M., and Vannitsem, S.: The EUPPBench postprocessing benchmark dataset v1.0, Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2022-465, in review, 2023.

**Please cite this article if you use (a part of) this code for a publication.**

# Method

In the **D-vine (drawable vine) copula based postprocessing**, a multivariate conditional copula $C$ is estimated using a pair-copula construction for the graphical D-vine structure according to [Kraus and Czado, (2016)](https://arxiv.org/pdf/1510.04161.pdf). D-vine copulas enable a flexible modelling of the dependence structure between the observation $y$ and the ensemble forecast $x_1, \ldots, x_m$ (see, e.g. [Möller et al., (2018)](https://arxiv.org/pdf/1811.02255.pdf)). The covariates $x_1, \ldots, x_m$ are selected by their predictive strength based on the conditional log-likelihood. Afterwards, the D-vine copula quantile regression (DVQR) allows to predict quantiles $\alpha\in (0,1)$ that represent the postprocessed forecasts via

$$F^{-1}_{y\vert x_1, \ldots, x_m}(\alpha\vert x_1(t), \ldots, x_m(t)):=F_y^{-1}\left(C^{-1}(\alpha\vert F_{x_1}(x_1(t)),\ldots, F_{x_m}(x_m(t)))\right),$$    

where $F_{x_i}$ denote the marginal distributions of $x_i$ for all $i=1,\ldots, m$, $F_{y}^{-1}$ the inverse marginal distribution of $y$ and $C^{-1}$ the conditional copula quantile function. The marginal distributions will be estimated by kernel densities in our case. DVQR is estimated separately for every station and lead time using a seasonal adaptive training period.

# Implementation details

The whole method is implemented in the programming language [R](https://www.r-project.org).

## Data 

- To get the ESSD benchmark data set you can use the [download script](https://github.com/EUPP-benchmark/ESSD-benchmark-datasets). 
- To construct data frames from the NetCDF files of the ESSD benchmark, you can use the the R-script `benchmarkdata.R`. With this script, you receive `benchmark_t2m_train_ext.Rdata`, which denotes the **training data** and `benchmark_t2m_test_ext.Rdata`, which contains the **test/validation data**. These constructed data sets will be used in the ESSD benchmark for this method.

The following table describes variable names that are referred to in training as well as test/validation data:

| Variable | Description |
| ---- | ----------- | 
| `init_date` | Initialization date of the ensemble forecasts. |
| `valid_date` | Valid date of the ensemble forecasts. |
| `leadtime` | Lead time of the ensemble forecasts in hours. |
| `id` | Station id. |
| `name` | Station name. |
| `station_altitude` | Station altitude in meter. |
| `station_longitude` | Station longitude. |
| `station_latitude` | Station latitude. |
| `station_land_usage` | Station land usage. |
| `model_altitude` | Model altitude in meter. |
| `model_longitude` | Model longitude. |
| `model_latitude` | Model latitude. |
| `model_land_usage` | Model land usage. |
| `observation` | Observation of t2m in Kelvin. |
| `ctrl` | Control ensemble forecast of t2m in Kelvin. |
| `ensi` | $i$-th member of t2m ensemble forecasts in Kelvin, $i = 1, ..., 10$ (training data) or $i = 1, ..., 50$ (validation data) |
| `mean` | Mean of the ensemble forecasts of t2m in Kelvin. |
| `q10` | 1st decile of the ensemble forecasts of t2m in Kelvin. |
| `median` | Median of the ensemble forecasts of t2m in Kelvin. |
| `q90` | 9th decile of the ensemble forecasts of t2m in Kelvin. |
| `sd` | Standard deviation of the ensemble forecasts of t2m in Kelvin. |
| `iqr` | Interquartil range of the ensemble forecasts of t2m in Kelvin. |
| `rand` | Random ensemble forecasts of the t2m ensemble forecasts in Kelvin. |

## D-vine copula based postprocessing

For the D-vine copula based postprocessing we applied the R-script `dvine_pp.R`.

- A **local postprocessing** is selected, i.e. the locations are postprocessed (for each lead time) separately.
- A **weekly training period** is used to account for seasonality, i.e. a symmetric interval of `nweeks` around the week where the forecast day belongs to is chosen. For example, if you would like to forecast any day in the 8th week of the year and you set `nweeks = 3`, you would take the data of the **5th, 6th, 7th**, <u>8th</u>, **9th, 10th, 11th** week in the **training data** `tdata` to estimate the model for the 8th week. This yields in total 52 estimated models for each station having a fixed lead time, i.e. for each week of the year you have one model (week 53 is set to week 52). For our final postprocessing models we select `nweeks = 2`, which is determined by the in-sample mean CRPS on the training data over all locations and `leadtimes`. Consequently we have for  each week, each location and each lead time a model which is used to produce the postprocessed quantile ensemble forecasts in the **validation data** `vdata`. 
- The **t2m observation** is selected as response variable for the training data, i.e. `tresponse = "observation"`.
- The **control and mean ensemble forecast of t2m** is provided as predictor variable for the training and validation data, i.e. `tpredictors = c("ctrl", "mean")` and `vpredictors = c("ctrl", "mean")`, respectively.
- **All copula families** in the R-package [vinereg](https://github.com/tnagler/vinereg) are available, i.e. `family_set = "all"`.
- As forward covariate selection criterion the **AIC-corrected conditional log-likelihood** is selected, i.e. `selcrit = "aic"`.
- The **order** of the covariates in the D-vine is automatically determined by the algorithm itself, i.e. `order = NA`.
- `ffd` denotes the **first forecast date** and `lfd` the **last forecast date**. 
- The parameter `cores` specifies the **number of cores used for code parallelization**. 

> **Remark:**
> The computational effort of this method could be enormously lessened by e.g. reducing the `family_set` to a smaller subset of copulas or by the pre-selection of the variable `order`. 

## Used R-packages
- [vinereg](https://github.com/tnagler/vinereg): For the D-vine copula quantile regression proposed by [Kraus and Czado, (2016)](https://arxiv.org/pdf/1510.04161.pdf).
- [lubridate](https://github.com/tidyverse/lubridate): For calculating the week in a year.
- [doMC](https://cran.r-project.org/web/packages/doMC/doMC.pdf): For multi-core usage.
- [ncdf4](https://cran.r-project.org/web/packages/ncdf4/index.html): For reading the NetCDF-files in R.

