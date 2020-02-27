# A Bayesian state-space model for seasonal growth of terrestrial vertebrates

### Erin R. Zylstra and Robert J. Steidl

### Ecological Modelling 420:108975 [10.1016/j.ecolmodel.2020.108975](https://doi.org/10.1016/j.ecolmodel.2020.108975)

### Code and Data DOI: [10.1016/j.ecolmodel.2020.108975](https://doi.org/10.1016/j.ecolmodel.2020.108975)

### Please contact the first author for questions about the code or data: Erin R. Zylstra (erinzylstra@gmail.com)
_______________________________________________________________________________________________________________________________________

## Abstract:
The rate of somatic growth determines when individuals transition between life stages, which along with survival and reproduction, are principal factors governing the rate of population change. For short-lived species that inhabit seasonally dynamic environments, accounting for fluctuations in somatic growth is necessary to make reliable inferences about population dynamics. We describe a Bayesian, state-space formulation of a von Bertalanffy growth model that integrates a sinusoidal model to allow for seasonal fluctuations in growth while also accounting for individual heterogeneity and measurement error. We use this model to describe post-metamorphic growth of canyon treefrogs, *Hyla arenicolor*, based on capture-recapture data from 404 individuals over a two-year period. For simulated data where we assumed growth varies seasonally, our model provides unbiased estimates of growth rate, mean asymptotic length, standard deviation of individual asymptotic lengths, and date of maximum growth. For field data from canyon treefrogs, we found strong evidence of seasonal variation in growth, with maximum growth during the summer monsoon season. Growth rate of females was 19% lower than males, although on average, females reached asymptotic lengths that were 15% greater than males. Ignoring systematic intra-annual variation in growth can bias inferences about population dynamics, particularly for fast-growing species that reproduce more than once per year or inhabit environments with strong seasonal signals. We present a straightforward approach for using repeated length measurements from individuals of unknown age to estimate growth while accounting for seasonality and individual heterogeneity, which are sources of variation common in many vertebrate populations.
## Code 
1. [Seasonal_Growth.R](Seasonal_Growth.R): R code used to model seasonal growth of canyon treefrogs (*Hyla arenicolor*) in southern Arizona based on three years of capture-recapture data.  Includes code used to reproduce Fig. 3 in the manuscript.
2. [Simulations_SeasonalGrowth.R](Simulations_SeasonalGrowth.R): R code used to simulate lengths of canyon treefrogs and estimate growth parameters based on von Bertalanffy models (seasonal and non-seasonal). Includes code used to reproduce Figs. 1 and 2 in the manuscript.

## Data
1. [AllFrogs.csv](AllFrogs.csv): Capture-recapture data for treefrogs captured in southern Arizona, 2014-2016.  Each row (n = 2778) represents one capture.  The columns are:
    - CapDate: date of capture
    - SVL: snout-vent length, in mm
    - Tag: unique alphanumeric tag number
    - Occasion: survey occasion, numbered 1-11 (three in spring 2014, two in fall 2014, three in spring 2015, two in fall 2015, and one in spring 2016).
    - Sex: sex assignment (f=female/m=male/u=unknown), using information from all captures of each individual (i.e., if sex was unknown when captured with SVL = 37, but assigned as female when later captured with SVL = 41, sex relabeled as female for earlier capture).
    - Sex40: sex assignment (f/m/u), with all individuals listed as “u” when SVL < 40.
