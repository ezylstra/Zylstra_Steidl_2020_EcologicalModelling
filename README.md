# A Bayesian state-space model for seasonal growth of terrestrial vertebrates

### Erin R. Zylstra and Robert J. Steidl

### In Revision

### Code and Data DOI: 
Will add once accepted

### Please contact the first author for questions about the code or data: Erin R. Zylstra (erinzylstra@gmail.com)
_______________________________________________________________________________________________________________________________________

## Abstract:
Will add once accepted
## Code 
1. [Seasonal_Growth.R](Seasonal_Growth.R): R code used to model seasonal growth of canyon treefrogs (*Hyla arenicolor*) in southern Arizona based on 3 years of capture-recapture data.  Includes code used to reproduce Fig. 3 in the manuscript.
2. [Simulations_SeasonalGrowth.R](Simulations_SeasonalGrowth.R): R code used to simulate lengths of canyon treefrogs and estimate growth parameters based on von Bertalanffy models (seasonal and non-seasonal). Includes code used to reproduce Figs. 1 and 2 in the manuscript.  

## Data
1. [AllFrogs.csv](AllFrogs.csv): Capture-recapture data for treefrogs captured in southern Arizona, 2014-2016.  Each row (n = 2778) represents one capture.  The columns are:
    - CapDate: date of capture
    - SVL: snout-vent length, in mm
    - Tag: unique alphanumeric tag number
    - Occasion: survey occasion, numbered 1-11 (three in spring 2014, two in fall 2014, three in spring 2015, two in fall 2015, and one in spring 2016).
    - Sex: sex assignment (f=female/m=male/u=unknown), using information from all captures of each individual (i.e., if sex was unknown when captured with SVL = 37, but assigned as female when later captured with SVL = 41, sex relabeled as female for earlier capture).
    - Sex40: sex assignment (f/m/u), with all individuals listed as “u” when SVL < 40.
