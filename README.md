# Metabolic Subphenotype Predictor

This repository includes code used in the Type 2 Diabetes Metabolic Phenotyping study. The project aims to address the heterogeneity in glucose metabolism by classifying individuals into distinct metabolic sub-phenotypes (muscle insulin resistance, beta-cell dysfunction, impaired incretin action, and hepatic insulin resistance) based on data from oral glucose tolerance tests (OGTTs) performed with continuous glucose monitoring (CGM) at the clinical research unit and at-home.


![](/summary_figure.png)


## This repository contains code for:
- Inference of T2D metabolic subphenotypes (MuscleIR, Beta-cell Function, Incretin Effect, Hepatic IR).
- Identification of dominant metabolic subphenotype.
- Feature extraction from glucose tiemseries.
- Extraction of reduced representation of glucose tiemseries
- Visualization of metabolic phenotypes based on various glucose-related metrics.
- Concordance between CGM and Venous glucose values from at home and at clinical setting
- Classification of metabolic subphenotypes.

## Prerequisites
#### R (version >= 4.4) and the following R packages:
- reshape2, ggplot2, data.table, dplyr, MESS, imputeTS, plyr, pheatmap, corrplot, Hmisc, corrr, PerformanceAnalytics, ggrepel, goeveg, scales

#### Python (version >= 3.7) and the following Python packages:
- numpy, pandas, matplotlib, scikit-learn, os, datetime, random, math, skfda, xgboost



## License
This project is licensed under the MIT License - see the LICENSE file for details.

For any inquiry, please contact ametwall@stanford.edu 



