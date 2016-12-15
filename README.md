# LDL analysis for MR Base methods paper


## Run main MR analysis

Perform analysis of LDL, trigs, Lp(a) on all best outcomes

```
cd scripts
Rscript run_mr.R
```

This saves a file in `results/lpa_ldl_trigs.RData` which contains the instruments published by Willer et al, and those variants extracted from a large set of disease and risk factor GWASs.


## Run analysis of drug proxies on outcomes

```
Rscript run_mr_drug_proxies.R
```

Saves the results to `results/drug_proxies.RData`


## Coordinate RCTs and MR drug prediction

```
Rscript drug_prediction.R
```

Saves the results to `results/drug_prediction.RData`

## Generate graphs

Drug predictions:

```
Rscript figure_drug_prediction.R
```

Forest plots:

```
Rscript figure_forest_plots.R
```

Sensitivity analyses for all exposure - outcome pairs with p < 0.05:

```
Rscript figure_sensitivity_boxes.R
```

Volcano plots:

```
Rscript figure_volcano_plots.R
```

## Create tables

```
Rscript tables.R
```

Generates

- Drug prediction vs RCT
- Full MR results
- Instrument data

## Notes

- Results of RCTs for drugs on CHD and T2D are in `data/drug_trials.csv`
- List of GLGC SNPs for all lipids are in `data/glgc_snps.txt`

