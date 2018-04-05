# Analysis for MR Base methods paper

## Run main MR analysis

Perform analysis of LDL, trigs, Lp(a) on all best outcomes

```
cd scripts
Rscript run_mr.R
```

This saves a file in `results/lpa_ldl_trigs.RData` which contains the instruments published by Willer et al, and those variants extracted from a large set of disease and risk factor GWASs.

## Generate graphs

LDL-CHD plots

```
Rscript ldl-chd.R
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

