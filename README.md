# Data from heart single-cell RNA-Seq study (Skelly et al. 2018 Cell Reports)

This repo is for data from 
Skelly et al. (2018) Cell Reports [doi: 10.1016/j.celrep.2017.12.072](http://dx.doi.org/10.1016/j.celrep.2017.12.072).
Raw data are available at ArrayExpress [E-MTAB-6173](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6173/).

`full_count_matrix.txt` is the read count matrix 
with 27,998 rows (genes) and 10,563 columns (cells)
and can be downloaded from the ArrayExpress link above 
(file named `E-MTAB-6173.processed.1.zip`).

To read the full count matrix in `R` (after unzipping the file above):
```r
library(data.table)
x <- fread("full_count_matrix.txt")
```

Of all sequenced cells, 10,519 passed basic quality control and filtering 
(see manuscript for details). 
File `cluster_labels.tsv` gives the major cluster to which each 
of these QC-passing cells was assigned.

```r
labs <- fread("cluster_labels.tsv")
```

