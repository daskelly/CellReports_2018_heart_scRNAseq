## Brief overview of cell-cell communication analysis

In the main text we have the following:

> To define intercellular communication networks within the cardiac cellulome,
> we used a dataset of human ligand-receptor pairs (Ramilowski et al., 2015)
> to develop a list of mouse orthologs comprising 2,009 ligand-receptor pairs
> (Table S3; Experimental Procedures).

The supplement includes the following text:

> To quantify the potential for cell-cell communication in the cardiac cellulome,
> we obtained mouse orthologs of human ligand-receptor pairs compiled by Ramilowski et al.
> (Ramilowski et al., 2015) ... We defined a ligand or receptor as "expressed" in a
> particular cell type if 20% of the cells of that type had non-zero read counts for
> the gene encoding the ligand/receptor. To define networks of cell-cell communication,
> we linked any two cell types where the ligand was expressed in the former cell type
> and the receptor in the latter.

Please note the caveats which are either raised directly or implied in the paper,
including:

 * an inability to model anatomical barriers between cell types
 * 20% of cells having non-zero read counts is just one possible way to define receptors/ligands as "expressed"
 * transcription does not necessarily indicate translation and/or functional/properly localized protein

### Overview of basic analysis steps in `R`

First, compile a list of receptor-ligand pairs. For example, you may want to make a 
`data.frame` with columns `receptor` and `ligand`, where each row is a receptor-ligand 
pair. If you have *N* receptor-ligand pairs, this `data.frame` will have *N* rows.
Let's call this `data.frame` of receptor-ligand pairs `pairs`.
We are going to work with the columns `pairs$receptor` and `pairs$ligand`
below. Note that many receptors and ligands are likely to be duplicated
in `pairs$receptor` and in `pairs$ligand`. 
However, the receptor-ligand pair `pairs$receptor[i]`--`pairs$ligand[i]`
should be unique.
If we have a matrix of raw UMI counts with genes in rows and
cells in columns `raw_counts`:

```r
library(tidyverse)
library(assertthat)
raw_ligand <- raw_counts[pairs$ligand, ]
raw_receptor <- raw_counts[pairs$receptor, ]
```

Obtain boolean matrices indicating expression of each receptor/ligand
in each cell type. `pops` is a vector of cell type labels:
```r
expressed <- function(x, cell_pops, threshold=0.2) {
  assert_that(length(x) == length(cell_pops))
  tapply(x, cell_pops, function(y) mean(y > 0) > threshold)
}
ligand_expressed <- apply(raw_ligand, 1, expressed, cell_pops=pops)
receptor_expressed <- apply(raw_receptor, 1, expressed, cell_pops=pops)
```

Finally, tally up the number of interactions between ligands and
receptors expressed in each pair of cell types.
```r
interactions <- tcrossprod(ligand_expressed, receptor_expressed)

idat <- rownames_to_column(as.data.frame(interactions), "ligand") %>%
  gather("receptor", 'value', -ligand)
```

Make a basic heatmap plot to show the numbers.
```r
lev <- names(sort(rowMeans(interactions)))
idat$ligand <- factor(idat$ligand, levels=lev)
idat$receptor <- factor(idat$receptor, levels=lev)
gi <- ggplot(idat, aes(x=ligand, y=receptor)) +
  geom_tile(aes(fill=value)) +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_distiller(palette="Greens") + theme_bw(base_size=16) +
  theme(axis.text.x=element_text(angle=90))
```
