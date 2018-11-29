## Obtaining mouse orthologs of ligand-receptor pairs

```r
library(tidyverse)
library(readxl)
library(biomaRt)
```

Read in supp table from Ramilowski et al.
```r
dat <- read_excel("ncomms8866-s3.xlsx", sheet="All.Pairs") %>%
  dplyr::filter(!grepl("EXCLUDED", Pair.Evidence)) %>%
  dplyr::select(Pair.Name, Ligand.ApprovedSymbol, Receptor.ApprovedSymbol) %>%
  mutate(pairnum=1:n()) %>%
  dplyr::rename("pair"="Pair.Name", "ligand"="Ligand.ApprovedSymbol",
         "receptor"="Receptor.ApprovedSymbol") %>%
  gather(key="type", value="gene", ligand:receptor)
```

Get human-mouse orthologs
```r
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
mouse <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
genesmus <- getLDS(attributes="hgnc_symbol", filters="hgnc_symbol",
  values=dat$gene, mart=human, attributesL="mgi_symbol", martL=mouse,
  uniqueRows=TRUE)
```

Select 1:1 orthologs only
```r
genesmus2 <- rename(genesmus, "human"="HGNC.symbol",
  "mouse"="MGI.symbol") %>% group_by(human) %>%
  mutate(nmouse=n()) %>% group_by(mouse) %>%
  mutate(nhuman=n()) %>%
  filter(nmouse == 1, nhuman == 1) %>%
  dplyr::select(-nmouse, -nhuman)
```

Get ENSEMBL gene IDs for mouse
```r
ens <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "chromosome_name"),
  filter="mgi_symbol", values=genesmus2$mouse, mart=mouse)
genesmus3 <- inner_join(genesmus2, ens, by=c("mouse"="mgi_symbol")) %>%
  dplyr::filter(!grepl("PATCH", chromosome_name)) %>%
  dplyr::filter(!grepl("CHR_MM", chromosome_name)) %>%
  dplyr::select(-chromosome_name)

```

Compile a dataset of receptor-ligand pairs in mouse
```r
pairgenes <- sort(unique(genesmus3$ensembl_gene_id))
pairs <- inner_join(dat, genesmus3, by=c("gene"="human")) %>%
  mutate(mouse_ens=paste0(mouse, "~", ensembl_gene_id)) %>%
  dplyr::select(-mouse, -ensembl_gene_id) %>%
  group_by(pair) %>% filter(n() == 2) %>%
  dplyr::select(-gene) %>%
  spread(type, mouse_ens) %>%
  separate(ligand, into=c("ligand_symbol", "ligand_gene_id"), sep="~") %>%
  separate(receptor, into=c("receptor_symbol", "receptor_gene_id"), sep="~")
```
