<<<<<<< HEAD
# WARNING Do not use!!! Not fully functional
=======
<span style="color: red;">\*\*WARNING Do not use!!! Not fully
functional</span>
>>>>>>> 3d45abd0c95ce4177165719ba251a4de03d765aa

# SCAPE 0.1.0

SCAPE is a tool that will create a network like visualization of the
following interaction between 2 cells/cluster using scRNA-seq and
scATAC-seq data.

1.  Transcription factor (TF) to ligand gene
2.  Ligand to Rececptor

## Install

``` r
install.packages("devtools")
devtools::install_github("mpmorley/SCAPE")
```

### Requiremnts

cytoscape

## Process scRNA-seq using Seurat

``` r
library(Seurat)
```

## Process scATAC-Seq using Signac

``` r
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
```

## Run SCAPE

``` r
rna.object <- readRDS("Seurat.RDS")
atac.object <- readRDS("Signac.RDS")
ligand.cell=10 #Amp
receptor.cell=8 #MANC


recligDB <- readr::read_csv('~/dsdata/NGSshare/FANTOM5/Mm_PairsLigRec.csv')
tf <- readr::read_tsv('~/dsdata/NGSshare/AnimalTFDB/Mus_musculus_TF_v3.txt')
mm2hs <- readr::read_csv('~/dsdata/NGSshare/homologs/mouse_human.csv')


net <- runSCAPE(atac.object,rna.object,ligand.cell=10,receptor.cell = 8)
```

We Can now visualize the netwrok using the cytoscape program

``` r
library(RCy3)
cytoscapePing()
createNetworkFromIgraph(net,"myIgraph")
```

Apply a custom layout for SCAPE

``` r
addSCAPEStyle()
```
