# WARNING Do not use!!! Not fully functional

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

org<-'mouse' 
npcs<-40 #How many inital PC dimensions to compute. 
k=30 #This for nearest neighbors, 30 is default
source('~/dsdata/projects/R_Packages/scExtras/R/markerGenes_mouse.R')


scrna= RunQC(dir=outdir,org=org,name=projectname,files=input10x ,filter=T, doubletdetection = T,UpperMitoCutoff=10)
scrna = processExper(scrna ,ccscale = T, sc.transform = T)


scrna <- RunPCA(scrna, npcs = 50)
p.elbow <- ElbowPlot(scrna, ndims = 50)

dims<-1:20
scrna <- FindNeighbors(scrna, dims = dims)
scrna <- FindClusters(scrna, resolution = 0.6),algorithm = 2)
scrna <- RunUMAP(scrna, dims = dims)
```

## Process scATAC-Seq using Signac

``` r
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
indir<- 'Adult'
outdir<-'Adult/Signac'
org='mm'

dir.create(outdir)
plotdir <- paste0(outdir,'/Plots')
dir.create(plotdir)


# Load data ---------------------------------------------------------------

counts <- Read10X_h5(paste0(indir,"/outs/filtered_peak_bc_matrix.h5"))

meta <- read.csv(paste0(indir,"/outs/singlecell.csv"),
                 header = TRUE, row.names = 1,
                 stringsAsFactors = FALSE)


assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = paste0(indir,'/outs/fragments.tsv.gz'),
  min.cells = 1
)



atac <- CreateSeuratObject(counts = assay, 
                           assay = "peaks", 
                           project = "ATAC",
                           meta.data = meta)


# Create Annotations ------------------------------------------------------



annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(atac) <- annotations



# QC ----------------------------------------------------------------------

atac <- NucleosomeSignal(object = atac)


atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = atac, group.by = 'nucleosome_group', region = 'chr1-1-10000000')



atac <- TSSEnrichment(atac, fast = FALSE)

atac$high.tss <- ifelse(atac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(atac, group.by = 'high.tss') + NoLegend()


atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

VlnPlot(
  object = atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


atac <- subset(
  x = atac,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
atac




atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(object = atac)


DepthCor(atac)




atac <- RunUMAP(
  object = atac,
  reduction = 'lsi',
  dims = 2:30
)
atac <- FindNeighbors(
  object = atac,
  reduction = 'lsi',
  dims = 2:30
)
atac <- FindClusters(
  object = atac,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

DimPlot(object = atac, label = TRUE) + NoLegend()


# Motif Analysis ----------------------------------------------------------

DefaultAssay(atac) <- 'peaks'

pfm<- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection='CORE',all_versions = FALSE,tax_group='vertebrates')
)



# add motif information
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)


atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)




# Gene Activity Matrix ----------------------------------------------------



gene.activities <- GeneActivity(atac)

# add the gene activity matrix to the Seurat object as a new assay
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)



DefaultAssay(atac) <- 'RNA'
FeaturePlot(
  object = atac,
  features = c('Pdgfra','Pdgfrb'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)



# Integrate RNA/ATAC ------------------------------------------------------

scrna <- readRDS('~/dsdata/projects/Morrisey/Jarod/scRNA/JZ_lung_timeseries/Adult/seuratv4/Seurat.RDS')


scrna <- FindVariableFeatures(
  object = scrna,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = scrna,
  query = atac,
  normalization.method='SCT',
  reduction = 'cca',
  dims = 1:30
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = scrna$var_cluster,
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

atac <- AddMetaData(object = atac, metadata = predicted.labels)


plot1 <- DimPlot(scrna, group.by = 'var_cluster', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(atac, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2



saveRDS(atac,paste0(outdir,'/Signac.RDS'))
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
