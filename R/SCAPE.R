library(Signac)
library(Seurat)
library(tidyverse)

library(GenomeInfoDb)

library(ggplot2)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)

library(ChIPpeakAnno)


library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)

library(RCy3)



#'RunScape
#' @param atac.object ATAC Seurat object with a peaks assay
#' @param rna.object RNASeurat object
#' @param signal.cell Ident value for cell that will singal
#' @param receive.cell Ident value for cell that will singal
#' @param min.pct
#' @param annoDB A list of tables of various annotation. See CreateAnnoDB
#' @import dplyr tidyr Seurat
#' @export

runSCAPE <- function(atac.object=NULL,
                     rna.object=NULL,
                     ligand.cell=NULL,
                     receptor.cell=NULL,
                     min.pct=.25,
                     annoDB=NULL
){




# Compute Diffexpression --------------------------------------------------



diffexp <- FindMarkers(object = rna.object,
            ident.1 =  ligand.cell,
            ident.2 = receptor.cell,
            test.use='MAST') %>%
  rownames_to_column('gene') %>%
  mutate(absdiff = abs(pct.1 - pct.2))

#Create a list of "expresss genes for each singal and reciever cells
ExpressGenes <- list(
  'ligand.cell' = GetExpressedGenes(rna.object,group =ligand.cell ,min.pct=min.pct),
  'receptor.cell' = GetExpressedGenes(rna.object,group =receptor.cell ,min.pct=min.pct)
)


tflist <- diffexp  %>% dplyr::filter(gene %in% tf$Symbol)
#whichdiff peaks to use.. perhaps caompare to all cells.
DefaultAssay(atac.object) <- 'peaks'
peaks <- FindMarkers(
  object = atac.object,
  ident.1 =  ligand.cell,
  ident.2 = receptor.cell,
  min.pct = .25,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)



sigpeaks <- peaks %>% tibble::rownames_to_column('peak') %>% dplyr::filter(p_val_adj < 0.05) %>% pull(peak)


peaks.ann <- annotatePeakInBatch(StringToGRanges(sigpeaks, sep=c(":","-")),
                                     AnnotationData=Annotation(atac.object),
                                     #output="nearestBiDirectionalPromoters",
                                     output="overlapping",
                                     #output="upstream&inside",
                                     #maxgap=2000,
                                     bindingType='fullRange',
                                     bindingRegion=c(-2500, 1)) %>%
unname()%>% as.data.frame(.) %>% mutate(peakid = paste0(seqnames,':',start,'-',end))




motifs <- FindMotifs(object = atac.object,
                     features =sigpeaks) %>%
  mutate(p_val_adj = p.adjust(pvalue),motif.name = gsub('\\(var.2\\)','',motif.name)) %>%
  tidyr::separate(motif.name,'motif.name',remove=T) %>%
  left_join(., mm2hs, by=c('motif.name' = 'human_name')) %>%
  mutate(motif.name=ifelse(is.na(mouse_name),motif.name,mouse_name)) %>% dplyr::select(-human_id:-mouse_id) %>%
  dplyr::filter(motif.name %in%  ExpressGenes$ligand.cell & p_val_adj < 0.05 & fold.enrichment > 1.5) %>%
  group_by(motif.name) %>% top_n(1)

nodes <- list()
edges <- list()



#make Nodes for the cells.
nodes$cell <- map(c(ligand.cell,receptor.cell),function(c){
  data.frame(name=as.character(c),type='celltype')
}) %>% bind_rows()



# Cell To TF --------------------------------------------------------------

#Nodes for TFs
nodes$tf <- pmap(motifs,function(motif,motif.name,...){
  data.frame(name=motif.name,type='TF')
}
) %>% bind_rows()

#Cell to TF edges
edges$cell2tf <- pmap(motifs,function(motif.name,...){
  data.frame(from=as.character(ligand.cell),to=motif.name,edgetype='cell2TF')
}
) %>% bind_rows()



# TFs to Ligands ----------------------------------------------------------

edges$tf2lig <- pmap(motifs,function(motif,motif.name,...){
  motif.peaks <- getMotifPeaks(atac.object,motif)
  peaks.ann %>% dplyr::filter(peakid %in% motif.peaks & !is.na(gene_name) & gene_name %in% recligDB$ligand & gene_name %in% ExpressGenes$ligand.cell) %>%
    mutate(from=motif.name,to=gene_name,edgetype='TF2lig') %>% dplyr::select(from,to,edgetype) %>% distinct()
}
) %>% bind_rows()


nodes$lig <- pmap(edges$tf2lig,function(to,...){
  data.frame(name=to,type='ligand')
}) %>% bind_rows() %>% distinct()



# Lig to Rec --------------------------------------------------------------


edges$lig2rec <-
  pmap(nodes$lig,function(name,...){
 name <- enquo(name)
 recligDB %>%  dplyr::filter(ligand %in% !!name & receptor %in% ExpressGenes$receptor.cell) %>%
   mutate(from=ligand,to=receptor,edgetype='lig2rec') %>%
   dplyr::select(from,to,edgetype) %>% distinct()

}
) %>% bind_rows()

nodes$receptor <- pmap(edges$lig2rec,function(to,...){
  data.frame(name=to,type='receptor')
}) %>% bind_rows() %>% distinct()



# Rec to cell -------------------------------------------------------------
edges$rec2cell <-
  pmap(nodes$receptor,function(name,...){
    data.frame(from=name,to=as.character(receptor.cell),edgetype='rec2cell')
  }
  ) %>% bind_rows()

E <- bind_rows(edges) %>% distinct()
N <- bind_rows(nodes) %>% distinct()


net <- igraph::graph_from_data_frame(d=E, vertices=N, directed=T)
net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)

return(net)

}


cytoscapePing()
createNetworkFromIgraph(net,"myIgraph")

addSCAPEStyle <- function(node.col=c('#8804e1','#e1210c','#2074f0','#ff9800')){
  style.name = "myStyle"
  defaults <- list(
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="N,S,c,0.00,0.00")
  nodeLabels <- mapVisualProperty('node label','id','p')
  nodeColors <- mapVisualProperty('node fill color','type','d',c('celltype','TF','ligand','receptor'), node.col)
  nodeShape <- mapVisualProperty('node shape','type','d',c('celltype','TF','ligand','receptor'),c('ELLIPSE','OCTAGON','DIAMOND','VEE'))


#and then create the style
createVisualStyle(style.name, defaults, list(nodeLabels,nodeColors,nodeShape))

#finsh by applying the style
setVisualStyle(style.name)
}





# getMotifPeaks -----------------------------------------------------------

getMotifPeaks <- function(object,motif){
  p <-GetMotifData(object)[,motif]==1
  sub('-',':',names(p)[which(p)])
}


GetExpressedGenes <- function(object,group,min.pct){
  cells <- WhichCells(object,idents=group)
  m <- Matrix::rowSums(x = object@assays$RNA@counts[,cells, drop = FALSE] > 0) /length(x = cells)
  names(m)[which(m>min.pct)]
}



