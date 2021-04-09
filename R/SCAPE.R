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


pal <- c('red','blue','green','orange')

#rectf <- readr::read_csv('~/dsdata/NGSshare/FANTOM5/Mm_PairsRecTF.csv')
recligDB <- readr::read_csv('~/dsdata/NGSshare/FANTOM5/Mm_PairsLigRec.csv')

tf <- readr::read_tsv('~/dsdata/NGSshare/AnimalTFDB/Mus_musculus_TF_v3.txt')
mm2hs <- readr::read_csv('~/dsdata/NGSshare/homologs/mouse_human.csv')



rna.object <- readRDS("~/dsdata/projects/Morrisey/Jarod/scRNA/JZ_lung_timeseries/Adult/seuratv4/Seurat.RDS")
atac.object <- readRDS("~/dsdata/projects/Morrisey/JohnLeach/scATAC/Adult/Signac/Signac.RDS")
ligand.cell=10
receptor.cell=8
min.pct = .2



runSCAPE <- function(atac.object=NULL,
                     rna.object=NULL,
                     ligand.cell=NULL,
                     receptor.cell=NULL,
                     min.pct=.25
){



}



# Compute Diff expression for peaks ---------------------------------------


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


<<<<<<< HEAD




=======




>>>>>>> de93b641d163a4fb71c34dbe9fbfffd16f487b28
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
  data.frame(name=as.character(c),type='celltype',size=10,color=pal[1])
}) %>% bind_rows()



#Nodes for TFs
nodes$tf <- pmap(motifs,function(motif,motif.name,...){
  data.frame(name=motif.name,type='TF',size=8,color=pal[2])
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

view(edges$tf2lig)

nodes$lig <- pmap(edges$tf2lig,function(to,...){
  data.frame(name=to,type='ligand',size=8,color=pal[2])
}) %>% bind_rows() %>% distinct()



# Lig to Rec --------------------------------------------------------------
<<<<<<< HEAD


#Perhaps add filter for Source and Evidence.

edges$lig2rec <-
  pmap(nodes$lig,function(name,...){
 name <- enquo(name)
 recligDB %>%  dplyr::filter(ligand %in% !!name & receptor %in% ExpressGenes$receptor.cell) %>%
   mutate(from=ligand,to=receptor,edgetype='lig2rec') %>%
   dplyr::select(from,to,edgetype) %>% distinct()

}
) %>% bind_rows()

nodes$receptor <- pmap(edges$lig2rec,function(to,...){
  data.frame(name=to,type='receptor',size=8,color=pal[4])
}) %>% bind_rows() %>% distinct()



=======


#Perhaps add filter for Source and Evidence.

edges$lig2rec <-
  pmap(nodes$lig,function(name,...){
 name <- enquo(name)
 recligDB %>%  dplyr::filter(ligand %in% !!name & receptor %in% ExpressGenes$receptor.cell) %>%
   mutate(from=ligand,to=receptor,edgetype='lig2rec') %>%
   dplyr::select(from,to,edgetype) %>% distinct()

}
) %>% bind_rows()

nodes$receptor <- pmap(edges$lig2rec,function(to,...){
  data.frame(name=to,type='receptor',size=8,color=pal[3])
}) %>% bind_rows() %>% distinct()



>>>>>>> de93b641d163a4fb71c34dbe9fbfffd16f487b28
# Rec to cell -------------------------------------------------------------
edges$rec2cell <-
  pmap(nodes$receptor,function(name,...){
    data.frame(from=name,to=as.character(receptor.cell),edgetype='rec2cell')
  }
  ) %>% bind_rows()
<<<<<<< HEAD





=======





>>>>>>> de93b641d163a4fb71c34dbe9fbfffd16f487b28
####################################### OLD CODE #########################

  # ### Add receptors  to cells
  # nodes <- rbind(nodes, mm %>% dplyr::filter(receptor %in% genesExprss) %>% mutate(name=receptor, type='receptor',size=8,color=pal[4]) %>% dplyr::select(name,type,size,color))
  # edges <- rbind(edges,mm %>% dplyr::filter(receptor %in% genesExprss) %>% mutate(from=receptor, to=motifRes[[n]]$celltype, edgetype='rec2cell') %>% dplyr::select(from,to,edgetype))
  #
  # rec2tf <- rectf %>% dplyr::filter(receptor %in% genesExprss & tf %in% genesExprss) %>%
  #   mutate(from=receptor, to=tf, edgetype='rec2tf')  %>%
  #   dplyr::select(from,to,edgetype) %>%
  #   distinct()
  #
  # for(tf in rec2tf$to){
  #   m <- motifRes[[n]]$motifs %>% dplyr::filter(motif.name==tf) %>% pull(motif)
  #   if(length(m)==0){next}
  #   t1 <- motif[,m]==1
  #   if(length(m)==1){
  #     geneid <- motifRes[[n]]$peaks.ann %>% dplyr::filter(peakid %in% names(t1)[which(t1==TRUE)]) %>% pull(gene_id)
  #   }else{
  #     geneid <- motifRes[[n]]$peaks.ann %>% dplyr::filter(peakid %in% rownames(t1)[(Matrix::rowSums(t1) > 0L)==TRUE]) %>% pull(gene_id)
  #   }
  #
  #   #nodes <- rbind(nodes, egoBP@result %>% mutate(name=Description, type='GO',size=8,color=pal[4]) %>% dplyr::select(name,type,size,color))
  #   #edges <- rbind(edges,egoBP@result %>% mutate(from=tf, to=Description, edgetype='tf2GO') %>% dplyr::select(from,to,edgetype))
  #
  # }
  #
  #
  # nodes <- rbind(nodes, rec2tf %>% mutate(type='receptor',size=8,color=pal[4]) %>% dplyr::rename(name=from) %>% dplyr::select(name,type,size,color))
  # nodes <- rbind(nodes, rec2tf %>% mutate(type='TF',size=8,color=pal[2]) %>% dplyr::rename(name=to) %>% dplyr::select(name,type,size,color))
  # edges <- rbind(edges,rec2tf)
  # edges <- rbind(edges,rec2tf %>% mutate(from=to, to=motifRes[[n]]$celltype, edgetype='tf2cell'))


######################################################################

E <- bind_rows(edges) %>% distinct()
N <- bind_rows(nodes) %>% distinct()


net <- igraph::graph_from_data_frame(d=E, vertices=N, directed=T)
net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)

cytoscapePing()
createNetworkFromIgraph(net,"myIgraph")

style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')

#and then create the style
createVisualStyle(style.name, defaults, list(nodeLabels))
<<<<<<< HEAD

#finsh by applying the style
setVisualStyle(style.name)

=======

#finsh by applying the style
setVisualStyle(style.name)

>>>>>>> de93b641d163a4fb71c34dbe9fbfffd16f487b28
setNodeColorMapping('color',mapping.type = 'p')

getLayoutNames()
layoutNetwork('force-directed defaultSpringLength=70 defaultSpringCoefficient=0.000003')
<<<<<<< HEAD

layoutNetwork('hierarchical')
=======
>>>>>>> de93b641d163a4fb71c34dbe9fbfffd16f487b28


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



