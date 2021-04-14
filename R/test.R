

rna.object <- readRDS("~/dsdata/projects/Morrisey/Jarod/scRNA/JZ_lung_timeseries/Adult/seuratv4/Seurat.RDS")
atac.object <- readRDS("~/dsdata/projects/Morrisey/JohnLeach/scATAC/Adult/Signac/Signac.RDS")
ligand.cell=10 #Amp
receptor.cell=8 #MANC


recligDB <- readr::read_csv('~/dsdata/NGSshare/FANTOM5/Mm_PairsLigRec.csv')
tf <- readr::read_tsv('~/dsdata/NGSshare/AnimalTFDB/Mus_musculus_TF_v3.txt')
mm2hs <- readr::read_csv('~/dsdata/NGSshare/homologs/mouse_human.csv')



runSCAPE()
