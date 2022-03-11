#!/app/easybuild/software/R/3.5.1-foss-2016b-fh1/bin/Rscript

#Jenny Smith 
#8/28/18

library(methods)

setwd('/fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2017.10.25_1031_NUP98-NSD1/')
source("~/scripts/RNAseq_Analysis/DifferentialExpn_PathwayAnalysis/GAGE_GSEA_Function.r")

#For most recent analysis (49 NSD1 fusion pos patients)
# NSD1 <- get(load("TARGET_AML_0531_1031_NUP98.NSD1_vs_otherAML_IncludingOtherNUPs_DEGs.RDS"))

#For Sanne's Analysis (47 NSD1 fusion pos patients)
# NSD1 <- get(load("DEGs/Voom_ForSanne/TARGET_AML_0531_1031_NUP98.NSD1_vs_otherAML_DEGs_list.RDS"))
# NSD1 <- NSD1[sapply(NSD1, function(x) length(x) > 3)]


#For FLT3-ITD 
NSD1.ITDvsNSD1 <- readRDS("TARGET_AML_1031_NUP98.NSD1.FLT3-ITD_vs_NSD1.alone_DEGs.RDS")
NSD1.ITDvsFLT3 <- readRDS("TARGET_AML_1031_NUP98.NSD1.FLT3-ITD_vs_FLT3.alone_DEGs.RDS")
NSD1 <- list("NSD1.ITDvsNSD1"=NSD1.ITDvsNSD1,"NSD1.ITDvsFLT3"=NSD1.ITDvsFLT3)


filename <- "TARGET_AML_1031_NSD1.ITDvsFLT3_NSD1.ITDvsNSD1_GAGE"
C2.KEGG <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.cp.kegg.v6.0.symbols.RDS")



print("starting1")

GSA <- lapply(NSD1, gage_from_pipeline,
              method="voom",
              type="expn",
              geneset=C2.KEGG)
  saveRDS(GSA,file=paste0(filename, "_expn_C2.KEGG.RDS"))
  rm(GSA)
  gc()

print("done1")


print("starting2")


C2.All <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c2.all.v6.0.symbols.RDS")

GSA.C2.All <- lapply(NSD1, gage_from_pipeline, 
                     method="voom",
                     type="expn",
                     geneset=C2.All)
saveRDS(GSA.C2.All,file=paste0(filename, "_expn_C2.All.RDS"))
rm(GSA.C2.All)
gc()

print("done2")


print("starting3")

GSA.KEGG <- lapply(NSD1, gage_from_pipeline,
                   method="voom",
                   type="expn",
                   geneset=NULL)
  saveRDS(GSA.KEGG,file=paste0(filename, "_expn_HSA.KEGG.RDS"))
  rm(GSA.KEGG)
  gc()

print("done3")

C5 <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c5_list_SetSize_GE.50_LE.300_v6.1.symbols.RDS")

print("starting4")
GSA.GO.BioProcess <- lapply(NSD1, gage_from_pipeline,
                            method="voom",
                            type="expn", 
                            geneset=C5[["c5.bp"]])
saveRDS(GSA.GO.BioProcess, file=paste0(filename, "_expn_C5.BioProcess_SetSize50to300.RDS"))
rm(GSA.GO.BioProcess)
gc()
print("done4")



print("starting5")
GSA.GO.CellComp <- lapply(NSD1, gage_from_pipeline, 
                          method="voom",
                          type="expn", 
                          geneset=C5[["c5.cc"]])
saveRDS(GSA.GO.CellComp, file=paste0(filename, "_expn_C5.CellComp_SetSize50to300.RDS"))
rm(GSA.GO.CellComp)
gc()
print("done5")    



print("starting6")
GSA.GO.MolFunc <- lapply(NSD1, gage_from_pipeline, 
                         method="voom",
                         type="expn",
                         geneset=C5[["c5.mf"]])
saveRDS(GSA.GO.MolFunc, file=paste0(filename, "_expn_C5.MolFunc_SetSize50to300.RDS"))
rm(GSA.GO.MolFunc)
gc()
print("done6")




# print("starting7")
# 
# C3 <- readRDS("~/RNA_seq_Analysis/0000.00.01_GSEA_geneSets_gmt/c3.tft.v6.2.symbols.RDS")
# 
# print("starting11")
# GSA.TFB <- lapply(NSD1, gage_from_pipeline,
#                             method="voom",
#                             type="expn", 
#                             geneset=C3)
# saveRDS(GSA.TFB, file=paste0(filename, "_expn_C3_TFbindingSites.RDS"))
# rm(GSA.TFB)
# gc()
# print("done7")







