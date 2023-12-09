# CNPS.cycle
This is an R package for element cycle analysis using metagenomic data.
![GF](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/183e531f-31ff-4bb0-9504-0635b67422a7)
# Citation
Zhengfu Yue, Jing Zhang, Rui Kou, Tianshun Liu, Ye Tao, Liang Zeng, Zelong Zhao, Shoushan Sheng, Qinfen Li, Jing Zhang, Yukun Zou.2023.CNPS.cycle: An R package for element cycle gene and functional microbial analysis based on shotgun metagenomic data
# Install
CNPS.cycle is available on Github, you can install it by:
```{r}
library(devtools) 
install_github("yuezhengfu/CNPS.cycle")
```
# Usage
Carbon cycle analysis
```{r}
#KO number associated with C cycle
C <- c("K00855","K01061","K01602","K08684","K02256","K02262","K02274","K02276",
       "K00174","K00175","K00244","K01648","K00194","K00197","K03518","K03519",
       "K03520","K00016","K00400","K00401")
#The KO abundance information table related to C cycle
C.ko <- ko[rownames(ko) %in% C,]
write.table(C.ko,file = "Results/Carbon/Gene/Abundance/C_cycle_ko_abun.txt",sep = "\t")
#The KO abundance information table was combined according to 7 carbon cycle processes
C.abundance <- Ccyc.abundance(ko)
write.table(C.abundance,file = "Results/Carbon/Gene/Abundance/C_cycle_gene_abun.txt",
            sep = "\t",row.names = FALSE)
#The carbon cycle abundance heatmap and the difference test were drawn according to the sample groups
dir.create("Results/Carbon/Gene/Heatmap")
result <- abun.heatmap.g(C.abundance,group,Group_numb)
pdf(file = "Results/Carbon/Gene/Heatmap/C_cycle_gene_abun_group.pdf",
    width = 2.3 + 0.6*Group_numb,height = 3)
result[[2]]
dev.off()
write.table(result[[1]],"Results/Carbon/Gene/Heatmap/Diff_C_gene_test.txt",sep = "\t",
            row.names = FALSE)
#Carbon cycle differential gene analysis and visualization
dir.create("Results/Carbon/Gene/Cycle image")
result <- fold.change(C.abundance,group)
write.table(result[[1]],"Results/Carbon/Gene/Cycle image/Gene_fold_change.txt",sep = "\t",
            row.names = FALSE)
pdf(file = "Results/Carbon/Gene/Cycle image/Gene_fold_change.pdf",
    width = 0.6*Group_numb + 2.3,height = 4.5)
result[[3]]
dev.off()
#Heatmap subsets of 7 carbon cycle categories were extracted
ACF <- result[[2]] + ylim("Aerobic C fixation")
ACH4O <- result[[2]] + ylim("Aerobic CH4 oxidation")
AR <- result[[2]] + ylim("Aerobic respiration")
AnCF <- result[[2]] + ylim("Anaerobic C fixation")
COo <- result[[2]] + ylim("CO oxidation")
Fer <- result[[2]] + ylim("Fermentation")
Meth <- result[[2]] + ylim("Methanogenesis")
#carbon cycle patterns
C.img <- system.file("data", "Ccyc.pdf", package = "CNPS.cycle")
pdf("Results/Carbon/Gene/Cycle image/C_cyc_fold_change.pdf",width = 13,height = 7)
gg <- ggplot()
ggbackground(gg,C.img)
if (rowSums(C.abundance[1,2:ncol(C.abundance)]) > 0) {
  print(ACF,vp = viewport(width = 0.035*Group_numb,height = 0.05,
                          x = 0.692,y = 0.87))
}
if (rowSums(C.abundance[2,2:ncol(C.abundance)]) > 0) {
  print(ACH4O,vp = viewport(width = 0.035*Group_numb,height = 0.05,
                            x = 0.5935,y = 0.23))
}
if (rowSums(C.abundance[3,2:ncol(C.abundance)]) > 0) {
  print(AR,vp = viewport(width = 0.035*Group_numb,height = 0.05,
                         x = 0.612,y = 0.64))
}
if (rowSums(C.abundance[4,2:ncol(C.abundance)]) > 0) {
  print(AnCF,vp = viewport(width = 0.035*Group_numb,height = 0.05,
                           x = 0.322,y = 0.86))
}
if (rowSums(C.abundance[5,2:ncol(C.abundance)]) > 0) {
  print(COo,vp = viewport(width = 0.035*Group_numb,height = 0.05,
                          x = 0.812,y = 0.385))
}
if (rowSums(C.abundance[6,2:ncol(C.abundance)]) > 0) {
  print(Fer,vp = viewport(width = 0.035*Group_numb,height = 0.05,
                          x = 0.36,y = 0.45))
}
if (rowSums(C.abundance[7,2:ncol(C.abundance)]) > 0) {
  print(Meth,vp = viewport(width = 0.035*Group_numb,height = 0.05,
                           x = 0.362,y = 0.34))
}
dev.off()

#Relative abundance statistics of microorganisms host to carbon cycling genes
#Aerobic C fixation宿主统计
if (rowSums(C.abundance[1,2:ncol(C.abundance)]) > 0 & Ccyc.h[1] == 1) {
  dir.create("Results/Carbon/Host_relative_Group/ACF")
  ACF <- ACF(Gene,tax,abundance,group)
  write.csv(ACF[[1]],file = "Results/Carbon/Host_relative_Group/ACF/ACF_phylum.csv")
  write.csv(ACF[[2]],file = "Results/Carbon/Host_relative_Group/ACF/ACF_class.csv")
  write.csv(ACF[[3]],file = "Results/Carbon/Host_relative_Group/ACF/ACF_order.csv")
  write.csv(ACF[[4]],file = "Results/Carbon/Host_relative_Group/ACF/ACF_family.csv")
  write.csv(ACF[[5]],file = "Results/Carbon/Host_relative_Group/ACF/ACF_genus.csv")
  write.csv(ACF[[6]],file = "Results/Carbon/Host_relative_Group/ACF/ACF_species.csv")
}

#Aerobic CH4 oxidation宿主统计
if (rowSums(C.abundance[2,2:ncol(C.abundance)]) > 0 & Ccyc.h[2] == 1) {
  dir.create("Results/Carbon/Host_relative_Group/ACH4O")
  ACH4O <- ACH4O(Gene,tax,abundance,group)
  write.csv(ACH4O[[1]],file = "Results/Carbon/Host_relative_Group/ACH4O/ACH4O_phylum.csv")
  write.csv(ACH4O[[2]],file = "Results/Carbon/Host_relative_Group/ACH4O/ACH4O_class.csv")
  write.csv(ACH4O[[3]],file = "Results/Carbon/Host_relative_Group/ACH4O/ACH4O_order.csv")
  write.csv(ACH4O[[4]],file = "Results/Carbon/Host_relative_Group/ACH4O/ACH4O_family.csv")
  write.csv(ACH4O[[5]],file = "Results/Carbon/Host_relative_Group/ACH4O/ACH4O_genus.csv")
  write.csv(ACH4O[[6]],file = "Results/Carbon/Host_relative_Group/ACH4O/ACH4O_species.csv")
}

#Aerobic respiration宿主统计
if (rowSums(C.abundance[3,2:ncol(C.abundance)]) > 0 & Ccyc.h[3] == 1) {
  dir.create("Results/Carbon/Host_relative_Group/AR")
  AR <- AR(Gene,tax,abundance,group)
  write.csv(AR[[1]],file = "Results/Carbon/Host_relative_Group/AR/AR_phylum.csv")
  write.csv(AR[[2]],file = "Results/Carbon/Host_relative_Group/AR/AR_class.csv")
  write.csv(AR[[3]],file = "Results/Carbon/Host_relative_Group/AR/AR_order.csv")
  write.csv(AR[[4]],file = "Results/Carbon/Host_relative_Group/AR/AR_family.csv")
  write.csv(AR[[5]],file = "Results/Carbon/Host_relative_Group/AR/AR_genus.csv")
  write.csv(AR[[6]],file = "Results/Carbon/Host_relative_Group/AR/AR_species.csv")
}

#Anaerobic C fixation宿主统计
if (rowSums(C.abundance[4,2:ncol(C.abundance)]) > 0 & Ccyc.h[4] == 1) {
  dir.create("Results/Carbon/Host_relative_Group/AnCF")
  AnCF <- AnCF(Gene,tax,abundance,group)
  write.csv(AnCF[[1]],file = "Results/Carbon/Host_relative_Group/AnCF/AnCF_phylum.csv")
  write.csv(AnCF[[2]],file = "Results/Carbon/Host_relative_Group/AnCF/AnCF_class.csv")
  write.csv(AnCF[[3]],file = "Results/Carbon/Host_relative_Group/AnCF/AnCF_order.csv")
  write.csv(AnCF[[4]],file = "Results/Carbon/Host_relative_Group/AnCF/AnCF_family.csv")
  write.csv(AnCF[[5]],file = "Results/Carbon/Host_relative_Group/AnCF/AnCF_genus.csv")
  write.csv(AnCF[[6]],file = "Results/Carbon/Host_relative_Group/AnCF/AnCF_species.csv")
}

#CO oxidation宿主统计
if (rowSums(C.abundance[5,2:ncol(C.abundance)]) > 0 & Ccyc.h[5] == 1) {
  dir.create("Results/Carbon/Host_relative_Group/COo")
  COo <- COo(Gene,tax,abundance,group)
  write.csv(COo[[1]],file = "Results/Carbon/Host_relative_Group/COo/COo_phylum.csv")
  write.csv(COo[[2]],file = "Results/Carbon/Host_relative_Group/COo/COo_class.csv")
  write.csv(COo[[3]],file = "Results/Carbon/Host_relative_Group/COo/COo_order.csv")
  write.csv(COo[[4]],file = "Results/Carbon/Host_relative_Group/COo/COo_family.csv")
  write.csv(COo[[5]],file = "Results/Carbon/Host_relative_Group/COo/COo_genus.csv")
  write.csv(COo[[6]],file = "Results/Carbon/Host_relative_Group/COo/COo_species.csv")
}

#Fermentation宿主统计
if (rowSums(C.abundance[6,2:ncol(C.abundance)]) > 0 & Ccyc.h[6] == 1) {
  dir.create("Results/Carbon/Host_relative_Group/Fer")
  Fer <- Fer(Gene,tax,abundance,group)
  write.csv(Fer[[1]],file = "Results/Carbon/Host_relative_Group/Fer/Fer_phylum.csv")
  write.csv(Fer[[2]],file = "Results/Carbon/Host_relative_Group/Fer/Fer_class.csv")
  write.csv(Fer[[3]],file = "Results/Carbon/Host_relative_Group/Fer/Fer_order.csv")
  write.csv(Fer[[4]],file = "Results/Carbon/Host_relative_Group/Fer/Fer_family.csv")
  write.csv(Fer[[5]],file = "Results/Carbon/Host_relative_Group/Fer/Fer_genus.csv")
  write.csv(Fer[[6]],file = "Results/Carbon/Host_relative_Group/Fer/Fer_species.csv")
}

#Methanogenesis宿主统计
#if (rowSums(C.abundance[7,2:ncol(C.abundance)]) > 0 & Ccyc.h[7] == 1) {
#  dir.create("Results/Carbon/Host_relative_Group/Meth")
#  Meth <- Meth(Gene,tax,abundance,group)
#  write.csv(Meth[[1]],file = "Results/Carbon/Host_relative_Group/Meth/Meth_phylum.csv")
#  write.csv(Meth[[2]],file = "Results/Carbon/Host_relative_Group/Meth/Meth_class.csv")
#  write.csv(Meth[[3]],file = "Results/Carbon/Host_relative_Group/Meth/Meth_order.csv")
#  write.csv(Meth[[4]],file = "Results/Carbon/Host_relative_Group/Meth/Meth_family.csv")
#  write.csv(Meth[[5]],file = "Results/Carbon/Host_relative_Group/Meth/Meth_genus.csv")
#  write.csv(Meth[[6]],file = "Results/Carbon/Host_relative_Group/Meth/Meth_species.csv")
#}

#Visualization of carbon cycle host statistics
#Aerobic C fixation
if (rowSums(C.abundance[1,2:ncol(C.abundance)]) > 0 & Ccyc.h[1] == 1) {
  title <- "Aerobic C fixation"
  aa <- ifelse(nrow(ACF[[1]]) > 5,
               max(str_length(rownames(ACF[[1]])[1:5])),
               max(str_length(rownames(ACF[[1]]))))
  pdf("Results/Carbon/Host_relative_Group/ACF/ACF_phylum.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACF[[1]],title)
  dev.off()
  aa <- ifelse(nrow(ACF[[2]]) > 5,
               max(str_length(rownames(ACF[[2]])[1:5])),
               max(str_length(rownames(ACF[[2]]))))
  pdf("Results/Carbon/Host_relative_Group/ACF/ACF_class.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACF[[2]],title)
  dev.off()
  aa <- ifelse(nrow(ACF[[3]]) > 5,
               max(str_length(rownames(ACF[[3]])[1:5])),
               max(str_length(rownames(ACF[[3]]))))
  pdf("Results/Carbon/Host_relative_Group/ACF/ACF_order.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACF[[3]],title)
  dev.off()
  aa <- ifelse(nrow(ACF[[4]]) > 5,
               max(str_length(rownames(ACF[[4]])[1:5])),
               max(str_length(rownames(ACF[[4]]))))
  pdf("Results/Carbon/Host_relative_Group/ACF/ACF_family.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACF[[4]],title)
  dev.off()
  aa <- ifelse(nrow(ACF[[5]]) > 5,
               max(str_length(rownames(ACF[[5]])[1:5])),
               max(str_length(rownames(ACF[[5]]))))
  pdf("Results/Carbon/Host_relative_Group/ACF/ACF_genus.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACF[[5]],title)
  dev.off()
  aa <- ifelse(nrow(ACF[[6]]) > 5,
               max(str_length(rownames(ACF[[6]])[1:5])),
               max(str_length(rownames(ACF[[6]]))))
  pdf("Results/Carbon/Host_relative_Group/ACF/ACF_species.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACF[[6]],title)
  dev.off()
}

#Aerobic CH4 oxidation
if (rowSums(C.abundance[2,2:ncol(C.abundance)]) > 0 & Ccyc.h[2] == 1) {
  title <- expression(paste("Aerobic CH"["4"]," oxidation"))
  aa <- ifelse(nrow(ACH4O[[1]]) > 5,
               max(str_length(rownames(ACH4O[[1]])[1:5])),
               max(str_length(rownames(ACH4O[[1]]))))
  pdf("Results/Carbon/Host_relative_Group/ACH4O/ACH4O_phylum.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACH4O[[1]],title)
  dev.off()
  aa <- ifelse(nrow(ACH4O[[2]]) > 5,
               max(str_length(rownames(ACH4O[[2]])[1:5])),
               max(str_length(rownames(ACH4O[[2]]))))
  pdf("Results/Carbon/Host_relative_Group/ACH4O/ACH4O_class.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACH4O[[2]],title)
  dev.off()
  aa <- ifelse(nrow(ACH4O[[3]]) > 5,
               max(str_length(rownames(ACH4O[[3]])[1:5])),
               max(str_length(rownames(ACH4O[[3]]))))
  pdf("Results/Carbon/Host_relative_Group/ACH4O/ACH4O_order.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACH4O[[3]],title)
  dev.off()
  aa <- ifelse(nrow(ACH4O[[4]]) > 5,
               max(str_length(rownames(ACH4O[[4]])[1:5])),
               max(str_length(rownames(ACH4O[[4]]))))
  pdf("Results/Carbon/Host_relative_Group/ACH4O/ACH4O_family.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACH4O[[4]],title)
  dev.off()
  aa <- ifelse(nrow(ACH4O[[5]]) > 5,
               max(str_length(rownames(ACH4O[[5]])[1:5])),
               max(str_length(rownames(ACH4O[[5]]))))
  pdf("Results/Carbon/Host_relative_Group/ACH4O/ACH4O_genus.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACH4O[[5]],title)
  dev.off()
  aa <- ifelse(nrow(ACH4O[[6]]) > 5,
               max(str_length(rownames(ACH4O[[6]])[1:5])),
               max(str_length(rownames(ACH4O[[6]]))))
  pdf("Results/Carbon/Host_relative_Group/ACH4O/ACH4O_species.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(ACH4O[[6]],title)
  dev.off()
}

#Anaerobic C fixation
if (rowSums(C.abundance[4,2:ncol(C.abundance)]) > 0 & Ccyc.h[3] == 1) {
  title <- "Anaerobic C fixation"
  aa <- ifelse(nrow(AnCF[[1]]) > 5,
               max(str_length(rownames(AnCF[[1]])[1:5])),
               max(str_length(rownames(AnCF[[1]]))))
  pdf("Results/Carbon/Host_relative_Group/AnCF/AnCF_phylum.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AnCF[[1]],title)
  dev.off()
  aa <- ifelse(nrow(AnCF[[2]]) > 5,
               max(str_length(rownames(AnCF[[2]])[1:5])),
               max(str_length(rownames(AnCF[[2]]))))
  pdf("Results/Carbon/Host_relative_Group/AnCF/AnCF_class.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AnCF[[2]],title)
  dev.off()
  aa <- ifelse(nrow(AnCF[[3]]) > 5,
               max(str_length(rownames(AnCF[[3]])[1:5])),
               max(str_length(rownames(AnCF[[3]]))))
  pdf("Results/Carbon/Host_relative_Group/AnCF/AnCF_order.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AnCF[[3]],title)
  dev.off()
  aa <- ifelse(nrow(AnCF[[4]]) > 5,
               max(str_length(rownames(AnCF[[4]])[1:5])),
               max(str_length(rownames(AnCF[[4]]))))
  pdf("Results/Carbon/Host_relative_Group/AnCF/AnCF_family.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AnCF[[4]],title)
  dev.off()
  aa <- ifelse(nrow(AnCF[[5]]) > 5,
               max(str_length(rownames(AnCF[[5]])[1:5])),
               max(str_length(rownames(AnCF[[5]]))))
  pdf("Results/Carbon/Host_relative_Group/AnCF/AnCF_genus.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AnCF[[5]],title)
  dev.off()
  aa <- ifelse(nrow(AnCF[[6]]) > 5,
               max(str_length(rownames(AnCF[[6]])[1:5])),
               max(str_length(rownames(AnCF[[6]]))))
  pdf("Results/Carbon/Host_relative_Group/AnCF/AnCF_species.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AnCF[[6]],title)
  dev.off()
}

#Aerobic respiration
if (rowSums(C.abundance[3,2:ncol(C.abundance)]) > 0 & Ccyc.h[4] == 1) {
  title <- "Aerobic respiration"
  aa <- ifelse(nrow(AR[[1]]) > 5,
               max(str_length(rownames(AR[[1]])[1:5])),
               max(str_length(rownames(AR[[1]]))))
  pdf("Results/Carbon/Host_relative_Group/AR/AR_phylum.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AR[[1]],title)
  dev.off()
  aa <- ifelse(nrow(AR[[2]]) > 5,
               max(str_length(rownames(AR[[2]])[1:5])),
               max(str_length(rownames(AR[[2]]))))
  pdf("Results/Carbon/Host_relative_Group/AR/AR_class.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AR[[2]],title)
  dev.off()
  aa <- ifelse(nrow(AR[[3]]) > 5,
               max(str_length(rownames(AR[[3]])[1:5])),
               max(str_length(rownames(AR[[3]]))))
  pdf("Results/Carbon/Host_relative_Group/AR/AR_order.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AR[[3]],title)
  dev.off()
  aa <- ifelse(nrow(AR[[4]]) > 5,
               max(str_length(rownames(AR[[4]])[1:5])),
               max(str_length(rownames(AR[[4]]))))
  pdf("Results/Carbon/Host_relative_Group/AR/AR_family.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AR[[4]],title)
  dev.off()
  aa <- ifelse(nrow(AR[[5]]) > 5,
               max(str_length(rownames(AR[[5]])[1:5])),
               max(str_length(rownames(AR[[5]]))))
  pdf("Results/Carbon/Host_relative_Group/AR/AR_genus.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AR[[5]],title)
  dev.off()
  aa <- ifelse(nrow(AR[[6]]) > 5,
               max(str_length(rownames(AR[[6]])[1:5])),
               max(str_length(rownames(AR[[6]]))))
  pdf("Results/Carbon/Host_relative_Group/AR/AR_species.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(AR[[6]],title)
  dev.off()
}

#CO oxidation
if (rowSums(C.abundance[5,2:ncol(C.abundance)]) > 0 & Ccyc.h[5] == 1) {
  title <- "CO oxidation"
  aa <- ifelse(nrow(COo[[1]]) > 5,
               max(str_length(rownames(COo[[1]])[1:5])),
               max(str_length(rownames(COo[[1]]))))
  pdf("Results/Carbon/Host_relative_Group/COo/COo_phylum.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(COo[[1]],title)
  dev.off()
  aa <- ifelse(nrow(COo[[2]]) > 5,
               max(str_length(rownames(COo[[2]])[1:5])),
               max(str_length(rownames(COo[[2]]))))
  pdf("Results/Carbon/Host_relative_Group/COo/COo_class.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(COo[[2]],title)
  dev.off()
  aa <- ifelse(nrow(COo[[3]]) > 5,
               max(str_length(rownames(COo[[3]])[1:5])),
               max(str_length(rownames(COo[[3]]))))
  pdf("Results/Carbon/Host_relative_Group/COo/COo_order.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(COo[[3]],title)
  dev.off()
  aa <- ifelse(nrow(COo[[4]]) > 5,
               max(str_length(rownames(COo[[4]])[1:5])),
               max(str_length(rownames(COo[[4]]))))
  pdf("Results/Carbon/Host_relative_Group/COo/COo_family.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(COo[[4]],title)
  dev.off()
  aa <- ifelse(nrow(COo[[5]]) > 5,
               max(str_length(rownames(COo[[5]])[1:5])),
               max(str_length(rownames(COo[[5]]))))
  pdf("Results/Carbon/Host_relative_Group/COo/COo_genus.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(COo[[5]],title)
  dev.off()
  aa <- ifelse(nrow(COo[[6]]) > 5,
               max(str_length(rownames(COo[[6]])[1:5])),
               max(str_length(rownames(COo[[6]]))))
  pdf("Results/Carbon/Host_relative_Group/COo/COo_species.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(COo[[6]],title)
  dev.off()
}

#Fermentation
if (rowSums(C.abundance[6,2:ncol(C.abundance)]) > 0 & Ccyc.h[6] == 1) {
  title <- "Fermentation"
  aa <- ifelse(nrow(Fer[[1]]) > 5,
               max(str_length(rownames(Fer[[1]])[1:5])),
               max(str_length(rownames(Fer[[1]]))))
  pdf("Results/Carbon/Host_relative_Group/Fer/Fer_phylum.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(Fer[[1]],title)
  dev.off()
  aa <- ifelse(nrow(Fer[[2]]) > 5,
               max(str_length(rownames(Fer[[2]])[1:5])),
               max(str_length(rownames(Fer[[2]]))))
  pdf("Results/Carbon/Host_relative_Group/Fer/Fer_class.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(Fer[[2]],title)
  dev.off()
  aa <- ifelse(nrow(Fer[[3]]) > 5,
               max(str_length(rownames(Fer[[3]])[1:5])),
               max(str_length(rownames(Fer[[3]]))))
  pdf("Results/Carbon/Host_relative_Group/Fer/Fer_order.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(Fer[[3]],title)
  dev.off()
  aa <- ifelse(nrow(Fer[[4]]) > 5,
               max(str_length(rownames(Fer[[4]])[1:5])),
               max(str_length(rownames(Fer[[4]]))))
  pdf("Results/Carbon/Host_relative_Group/Fer/Fer_family.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(Fer[[4]],title)
  dev.off()
  aa <- ifelse(nrow(Fer[[5]]) > 5,
               max(str_length(rownames(Fer[[5]])[1:5])),
               max(str_length(rownames(Fer[[5]]))))
  pdf("Results/Carbon/Host_relative_Group/Fer/Fer_genus.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(Fer[[5]],title)
  dev.off()
  aa <- ifelse(nrow(Fer[[6]]) > 5,
               max(str_length(rownames(Fer[[6]])[1:5])),
               max(str_length(rownames(Fer[[6]]))))
  pdf("Results/Carbon/Host_relative_Group/Fer/Fer_species.pdf",
      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
  host.ratio(Fer[[6]],title)
  dev.off()
}

#Methanogenesis
#if (rowSums(C.abundance[7,2:ncol(C.abundance)]) > 0 & Ccyc.h[7] == 1) {
#  title <- "Methanogenesis"
#  aa <- ifelse(nrow(Meth[[1]]) > 5,
#             max(str_length(rownames(Meth[[1]])[1:5])),
#               max(str_length(rownames(Meth[[1]]))))
#  pdf("Results/Carbon/Host_relative_Group/Meth/Meth_phylum.pdf",
#      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
#  host.ratio(Meth[[1]],title)
#  dev.off()
#  aa <- ifelse(nrow(Meth[[2]]) > 5,
#             max(str_length(rownames(Meth[[2]])[1:5])),
#               max(str_length(rownames(Meth[[2]]))))
#  pdf("Results/Carbon/Host_relative_Group/Meth/Meth_class.pdf",
#      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
#  host.ratio(Meth[[2]],title)
#  dev.off()
#  aa <- ifelse(nrow(Meth[[3]]) > 5,
#             max(str_length(rownames(Meth[[3]])[1:5])),
#               max(str_length(rownames(Meth[[3]]))))
#  pdf("Results/Carbon/Host_relative_Group/Meth/Meth_order.pdf",
#      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
#  host.ratio(Meth[[3]],title)
#  dev.off()
#  aa <- ifelse(nrow(Meth[[4]]) > 5,
#             max(str_length(rownames(Meth[[4]])[1:5])),
#               max(str_length(rownames(Meth[[4]]))))
#  pdf("Results/Carbon/Host_relative_Group/Meth/Meth_family.pdf",
#      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
#  host.ratio(Meth[[4]],title)
#  dev.off()
#  aa <- ifelse(nrow(Meth[[5]]) > 5,
#             max(str_length(rownames(Meth[[5]])[1:5])),
#               max(str_length(rownames(Meth[[5]]))))
#  pdf("Results/Carbon/Host_relative_Group/Meth/Meth_genus.pdf",
#      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
#  host.ratio(Meth[[5]],title)
#  dev.off()
#  aa <- ifelse(nrow(Meth[[6]]) > 5,
#             max(str_length(rownames(Meth[[6]])[1:5])),
#               max(str_length(rownames(Meth[[6]]))))
#  pdf("Results/Carbon/Host_relative_Group/Meth/Meth_species.pdf",
#      width = ifelse(aa > 40,aa*0.1 + 0.5*Group_numb,3.2 + 0.5*Group_numb),height = 3.5)
#  host.ratio(Meth[[6]],title)
#  dev.off()
#}

#### Legend
df <- data.frame(A = c(0,0.1,0.2,0.3,0.4,0.5),
                 B = c(0.5,0.6,0.7,0.8,0.9,1))
df <- t(df)
df1 <- data.frame(A = c("0%","10%","20%","30%","40%","50%"),
                  B = c("50%","60%","70%","80%","90%","100%"))
df1 <- t(df1)
bk <- c(seq(0,1,by = 0.01))
pheatmap(df,fontsize = 30,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
         cellwidth = 80,cellheight = 50,legend = FALSE,breaks = bk,show_rownames = FALSE,
         color = colorRampPalette(c("white","Red"))(100),show_colnames = FALSE,
         display_numbers = df1,number_color = "black",border_color = "black",
         filename = "Results/Carbon/Host_relative_Group/Legend_relative_Group.pdf",width = 7,height = 2)

#Analysis and visualization of carbon cycle β-diversity
cbbPalette <- c("#B2182B","#56B4E9","#E69F00","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",
                "#ADD1E5")

dir.create("Results/Carbon/Beta diversity")
dir.create("Results/Carbon/Beta diversity/Distance")
dir.create("Results/Carbon/Beta diversity/PCA")
dir.create("Results/Carbon/Beta diversity/PCoA")
dir.create("Results/Carbon/Beta diversity/NMDS")

result <- pcoa.arg(C.abundance[,2:ncol(C.abundance)],group)
write.table(as.matrix(result[[1]]),"Results/Carbon/Beta diversity/Distance/distance_bray_curtis.txt",
            sep = "\t")
write.table(result[[2]],"Results/Carbon/Beta diversity/Distance/diff_test.txt",sep = "\t",row.names = FALSE)
write.table(result[[3]],"Results/Carbon/Beta diversity/PCoA/C_pcoa.txt",sep = "\t",row.names = FALSE)
pdf(file = "Results/Carbon/Beta diversity/PCoA/C_pcoa_group.pdf",width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Carbon/Beta diversity/PCoA/C_pcoa_ellipse.pdf",width = 7.5,height = 5.4)
result[[5]]
dev.off()
pdf(file = "Results/Carbon/Beta diversity/PCoA/C_pcoa_label.pdf",width = 7.5,height = 5.4)
result[[6]]
dev.off()

result <- pca.arg(C.abundance[,2:ncol(C.abundance)],group)
write.table(result[[1]],"Results/Carbon/Beta diversity/PCA/C_pca.txt",sep = "\t",row.names = FALSE)
pdf(file = "Results/Carbon/Beta diversity/PCA/C_pca_group.pdf",width = 7.5,height = 5.4)
result[[2]]
dev.off()
pdf(file = "Results/Carbon/Beta diversity/PCA/C_pca_ellipse.pdf",width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Carbon/Beta diversity/PCA/C_pca_label.pdf",width = 7.5,height = 5.4)
result[[4]]
dev.off()

result <- nmds.arg(C.abundance[,2:ncol(C.abundance)],group)
write.table(result[[1]],"Results/Carbon/Beta diversity/NMDS/C_stress.txt",sep = "\t",row.names = FALSE)
write.table(result[[2]],"Results/Carbon/Beta diversity/NMDS/C_nmds.txt",sep = "\t",row.names = FALSE)
pdf(file = "Results/Carbon/Beta diversity/NMDS/C_nmds_group.pdf",width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Carbon/Beta diversity/NMDS/C_nmds_ellipse.pdf",width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Carbon/Beta diversity/NMDS/C_nmds_label.pdf",width = 7.5,height = 5.4)
result[[5]]
dev.off()
```
![Fig2](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/0dcc4dae-95ca-43a2-a1f7-64b31cd6920d)

