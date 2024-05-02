library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggraph)
library(igraph)
library(readr)
library(cowplot)
library(readxl)
library(rlang)
library(hexbin)
library(lattice)
library(scRepertoire)
library(ggpubr)
library(rstatix)

load("./dat.3_rename.RData")


Idents(dat.3_rename) <- 'seurat_clusters'
p1 <- DimPlot(dat.3_rename, reduction = "umap", label = TRUE, repel = TRUE) 
print(p1)

colnames(dat.3_rename@meta.data)
p2 <- DimPlot(dat.3_rename, reduction = "umap", split.by = "patient",label=TRUE)
print(p2)

p2 <- DimPlot(dat.3_rename, reduction = "umap", split.by = "treatment",label=TRUE)
print(p2)


p2 <- DimPlot(dat.3_rename, reduction = "umap", split.by = "response",label=TRUE)
print(p2)

unique(dat.3_rename$timepoint)
unique(dat.3_rename$treatment)
unique(dat.3_rename$response)
unique(dat.3_rename$disease)
unique(dat.3_rename$label.fine)
unique(dat.3_rename$cell.types)
unique(dat.3_rename$orig.ident)
unique(dat.3_rename$patient)
unique(dat.3_rename[which(dat.3_rename$orig.ident=="22_1"),]$timepoint)
colnames(dat.3_rename@meta.data)
table(dat.3_rename$patient)
dat.3_rename$group <- "Group B"
dat.3_rename@meta.data[which(dat.3_rename$patient==2),]$group <- "Group A"
dat.3_rename@meta.data[which(dat.3_rename$patient==6),]$group <- "Group A"
dat.3_rename@meta.data[which(dat.3_rename$patient==11),]$group <- "Group A"
dat.3_rename@meta.data[which(dat.3_rename$patient==16),]$group <- "Group A"
dat.3_rename@meta.data[which(dat.3_rename$patient==17),]$group <- "Group A"
dat.3_rename@meta.data[which(dat.3_rename$patient==23),]$group <- "Group A"
table(dat.3_rename$group,dat.3_rename$patient)

asdf <- colnames(dat.3_rename[,which(dat.3_rename$orig.ident=="22_1")])
tail(asdf,10)
which(str_detect(colnames(dat.3_rename),"-1"))
tcr_dir <- c("./tcr_seq/TCR-seq data/")



each_file <- list.files(tcr_dir)

each_file <- toupper(each_file)
i<-each_file[7]
vdj.list <- list()
total_cell <-0
for(i in each_file[1:length(each_file)]){
  tcr_file <- list.files(paste0(tcr_dir,i))  
  
  file_extract <- read.csv(paste0(tcr_dir,i,"/",tcr_file))
  
  vdj.list <- c(vdj.list, list(file_extract))
  total_cell <- total_cell + length(rownames(file_extract))
  
}

#for group
#each_file <- each_file[-c(7,30)]
#vdj.list[[7]] <- NULL
#vdj.list[[29]] <- NULL

patient_id <- c("15","15","17","17","22","22","9","20","20","1","1","2","2","21","21",
                "4","4","6","6","16","16","5","5","10","10","23","23","11","11","8")

#for group
#patient_id <- c("15","15","17","17","22","22","20","20","1","1","2","2","21","21",
#                "4","4","6","6","16","16","5","5","10","10","23","23","11","11")

combined <- combineTCR(vdj.list, samples = each_file, 
                        removeNA = F, removeMulti = F)

combined <- addVariable(combined, 
                        variable.name = "patient", 
                        variables = patient_id)

groupab <- c("Group B","Group B","Group A","Group A","Group B",
             "Group B","Nothing","Group B","Group B","Group B",
             "Group B","Group A","Group A","Group B","Group B",
             "Group B","Group B","Group A","Group A","Group A",
             "Group A","Group B","Group B","Group B","Group B",
             "Group A","Group A","Group A","Group A","Nothing")

#for group
#groupab <- c("Group B","Group B","Group A","Group A","Group B",
#             "Group B","Group B","Group B","Group B",
#             "Group B","Group A","Group A","Group B","Group B",
#             "Group B","Group B","Group A","Group A","Group A",
#             "Group A","Group B","Group B","Group B","Group B",
#             "Group A","Group A","Group A","Group A")

combined <- addVariable(combined, 
                        variable.name = "group", 
                       variables = groupab)

time <- c("C1D-7","C1D+1","C1D-7","C1D+1","C1D-7",
          "C1D+1","C1D-7","C1D+1","C1D-7","C1D+1",
          "C1D-7","C1D+1","C1D-7","C1D+1","C1D-7",
          "C1D+1","C1D-7","C1D+1","C1D-7","C1D+1",
          "C1D-7","C1D+1","C1D-7","C1D+1","C1D-7",
          "C1D+1","C1D-7","C1D+1","C1D-7","C1D+1")

#for group
#time <- c("C1D-7","C1D+1","C1D-7","C1D+1","C1D-7",
#          "C1D+1","C1D+1","C1D-7","C1D+1",
#          "C1D-7","C1D+1","C1D-7","C1D+1","C1D-7",
#          "C1D+1","C1D-7","C1D+1","C1D-7","C1D+1",
#          "C1D-7","C1D+1","C1D-7","C1D+1","C1D-7",
#          "C1D+1","C1D-7","C1D+1","C1D-7")

combined <- addVariable(combined, 
                        variable.name = "timepoint", 
                        variables = time)
treat <- c("Baseline","Ibrutinib","Baseline","Ibrutinib","Baseline",
           "Ibrutinib","Baseline","Ibrutinib","Baseline","Ibrutinib",
           "Baseline","Ibrutinib","Baseline","Ibrutinib","Baseline",
           "Ibrutinib","Baseline","Ibrutinib","Baseline","Ibrutinib",
           "Baseline","Ibrutinib","Baseline","Ibrutinib","Baseline",
           "Ibrutinib","Baseline","Ibrutinib","Baseline","Ibrutinib")

#for group
#treat <- c("Baseline","Ibrutinib","Baseline","Ibrutinib","Baseline",
#           "Ibrutinib","Ibrutinib","Baseline","Ibrutinib",
#           "Baseline","Ibrutinib","Baseline","Ibrutinib","Baseline",
#           "Ibrutinib","Baseline","Ibrutinib","Baseline","Ibrutinib",
#           "Baseline","Ibrutinib","Baseline","Ibrutinib","Baseline",
#           "Ibrutinib","Baseline","Ibrutinib","Baseline")

combined <- addVariable(combined, 
                        variable.name = "treatment", 
                        variables = treat)


tot_cell <- 0
for (i in seq_along(combined)) {
  tot_cell <- tot_cell + length(rownames(combined[[i]]))
}

for (i in seq_along(combined)) {
  combined[[i]]$barcode <- t(as.data.frame(strsplit(combined[[i]]$barcode,split="-")))[,1]
  combined[[i]]$barcode <- sub("_","-",combined[[i]]$barcode)
  
}

ga <- list()
gb <- list()
for (i in seq_along(combined)) {
  if(combined[[i]]$group[1]=="Group A"){
    ga <- c(ga, list(combined[[i]]))
  } else{
    gb <-c(gb, list(combined[[i]]))
  }
  
}


pre_a <- list()
post_a <- list()
for (i in seq_along(ga)) {
  if(ga[[i]]$timepoint[1]=="C1D-7"){
    pre_a <- c(pre_a, list(ga[[i]]))
  } else{
    post_a <-c(post_a, list(ga[[i]]))
  }
  
}


pre_b <- list()
post_b <- list()
for (i in seq_along(gb)) {
  if(gb[[i]]$timepoint[1]=="C1D-7"){
    pre_b <- c(pre_b, list(gb[[i]]))
  } else{
    post_b <-c(post_b, list(gb[[i]]))
  }
  
}

pre_ab <- list()
post_ab <- list()
for (i in seq_along(combined)) {
  if(combined[[i]]$timepoint[1]=="C1D-7"){
    pre_ab <- c(pre_ab, list(combined[[i]]))
  } else{
    post_ab <-c(post_ab, list(combined[[i]]))
  }
  
}


save(pre,post,combined, file="./tcr_seq_mdsc_all_patient.RData")
save(pre_ab,post_ab,pre_a,post_a,pre_b,post_b,ga,gb,combined, file="./tcr_seq_mdsc_all_groupab.RData")
load("./tcr_seq_mdsc_all_patient.RData")
load("./tcr_seq_mdsc_all_groupab.RData")
library(vegan)
com <- data.frame()
for (i in seq_along(combined)) {
  com <- rbind(com,combined[[i]])
}





pdf("./mdsc_tcr_analysis_4.pdf")
#all patient
a<-clonalDiversity(pre, 
                      cloneCall = "gene", 
                      x.axis = "timepoint", 
                      group.by = "patient",
                      metrics = c("shannon"),
                      #skip.boots = TRUE,
                      n.boots = 100)+ylim(c(3,7))+ggtitle("All patients")

b<-clonalDiversity(post, 
                      cloneCall = "gene", 
                      x.axis = "timepoint", 
                      group.by = "patient",
                      metrics = c("shannon"),
                      #skip.boots = TRUE,
                      n.boots = 100)+ylim(c(3,7))+ggtitle("All patients")


a$data$con <- "C1D-7"
b$data$con <- "C1D+1"
db <- rbind(a$data,b$data)

pv <- wilcox.test(a$data$value, y = b$data$value)$p.value

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(3,7))+ggtitle(paste0("All patient - Shannon index P-val: ",pv))+
  xlab("Timepoint") + ylab("Index score")
print(pic) 

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(3,7))+ggtitle(paste0("All patient - Shannon index"))+
  stat_compare_means()+
  xlab("Timepoint") + ylab("Index score")
print(pic) 
#a group
a<-clonalDiversity(pre_a, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("shannon"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ylim(c(3,7))+ggtitle("Group A")

b<-clonalDiversity(post_a, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("shannon"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ylim(c(3,7))+ggtitle("Group A")

a$data$con <- "C1D-7"
b$data$con <- "C1D+1"
db <- rbind(a$data,b$data)

pv <- wilcox.test(a$data$value, y = b$data$value)$p.value

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(3,7))+ggtitle(paste0("Group A - Shannon index P-val: ",pv))+
  xlab("Timepoint") + ylab("Index score")
print(pic) 

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(3,7))+ggtitle(paste0("Group A - Shannon index"))+
  xlab("Timepoint") + ylab("Index score")+
  stat_compare_means()
print(pic) 
#b group
c<-clonalDiversity(pre_b, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("shannon"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ylim(c(3,7))+ggtitle("Group B")

d<-clonalDiversity(post_b, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("shannon"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ylim(c(3,7))+ggtitle("Group B")

c$data$con <- "C1D-7"
d$data$con <- "C1D+1"
db <- rbind(c$data,d$data)

pv <- wilcox.test(c$data$value, y = d$data$value)$p.value
pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(3,7))+ggtitle(paste0("Group B - Shannon Index P-val: ",pv))+
  xlab("Timepoint") + ylab("Index score")
print(pic) 

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(3,7))+ggtitle(paste0("Group B - Shannon Index"))+
  xlab("Timepoint") + ylab("Index score")+
  stat_compare_means()
print(pic) 

a$data$group <- "Group A"
b$data$group <- "Group A"
c$data$group <- "Group B"
d$data$group <- "Group B"


a$data$gp <- "C1D-7 Group A"
b$data$gp <- "C1D+1 Group A"
c$data$gp <- "C1D-7 Group B"
d$data$gp <- "C1D+1 Group B"

db2 <- rbind(a$data,b$data,c$data,d$data)

pv<-pairwise.wilcox.test(db2$value, db2$gp, p.adjust.method="bonf")$p.value

pic <- ggplot(db2 %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+
  geom_boxplot(aes(fill=group))+
  geom_point(aes(y=value,color=patient))+
  ylim(c(3,7))+ggtitle(paste0("Group AB - Shanon Index"))+
  xlab("Timepoint") + ylab("Index score") 

print(pic) 

pic <- ggplot(db2 %>% mutate(gp = fct_relevel(gp, "C1D-7 Group A","C1D-7 Group B", "C1D+1 Group A","C1D+1 Group B")),aes(x=as.factor(gp),y=value))+
  geom_boxplot()+
  geom_point(aes(y=value,color=patient))+
  ylim(c(3,10))+ggtitle(paste0("Group AB - Shanon Index: adjusted p-value(pair)"))+
  xlab("Timepoint") + ylab("Index score")
#print(pic)


my_comp <- list( c("C1D-7 Group A", "C1D+1 Group B"), c("C1D-7 Group A", "C1D-7 Group B"),
                 c("C1D-7 Group A", "C1D+1 Group A"), c("C1D-7 Group B", "C1D+1 Group B"),
                 c("C1D-7 Group B", "C1D+1 Group A"), c("C1D+1 Group A", "C1D+1 Group B"))


stat.test <- db2 %>%
  wilcox_test(value ~ gp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


stat.test <- stat.test %>%
  mutate(y.position = c(6.4,6.9,7.4,7.9,8.4,8.9))
print(pic + stat_pvalue_manual(stat.test, label = "p.adj") +  stat_compare_means(method="kruskal.test"))


#all patient gini simpson
a<-clonalDiversity(pre, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("gini.simpson"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ggtitle("All patients")+ylim(c(0.8,1.2))

b<-clonalDiversity(post, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("gini.simpson"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ggtitle("All patients")+ylim(c(0.8,1.2))

a$data$con <- "C1D-7"
b$data$con <- "C1D+1"
db <- rbind(a$data,b$data)

pv<-wilcox.test(a$data$value, y = b$data$value)$p.value
pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.2))+ggtitle(paste0("All patients - Gini-Simpson P-val: ",pv))+
  xlab("Timepoint") + ylab("Index score")
print(pic) 

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.2))+ggtitle(paste0("All patients - Gini-Simpson"))+
  xlab("Timepoint") + ylab("Index score")+
  stat_compare_means()
print(pic) 

#A group gini simpson
a<-clonalDiversity(pre_a, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("gini.simpson"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ggtitle("Group A")+ylim(c(0.8,1.2))

b<-clonalDiversity(post_a, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("gini.simpson"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ggtitle("Group A")+ylim(c(0.8,1.2))
a$data$con <- "C1D-7"
b$data$con <- "C1D+1"
db <- rbind(a$data,b$data)


pv<-wilcox.test(a$data$value, y = b$data$value)$p.value
pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.2))+ggtitle(paste0("Group A - Gini-Simpson P-val: ",pv))+
  xlab("Timepoint") + ylab("Index score")
print(pic) 

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.2))+ggtitle(paste0("Group A - Gini-Simpson"))+
  xlab("Timepoint") + ylab("Index score")+
  stat_compare_means()
print(pic) 

#B group
c<-clonalDiversity(pre_b, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("gini.simpson"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ggtitle("Group B")+ylim(c(0.8,1.2))

d<-clonalDiversity(post_b, 
                   cloneCall = "gene", 
                   x.axis = "timepoint", 
                   group.by = "patient",
                   metrics = c("gini.simpson"),
                   #skip.boots = TRUE,
                   n.boots = 100)+ggtitle("Group B")+ylim(c(0.8,1.2))
c$data$con <- "C1D-7"
d$data$con <- "C1D+1"
db <- rbind(c$data,d$data)

pv <- wilcox.test(c$data$value, y = d$data$value)$p.value
pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.2))+ggtitle(paste0("Group B - Gini-Simpson P-val: ",pv))+
  xlab("Timepoint") + ylab("Index score")
print(pic) 

pic <- ggplot(db %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+geom_boxplot()+geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.2))+ggtitle(paste0("Group B - Gini-Simpson"))+
  xlab("Timepoint") + ylab("Index score")+
  stat_compare_means()
print(pic) 

a$data$group <- "Group A"
b$data$group <- "Group A"
c$data$group <- "Group B"
d$data$group <- "Group B"

a$data$gp <- "C1D-7 Group A"
b$data$gp <- "C1D+1 Group A"
c$data$gp <- "C1D-7 Group B"
d$data$gp <- "C1D+1 Group B"

db2 <- rbind(a$data,b$data,c$data,d$data)
pv <- pairwise.wilcox.test(db2$value, db2$gp, p.adjust.method="bonf")$p.value


my_comp <- list( c("C1D-7 Group A", "C1D+1 Group B"), c("C1D-7 Group A", "C1D-7 Group B"),
                 c("C1D-7 Group A", "C1D+1 Group A"), c("C1D-7 Group B", "C1D+1 Group B"),
                 c("C1D-7 Group B", "C1D+1 Group A"), c("C1D+1 Group A", "C1D+1 Group B"))

pic <- ggplot(db2 %>% mutate(con = fct_relevel(con, "C1D-7", "C1D+1")),aes(x=as.factor(con),y=value))+
  geom_boxplot(aes(fill=group))+
  geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.2))+ggtitle(paste0("Group AB - Gini-Simpson"))+
  xlab("Timepoint") + ylab("Index score")

print(pic) 


my_comp <- list( c("C1D-7 Group A", "C1D+1 Group B"), c("C1D-7 Group A", "C1D-7 Group B"),
                 c("C1D-7 Group A", "C1D+1 Group A"), c("C1D-7 Group B", "C1D+1 Group B"),
                 c("C1D-7 Group B", "C1D+1 Group A"), c("C1D+1 Group A", "C1D+1 Group B"))
pic <- ggplot(db2 %>% mutate(gp = fct_relevel(gp, "C1D-7 Group A","C1D-7 Group B", "C1D+1 Group A","C1D+1 Group B")),aes(x=as.factor(gp),y=value))+
  geom_boxplot()+
  geom_point(aes(y=value,color=patient))+
  ylim(c(0.8,1.3))+ggtitle(paste0("Group AB - Gini-Simpson: adjusted p-value(pair)"))+
  xlab("Timepoint") + ylab("Index score")
#print(pic) 



stat.test <- db2 %>%
  wilcox_test(value ~ gp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


stat.test <- stat.test %>%
  mutate(y.position = c(1.05,1.09,1.13,1.17,1.21,1.25))

print(pic + stat_pvalue_manual(stat.test, label = "p.adj") +  stat_compare_means(method="kruskal.test"))

#all patient
c <- clonalHomeostasis(pre, cloneCall = "gene",group.by = "timepoint")+ggtitle("All patients")
d <- clonalHomeostasis(post, cloneCall = "gene",group.by = "timepoint")+ggtitle("All patients")

print(c+d)

#group a
c<-clonalHomeostasis(pre_a, cloneCall = "gene",group.by = "timepoint")+ggtitle("Group A")
d<- clonalHomeostasis(post_a, cloneCall = "gene",group.by = "timepoint")+ggtitle("Group A")
print(c+d)

#group b
c<- clonalHomeostasis(pre_b, cloneCall = "gene",group.by = "timepoint")+ggtitle("Group B")
d<- clonalHomeostasis(post_b, cloneCall = "gene",group.by = "timepoint")+ggtitle("Group B")

print(c+d)

#together
c<- clonalHomeostasis(ga, cloneCall = "gene",group.by = "group")+ggtitle("Group A")
d<- clonalHomeostasis(gb, cloneCall = "gene",group.by = "group")+ggtitle("Group B")

print(c+d)
dev.off()

dat.3_rename$cell.types
Idents(dat.3_rename) <- 'cell.types'
p1 <- DimPlot(dat.3_rename, reduction = "umap", label = TRUE, repel = TRUE) 


pdf("./mdsc_cell_annotation.pdf")
print(p1)
dev.off()


#pdf("./tcr_seq_analysis.pdf")

querydata <- combineExpression(combined, dat.3_rename, cloneCall = "gene", group.by = "sample",proportion = FALSE, 
                               cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
save(querydata, file="./tcr_seq_seurat_integrated_mdsc.RData")

load("./tcr_seq_seurat_integrated_mdsc.RData")
colnames(querydata@meta.data)
length(unique(querydata$patient))

#total frequency of the clonetype
print(DimPlot(querydata, group.by = "patient"))

unique(querydata$patient)
unique(querydata$patient[which(querydata$orig.ident=="C6")])
print(clonalDiversity(querydata, 
                      cloneCall = "gene", 
                      x.axis = "sample",group.by="response", 
                      metrics = "shannon",
                      n.boots = 100))



pdf("./anno_single_cell_tcr_result.pdf")
print(clonalDiversity(combined2, cloneCall = "nt")+ggtitle("from nucleotide sequence"))
print(clonalHomeostasis(combined2, cloneCall = "nt")+ggtitle("from nucleotide sequence"))
print(clonalProportion(combined2, cloneCall = "nt")+ggtitle("from nucleotide sequence"))
print(clonalOverlap(combined2, cloneCall="aa", method="overlap")+ggtitle("from nucleotide sequence"))
dev.off()


tcr_data <- rbind(combined[[1]],combined[[2]],combined[[3]])
save(combined, file="./tcr_seq_newcell.RData")
