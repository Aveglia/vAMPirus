#!/usr/bin/env R

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")

if (!require("vegan", quietly = TRUE))
  BiocManager::install("vegan")

if (!require("rstatix", quietly = TRUE))
  BiocManager::install("rstatix")

if (!require("tidyverse", quietly = TRUE))
  BiocManager::install("tidyverse")

if (!require("scales", quietly = TRUE))
  BiocManager::install("scales")

if (!require("cowplot", quietly = TRUE))
  BiocManager::install("cowplot")

if (!require("dplyr", quietly = TRUE))
  BiocManager::install("dplyr")

##Loading data
#args = commandArgs(trailingOnly=TRUE)
#sample_name=args[1]
#sample_metadata=args[2]
#options(warn=-1)
sample_name="vAMPrun_otu.55_counts.csv"
sample_metadata="pvid_samples.csv"
data<- read.csv(sample_name) #PVID_protcounts5_15_20.csv

data2 <-as.data.frame(t(data))
data2$sample <- row.names(data2)
str(data2)
colnames(data2)<- as.matrix(data2[1,])
as.data.frame(data2)
data2 <- data2[-1,]

#X.OTU.ID for X.Sequence.
data2 <- data2 %>%
  rename(Sample=X.OTU.ID)
data2dim <- dim(data2)

# test to remove index name
#row.names(data2) <- NULL

##Loading metadata
samples <- read.csv(sample_metadata)

##Combining data and metadata
data3 <- merge(data2, samples, by="Sample")

dim_data3 <- dim(data3)
dim_samples <- dim(samples)
cols <- dim_data3[2]-dim_samples[2]+1
first <-colnames(data3)[2]
last <- colnames(data3)[cols]
data3[,2:cols] <- lapply(data3[,2:cols], as.character)
data3[,2:cols] <- lapply(data3[,2:cols], as.numeric)

str(data3)

#Calculate total reads per sample
data4 <- data3%>%
  mutate(sum=select(.,2:cols)%>%
           apply(1, sum, na.rm=TRUE))

##Figures with read depth (empty!)
readdepth <- ggplot(data4, aes(x=Sample, y=sum))+
  geom_point()+
  theme_classic()+
  theme(axis.text.y = element_text(size=6))+
  coord_flip()+
  labs(y="Total reads")
readdepth

##Filter samples with low reads
data5 <- data4 %>%
  filter(sum>10000)
data5dim <-dim(data5)
minreads<-min(data5$sum)

h <- data5 %>%
  filter(treatment == "H")
mean(h$sum)
min(h$sum)
max(h$sum)

##Rarefaction curves
rarefaction <- rarecurve(data5[,2:cols])
#rarefaction

##rarefied dataset
raredata <- as.data.frame(rrarefy(data5[,2:cols], sample=minreads))

#Rarefied diversity metrics
##Diversity analysis
index <-diversity(raredata, index= "shannon")
shannondata5 <- as.data.frame(index)

index <- diversity(raredata, index= "simpson")
simpsondata5 <- as.data.frame(index)

mind5<-min(data5$sum)
index2 <- rarefy(data5[,2:cols], sample=mind5)
rarerichnessdata5 <- as.data.frame(index2)

metadata <- data5[,(cols+1):data5dim[2]]

metadata$Sample <- data5$Sample
shannondata5$Sample<- data5$Sample
simpsondata5$Sample<- data5$Sample
rarerichnessdata5$Sample <-data5$Sample

shannondata5_2 <- merge(shannondata5, metadata, by="Sample")
simpsondata5_2 <- merge(simpsondata5, metadata, by="Sample")
rarerichnessdata5_2 <- merge(rarerichnessdata5, metadata, by="Sample")

shannonplot <- ggplot(shannondata5_2, aes(x=treatment,color=treatment,y=index))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  labs(y="Shannon diversty",title="Shannon", x="Treatment")+
  theme(legend.position = "none")
simpsonplot <- ggplot(simpsondata5_2, aes(x=treatment,color=treatment,y=index))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  labs(y="Simpson diversty",title="Simpson",x="Treatment")+
  theme(legend.position = "none")
richnessplot <- ggplot(rarerichnessdata5_2, aes(x=treatment,color=treatment,y=index2))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  labs(y="ASV richness",title="Richness",x="Treatment")+
  theme(legend.position = "none")
plot_grid(shannonplot, simpsonplot,richnessplot)

##Distance
intermediate <- raredata
bray <- vegdist(sqrt(intermediate), method="bray")

##Dispersion
disper <- betadisper(bray, group = raredata$treatment, type="centroid")
boxplot(disper)

##Figure of dispersion
df <- data.frame(Distance_to_centroid=disper$distances,Group=disper$group)
df$Sample <- data5$Sample
df2 <- merge(df, metadata, by="Sample")

p<- ggplot(data=df2,aes(x=treatment,y=Distance_to_centroid,colour=treatment))+
  geom_boxplot(outlier.alpha = 0)+
  theme_classic()+
  geom_point(alpha=1/2,position=position_dodge(width=0.75))+
  labs(y="Distance to centroid")
plot(p)

##NMDS
datax <- decostand(raredata,method="total") #method 'total' normalizes data to sum up to 1 --data5[,2:cols]

MDS <- metaMDS(sqrt(datax),
               distance = "bray",autotransform = FALSE,
               k = 2,
               maxit = 999,
               trymax = 500,
               wascores = TRUE)

#Plots for goodness of fit
goodness(MDS)
stressplot(MDS)
plot(MDS, "sites")
MDS
MDS$stress

#Converting to matrix for ggplot
data.scores <- as.data.frame(scores(MDS))
data.scores$Sample <- data5$Sample
data.scores.2 <- merge(data.scores, metadata, by="Sample")

##Plot nMDS
PVID_NMDS <- ggplot(data.scores.2, aes(x=NMDS1, y=NMDS2, fill=treatment,color=treatment))+
  theme_classic()+
  geom_point(size=3, aes(shape=colony))+
  geom_polygon(data=x1C, alpha=0.1)+
  geom_polygon(data=x1H, alpha=0.1)+
  geom_polygon(data=x2C, alpha=0.1)+
  geom_polygon(data=x2H, alpha=0.1)+
  geom_polygon(data=x3H, alpha=0.1)+
  geom_polygon(data=x3C, alpha=0.1)+
  geom_polygon(data=x4H, alpha=0.1)+
  geom_polygon(data=x4C, alpha=0.1)+
  geom_polygon(data=x5H, alpha=0.1)+
  geom_polygon(data=x5C, alpha=0.1)+
  scale_color_manual(values=colors)+
  #annotate("text", label="Stress=0.12", x=-0.35, y=-0.5)+
  scale_fill_manual(values=colors)+
  labs(title="")
PVID_NMDS

##Make dataframe tidy or long
dataz <- decostand(data5[,2:cols],method="total") #method 'total' normalizes data to sum up to 1
dataz$Sample <- data5$Sample
datay <- merge(dataz, metadata, by="Sample")
datalong <- datay %>%
  tidyr::gather(first:last, key=hit, value=reads)

##Barplot
spec_bar <- ggplot(datalong, aes(x=forcats::fct_reorder(timepoint,as.numeric(as.character(timepoint))), y=reads, fill=hit))+ #Consider sqrt-transforming data
  geom_bar(aes(), stat="identity", position="fill")+
  #coord_polar("y", start=0)+
  theme_classic()+
  facet_wrap(colony~treatment, nrow=5)+
  #coord_flip()+
  labs(x="timepoint")
spec_bar #+ theme(legend.position = "none")

datalong <- datalong %>%
  filter(reads>0)
asv_bar <- ggplot(datalong, aes(x=reorder(hit,reads), y=reads,fill=treatment))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=colors)+
  coord_flip()+
  theme_classic()
asv_bar

##Sum samples <1% from data 5 for Barplot
databar <-data5
rownames(databar) <- databar[,1]
databar<-databar[,-1]
databar<-as.data.frame(databar)
str(databar)
b<-dim(databar)[2]
c<-dim(metadata)[2]
d <- c-1
e <- b-d
databar[,1:e] <- lapply(databar[,1:e], as.character)
databar[,1:e] <- lapply(databar[,1:e], as.numeric)
databar[,1:e] <- decostand(databar[,1:e],method="total")

dataj <- t(databar)
dataj2 <- dataj[1:e,]
asvnames <- rownames(dataj2)

b<-dim(dataj2)[2]
str(dataj2)
dataj2<-as.data.frame(dataj2)
dataj2[,1:b] <- lapply(dataj2[,1:b], as.character)
dataj2<-as.data.frame(dataj2)
dataj2[,1:b] <- lapply(dataj2[,1:b], as.numeric)
dataj2<-as.data.frame(dataj2)
str(dataj2)

datak <- dataj2 %>%
  mutate(sum=select(.,1:b)%>%
           apply(1, sum, na.rm=TRUE))
rownames(datak) <- asvnames

total<-sum(datak$sum)
datak <- datak %>%
  rownames_to_column("asv")%>%
  mutate(avg = sum/total)%>%
  column_to_rownames("asv")

lowab <- datak%>%
  rownames_to_column("asv")%>%
  filter(avg<0.01)%>%
  column_to_rownames("asv")
str(lowab)

highab <- datak%>%
  rownames_to_column("asv")%>%
  filter(avg>0.01)%>%
  column_to_rownames("asv")
str(highab)

dimlowab <- dim(lowab)
lowab[dimlowab[1]+1,] <- colSums(lowab)
rownames(lowab)[dimlowab[1]+1] <- "Other < 1%"
dimhighab<-dim(highab)
highab[dimhighab[1]+1,] <- lowab[dimlowab[1]+1,]
View(highab)

highab2<-t(highab)
highab2 <- as.data.frame(highab2)
dm <-dim(metadata)[2]-1

highab2$Sample <- rownames(highab2)
highab3 <- merge(highab2, metadata, by="Sample")
totalcols <- dim(highab3)[2]
x1 <-colnames(highab3)[2]
xlast <-colnames(highab3)[totalcols-dm]
highablong <- highab3 %>%
  tidyr::gather(x1:xlast, key=hit, value=reads)
