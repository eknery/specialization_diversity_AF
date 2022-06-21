###
setwd("C:/Users/eduar/Documents/GitHub/specialization_endemism_AF")

install.packages("seqinr")

require(phangorn)
require(ape)
require(seqinr)
require(ggplot2)

##### carregando sequências e lista de terminais					

accd<-read.dna("0_data/loci/accd_align.fas", format = "fasta", skip = 0)
atpf<-read.dna("0_data/loci/atpf_align.fas", format = "fasta", skip = 0)
psbk<-read.dna("0_data/loci/psbk_align.fas", format = "fasta", skip = 0)
trns<-read.dna("0_data/loci/trns_align.fas", format = "fasta", skip = 0)
ets<-read.dna("0_data/loci/ets_align.fas", format = "fasta", skip = 0)
its<-read.dna("0_data/loci/its_align.fas", format = "fasta", skip = 0)

keep_list = read.table("0_data/loci/0_keep_list.csv", sep=",", h=F)

##### convertendo em caracteres
accd_char<-as.character(phyDat(accd))
atpf_char<-as.character(phyDat(atpf))
psbk_char<-as.character(phyDat(psbk))
trns_char<-as.character(phyDat(trns))
ets_char<-as.character(phyDat(ets))
its_char<-as.character(phyDat(its))

##### selecionando terminais
accd_clean<-accd_char[which(row.names(accd_char) %in% keep_list$V1),]
atpf_clean<-atpf_char[which(row.names(atpf_char) %in% keep_list$V1),]
psbk_clean<-psbk_char[which(row.names(psbk_char) %in% keep_list$V1),]
trns_clean<-trns_char[which(row.names(trns_char) %in% keep_list$V1),]
ets_clean<-ets_char[which(row.names(ets_char) %in% keep_list$V1),]
its_clean<-its_char[which(row.names(its_char) %in% keep_list$V1),]

accd_keep<-phyDat(accd_clean, type = "DNA")
atpf_keep<-phyDat(atpf_clean, type = "DNA")
psbk_keep<-phyDat(psbk_clean, type = "DNA")
trns_keep<-phyDat(trns_clean, type = "DNA")
ets_keep<-phyDat(ets_clean, type = "DNA")
its_keep<-phyDat(its_clean, type = "DNA")

# para comparar nomes
row.names(accd_clean)[order(row.names(accd_clean))]
row.names(psbk_clean)[order(row.names(psbk_clean))]

############################## analisando topologia #############################

## árvores de bootstrap
njs_accd<- bootstrap.phyDat(accd_keep,FUN=function(x)NJ(dist.ml(x,model="F81",k=6)), bs=100)
njs_atpf<- bootstrap.phyDat(atpf_keep,FUN=function(x)NJ(dist.ml(x,model="F81",k=6)), bs=100)
njs_psbk<- bootstrap.phyDat(psbk_keep,FUN=function(x)NJ(dist.ml(x,model="F81",k=6)), bs=100)
njs_trns<- bootstrap.phyDat(trns_keep,FUN=function(x)NJ(dist.ml(x,model="F81",k=6)), bs=100)
njs_ets<- bootstrap.phyDat(ets_keep,FUN=function(x)NJ(dist.ml(x,model="F81",k=6)), bs=100)
njs_its<- bootstrap.phyDat(its_keep,FUN=function(x)NJ(dist.ml(x,model="F81",k=6)), bs=100)

## comparando árvores 
all_njs<-c(njs_accd,njs_atpf,njs_psbk,njs_trns,njs_ets,njs_its)

rf<-RF.dist(all_njs, tree2 = NULL, normalize = T)
pc_rf<-pcoa(rf, correction="none", rn=NULL)
vec_pcrf<-pc_rf$vectors

locus<-c(rep("accd",length.out=100),rep("atpf",length.out=100),rep("psbk",length.out=100),rep("trns",length.out=100),rep("ets",length.out=100),rep("its",length.out=100))
df_rf<-data.frame(locus,vec_pcrf)


## plot

tiff("0_data/topology_dissimilarity.tiff", units="in", width=4, height=3, res=600)
ggplot(data = df_rf,aes(x=Axis.1,y=Axis.2,color=locus)) +
	geom_point(size = 1, alpha = 0.5) +
	scale_colour_manual(values=c("red4","orange2","blue4","cyan4","red2","hotpink1"))+
	labs(x="PCo1", y="PCo2")+
	theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=13,face="bold"))
dev.off()

####################### analisando comprimento de ramos ############################

lengths<-matrix(NA, nrow = 600, ncol = 1)

for (i in 1:length(all_njs))
	{
	lengths[i,]<-sum(all_njs[[i]]$edge.length)
	}

df_len<-as.data.frame(cbind(locus,lengths))
df_len$V2<-as.numeric(lengths)
colnames(df_len)<-c("locus","length")

boxplot(df_len$length~df_len$locus)

write.table(df_len,"total_lengths.csv",sep=",",quote=F,row.names=F)

##################### teste de modelos ##########################################

accd<-read.dna(file=file.choose(), format = "fasta", skip = 0)
atpf<-read.dna(file=file.choose(), format = "fasta", skip = 0)
psbk<-read.dna(file=file.choose(), format = "fasta", skip = 0)
trns<-read.dna(file=file.choose(), format = "fasta", skip = 0)
ets<-read.dna(file=file.choose(), format = "fasta", skip = 0)
its<-read.dna(file=file.choose(), format = "fasta", skip = 0)

pd_accd<- phyDat(accd, type = "DNA")
pd_atpf<- phyDat(atpf, type = "DNA")
pd_psbk<- phyDat(psbk, type = "DNA")
pd_trns<- phyDat(trns, type = "DNA")
pd_ets<- phyDat(ets, type = "DNA")
pd_its<- phyDat(its, type = "DNA")

mt_accd <- modelTest(pd_accd)
mt_atpf <- modelTest(pd_atpf)
mt_psbk <- modelTest(pd_psbk)
mt_trns <- modelTest(pd_trns)
mt_ets <- modelTest(pd_ets)
mt_its <- modelTest(pd_its)

mt_accd[mt_accd$AIC == min(mt_accd$AIC),]
mt_atpf[mt_atpf$AIC == min(mt_atpf$AIC),]
mt_psbk[mt_psbk$AIC == min(mt_psbk$AIC),]
mt_trns[mt_trns$AIC == min(mt_trns$AIC),]
mt_ets[mt_ets$AIC == min(mt_ets$AIC),]
mt_its[mt_its$AIC == min(mt_its$AIC),]

write.table(mt_its,"model_test_its.csv",sep=",",quote=F,row.names=F)


####################### analisando distância ###########################

d1 <- as.matrix(dist.dna(accd_keep, model="TN93",gamma=4))
d2 <- as.matrix(dist.dna(atpf_keep, model="TN93",gamma=4))
d3 <- as.matrix(dist.dna(psbk_keep, model="TN93",gamma=4))
d4 <- as.matrix(dist.dna(trns_keep, model="TN93",gamma=4))
d5 <- as.matrix(dist.dna(ets_keep, model="TN93",gamma=4))
d6 <- as.matrix(dist.dna(its_keep, model="TN93",gamma=4))

dmat<-rbind(d1,d2,d3,d4,d5,d6,d7)

CADM.global(dmat, nmat=7, n=26, nperm=9999, silent=F)
CADM.post(dmat, nmat=7, n=26, nperm=9999, mult="bonferroni", mantel=T, silent=F)


