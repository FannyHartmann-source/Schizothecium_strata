#LD plot#
library(LDheatmap)
library(reshape2)
library(tydiverse)

#Schizothecium tetrasporum sensu stricto
LD_esp1 <- read.table(pipe("awk '{print $2, $3, $5}' VCF.Octobre2022.CBS815.71Sp3.final.diploid.VCF_final.genotyped.1.SNPs.filtered.geno0.9.only.biallelic.pol.renamed.chrom1.espece_1.pol.geno1.maf0.2.LD.output.geno.ld"), header=T)
head(LD_esp1)
summary(LD_esp1)
dim(LD_esp1)

LD.m_esp1 <- dcast(LD_esp1, POS1~POS2, value.var= "R.2")
row.names(LD.m_esp1) <- LD.m_esp1$POS1
ld.dat_esp1 <- as.matrix(LD.m_esp1 [,-1])
print(dim(ld.dat_esp1))

#Looking for a SNP in the MAt locus in order to highlight the MAT locus position
#Mat locus position =  5,182,946
MAT = LD_esp1 %>% filter(POS1>=5181946)
MAT_1 = MAT %>% filter(POS1<=5183946)
MAT_2 = unique(MAT_1$POS1)
MAT_2
# [1] 5182901 5182916 5182962 5183006 5183077 5183321 5183336 5183371 5183496
# [10] 5183654 5183698 5183776 5183935

#We take the closest SNP = 5,182,962 for the MAT locus position
#We now look for a SNP located near the start of the non-recombining region (NRR) and one near the end
#NRR = 5011904 – 6390999
NRR_MIN = LD_esp1 %>% filter(POS1>=5010904)
NRR_MIN_1 = NRR_MIN %>% filter(POS1<=5012904)
NRR_MIN_2 = unique(NRR_MIN_1$POS1)
NRR_MIN_2

# [1] 5011188 5011474 5011863
# SNP inside the NRR closest to the region border : 5011188 )
#: CBS815.71p_contig_1 : CBS815.71p_g1452:g : 5010871 – 5012936

NRR_MAX= LD_esp1 %>% filter(POS1>=6389999)
NRR_MAX_1 = NRR_MAX %>% filter(POS1<=6391999)
NRR_MAX_2 = unique(NRR_MAX_1$POS1)
NRR_MAX_2

# [1] 6390148 6390271 6390802 6391060
#SNP inside the NRR closest to the region border: 6391060 

#Looking for these SNPs rank 
esp1_pos1 = unique(LD_esp1$POS1) 
#Total position number = lenght of esp1_pos1+1 = 7276+1=7277
which(esp1_pos1 == 5011188, arr.ind=TRUE)
#[1] 4178

which(esp1_pos1 ==5182962, arr.ind=TRUE)
#[1] 4358

which(esp1_pos1 ==6391060, arr.ind=TRUE)
#[1] 6238

#plot
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")
png("esp1_LDheatmap_chrom1_final.png")
LDheatmap(ld.dat_esp1, color=rgb.palette(36), add.map=F, genetic.distances=as.integer(row.names(ld.dat_esp1)),title="LD heatmap of espece1 chrom 1",flip=TRUE)
dev.off()

#S. tritetrasporum
LD_esp3 <- read.table(pipe("awk '{print $2, $3, $5}' VCF.Octobre2022.CBS815.71Sp3.final.diploid.VCF_final.genotyped.1.SNPs.filtered.geno0.9.only.biallelic.pol.renamed.chrom1.espece_3.pol.geno1.maf0.2.LD.output.geno.ld"), header=T)
head(LD_esp3)
summary(LD_esp3)
dim(LD_esp3)

LD.m_esp3 <- dcast(LD_esp3, POS1~POS2, value.var= "R.2")
row.names(LD.m_esp3) <- LD.m_esp3$POS1
ld.dat_esp3 <- as.matrix(LD.m_esp3 [,-1])
print(dim(ld.dat_esp3))

#Looking for a SNP in the MAt locus in order to highlight the MAT locus position
#MAT = 5,182,946
MAT_sp3 = LD_esp3 %>% filter(POS1>=5181946)
MAT_1_sp3 = MAT_sp3 %>% filter(POS1<=5183946)
MAT_2_sp3 = unique(MAT_1_sp3$POS1)
MAT_2_sp3
# [1] 5183198 5183199 5183231 5183234 5183600 5183735
#Let's take the closest value for the MAT locus position: 5,183,198 
#We now look for a SNP located near the start of the non-recombining region (NRR) and one near the end
#NRR = 5138212 - 6433654

NRR_MIN_sp3 = LD_esp3 %>% filter(POS1>=5137212)
NRR_MIN_1_sp3 = NRR_MIN_sp3%>% filter(POS1<=5139212)
NRR_MIN_2_sp3 = unique(NRR_MIN_1_sp3$POS1)
NRR_MIN_2_sp3 : 5137966

#We chose the closest SNP to the region border : 5137966


NRR_MAX_sp3= LD_esp3 %>% filter(POS1>=6432654)
NRR_MAX_1_sp3 = NRR_MAX_sp3 %>% filter(POS1<=6434654)
NRR_MAX_2_sp3 = unique(NRR_MAX_1_sp3$POS1)
NRR_MAX_2_sp3

# [1] 6432764 6432791 6432877 6432930 6434062 6434297

# #We chose the closest SNP to the region border :6434062 même gène 

#Looking for these SNPs rank 
Esp3_pos1 = unique(LD_esp3$POS1) 
#Total position number = length of esp3_pos1 = 6420 + 1 = 6421 
which(Esp3_pos1 == 5183198, arr.ind=TRUE)
# [1] 3911

which(Esp3_pos1 ==5137966, arr.ind=TRUE)
#[1] 3879

which(Esp3_pos1 ==6434062, arr.ind=TRUE)
#[1] 5545

rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")
png("esp3_LDheatmap_final.png")
LDheatmap(ld.dat_esp3, color=rgb.palette(36), add.map=F, genetic.distances=as.integer(row.names(ld.dat_esp3)),title="LD heatmap of espece3 chrom 1",flip=TRUE,SNP.name = c("","",""))
dev.off()

#Identifying SNPs near the strata limits
#S.tetrasporum : strata limits : 5370425 ; 5734231
#S.tritetrasporum : strata limits : 5651364

lim1_sp1=filter(LD_esp1,POS1>=5368425)
lim1_sp1_1=filter(lim1_sp1,POS1<=5372425)
lim_sp1_2=unique(lim1_sp1_1)
lim_sp1_2
# [1] 5372238 5372361 5372422
# closest SNP : 5372238

lim_2_sp1=filter(LD_esp1,POS1>=5733231)
lim_2_sp1_1=filter(lim_2_sp1,POS1<=5735231)
lim_2_sp1_2=unique(lim_2_sp1_1$POS1)
lim_2_sp1_2
# [1] 5733874 5734295 5734437 5734543 5734771 5734870
#closest SNP = 5734295

lim1_sp3=filter(LD_esp3,POS1>=5650364)
lim1_sp3_1=filter(lim1_sp3,POS1<=5652364)
lim1_sp3_2=unique(lim1_sp3_1$POS1)
lim1_sp3_2
#closest SNP
#5651371 le seul dans le même gène

#Looking for these SNPs rank
which(esp1_pos1==5372238,arr.ind=TRUE)
#[1] 4689

which(esp1_pos1==5734295,arr.ind=TRUE)
#[1] 5700

which(Esp3_pos1==5651371,arr.ind=TRUE)
#[1] 4977

#Let's plot 

#Mean R2 in each strata
#Strata limits S. tetrasporum 
#Strata 1 = 5372238 - 5734295
#Strata 2 = 5011188 - 5372238
#Strata 3 = 5734295 – 6391060

LD_esp1_Strate_1 = filter(LD_esp1,POS1>=5372238)
LD_esp1_Strate1_1=filter(LD_esp1_Strate_1,POS1<=5734295)
LD_esp1_Strate1_2 = filter(LD_esp1_Strate1_1,POS2>=5372238)
LD_esp1_Strate_1_3= filter(LD_esp1_Strate1_2,POS2<=5734295)
R2_mean_Strate1=mean(LD_esp1_Strate_1_3$R.2)
R2_mean_Strate1
#[1] 0.7322299

LD_esp1_Strate_2=filter(LD_esp1,POS1>=5011188)
LD_esp1_Strate_2_1=filter(LD_esp1_Strate_2,POS1<=5372238)
LD_esp1_Strate2_2=filter(LD_esp1_Strate_2_1,POS2>=5011188)
LD_esp1_Strate2_3=filter(LD_esp1_Strate2_2,POS2<=5372238)
min(LD_esp1_Strate2_3$POS1)
#[1] 5011188
min(LD_esp1_Strate2_3$POS2)
#[1] 5011474
max(LD_esp1_Strate2_3$POS2)
#[1] 5372238
max(LD_esp1_Strate2_3$POS1)
#[1] 5367926
mean(LD_esp1_Strate2_3$R.2)
#[1] 0.2747952

LD_esp1_Strate3=filter(LD_esp1,POS1>=5734295 )
LD_esp1_Strate3_1=filter(LD_esp1_Strate3,POS1<=6391060  )
LD_esp1_Strate3_2=filter(LD_esp1_Strate3_1,POS2>=5734295 )
LD_esp1_Strate3_3=filter(LD_esp1_Strate3_2,POS2<=6391060)
mean(LD_esp1_Strate3_3$R.2)
#[1] 0.06612436


LD_esp1_NRR = rbind(LD_esp1_Strate_1_3,LD_esp1_Strate2_3,LD_esp1_Strate3_3)

mean(LD_esp1_NRR$R.2)
#[1] 0.5335705


#S. tritetrasporum strata limits:
#Strata 1 = 5137966 - 5651371
#Strata 2 = 5651371 – 6434062

LD_esp3_Strate1 = filter(LD_esp3,POS1>=5137966 )
LD_esp3_Strate1_1=filter(LD_esp3_Strate1,POS1<=5651371)
LD_esp3_Strate1_2=filter(LD_esp3_Strate1_1,POS2>=5137966)
LD_esp3_Strate1_3=filter(LD_esp3_Strate1_2,POS2<=5651371)

mean(LD_esp3_Strate1_3$R.2)
#[1] 0.6711749

LD_esp3_Strate2 = filter(LD_esp3,POS1>=5651371 )
LD_esp3_Strate2_1= filter(LD_esp3_Strate2,POS1<=6434062)
LD_esp3_Strate2_2=filter(LD_esp3_Strate2_1,POS2>=5651371)
LD_esp3_Strate2_3=filter(LD_esp3_Strate2_2,POS2<=6434062)

mean(LD_esp3_Strate2_3$R.2)
#[1] 0.09692784
LD_esp3_NRR=rbind(LD_esp3_Strate1_3,LD_esp3_Strate2_3)

mean(LD_esp3_NRR$R.2)
#[1] 0.5498645