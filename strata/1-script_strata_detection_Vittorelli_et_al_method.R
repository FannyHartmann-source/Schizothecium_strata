#Strata detection based on Vittorelli et al, 2023 method#

library(dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# {r Function to code for strata, echo=TRUE}

var_het <- function(df){
  
  ### For the dataframe df containing a list of orthogroups, this list is sequentially divided 
  ### into two segments, with the limit sliding along orthogroups, and the difference of average 
  ### heterozygosity between the two segments is computed. The partitions containing a segment with less than
  ### 10 orthogroups are excluded from the analyses. 
  
  # the columns of the database
  vmean1 <- c()
  vmean2 <- c()
  pos_split <- c()
  vdim1 <- c()
  vdim2 <- c()
  
  # looping over the orthogroups
  for(i in 1:dim(df)[1]){
    
    # position of the split between the two segments
    pos_split <- c(pos_split,df$start_matp[i])
    
    # the 1st segment
    data_region1 = df[1:i,]
    vdim1 <- c(vdim1,dim(data_region1)[1])      # number of orthogroups in the 1st segment
    mean1 <- mean(data_region1$SNPs_par_pb ,na.rm=T) # average heterozygosity in the 1st segment
    vmean1 <- c(vmean1,mean1)
    
    # the 2nd segment
    data_region2 = df[(i+1):dim(df)[1],]
    vdim2 <- c(vdim2,dim(data_region2)[1])      # number of orthogroups in the 2nd segment
    mean2 <- mean(data_region2$SNPs_par_pb ,na.rm=T) # average heterozygosity in the 1st segment
    vmean2 <- c(vmean2,mean2) 
  }
  
  # the database
  datab <- data.frame(pos_split,vmean1,vmean2,vdim1,vdim2) %>%
    # computing the difference of heterozygosity between the two segments
    mutate(Difference_het=vmean1-vmean2) %>%
    # removing the partitions in which one of the two segments had less than 10 genes
    filter(vdim1 >= 20, vdim2 >= 20)
  
  return(datab)
}
path <- "$path_to_strains"

####Order_1####
#strains_order_1 = c("S1.PSN779.h","S1.CBS81571","S5.PSQ32","S6.PSN1067.h","S13.PSN1042","S15.PSN377","S3.PSN1070","S10.PSN684.h","S4.PSN936.h")
for (i in strains_order_1) {
 
  df <- read.table(paste0(path,"/VERIF_CORRECTION_souche_SNPsSyn_avecHET",i,".txt"),sep=" ",header=T)
  head(df)
  df$scaffold_matp=factor(df$CHROM)
  df$scaffold_matm=factor(df$CHROM)
  df$start_matp = as.numeric(df$POS_GENE)
  df <- df[order(df$POS_GENE),] 
  head(df)
  summary(df)
  dim(df)
  scaffold1 <- df 
  
  ### the limit of the putative non-recombining region
  
  NRR <- read.csv(paste0(path,"/SIZE_NRR.txt"),sep="")
  head(NRR) 
  summary(NRR)
  dim(NRR)
  
  NRR_souche <- subset(NRR, SOUCHE == as.character(i))
  head(NRR_souche)
  summary(NRR_souche)
  dim(NRR_souche)
  
  start_NRR <- as.numeric(as.character(NRR_souche$MIN))
  end_NRR <-  as.numeric(as.character(NRR_souche$MAX))
  
  putative_NRR <- scaffold1 %>%
    filter(start_matp>=start_NRR,
           start_matp<=end_NRR)
  
  summary(putative_NRR)
  dim(putative_NRR)
  
  
  #{r plot: heterozygosity between MAT1_1 and MAT1_2, echo=TRUE}
  het_plot <- ggplot(putative_NRR) +
    geom_point(aes(x=start_matp/1000000, y=SNPs_par_pb ), 
               alpha=0.3, size=0.1) +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=0, ymax= 0.01 , fill="red") +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    annotate("text",x=5.15,y=0.01,label="mating-type locus",angle=90, color="red") +
    ylab("Nber of het SNPs per gene") +
    xlim(4.9,6.5)
  
  #
  
  ### Testing for four strata
  
  #`{r Testing for 4 strata: boundary between strata 2 and 3, echo=TRUE}
  # selecting the orthogroups in the putative strata 2 and 3
  putative_strata_2_and_3 <- scaffold1 %>%
    filter(start_matp>=5200000,
           start_matp<=end_NRR)
  
  # creating the database with the variation of heterozygosity
  var_het_strata_2_3 <- var_het(putative_strata_2_and_3)
  summary(var_het_strata_2_3)
  
  # position of the boundary between strata 2 and 3
  max_var <- max(abs(var_het_strata_2_3$Difference_het))
  pos_max<-(abs(var_het_strata_2_3$Difference_het)==max_var)
  boundary_strata_2_3 <- var_het_strata_2_3$pos_split[pos_max]
  boundary_strata_2_3
  
  # plotting the results
  variation_het_plot_2_3 <- ggplot(var_het_strata_2_3) +
    geom_point(aes(x= pos_split/1000000, y=Difference_het), size=0.3) +
    geom_line(aes(x= pos_split/1000000, y=Difference_het), alpha=0.3)  +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=-1,ymax=1 , fill="red") +
    geom_vline(xintercept=boundary_strata_2_3/1000000, linetype="dashed", alpha=0.5) +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    ylab("Variation of SNPs Nber") +
    xlim(4.9,6.5)
  
  variation_het_plot_2_3 +
    ggtitle("Variation of heterozygosity between the haplotypes in the putative \n non-recombining region") +
    annotate("text",x=5.15,y=0,label="mating-type locus",angle=90, color="red") 
  
  ggsave(paste0(path,"/Variation-heterozygosity-CBS815.71_strata_2_3",i,".pdf"), width=6,height=4)
  
  #
  
  #{r Testing for 4 strata: boundary between strata 1 and 2, echo=TRUE}
  # selecting the orthogroups in the putative strata 1 and 2
  putative_strata_1_and_2 <- scaffold1 %>%
    filter(start_matp>=start_NRR,
           start_matp<=boundary_strata_2_3)
  
  # creating the database with the variation of heterozygosity
  var_het_strata_1_2 <- var_het(putative_strata_1_and_2)
  summary(var_het_strata_1_2)
  
  # position of the boundary between strata 1 and 2
  min_var <- max(abs(var_het_strata_1_2$Difference_het))
  pos_min<-(abs(var_het_strata_1_2$Difference_het)==min_var)
  boundary_strata_1_2 <- var_het_strata_1_2$pos_split[pos_min]
  boundary_strata_1_2
  
  # plotting the results
  variation_het_plot_1_2 <- ggplot(var_het_strata_1_2) +
    geom_point(aes(x= pos_split/1000000, y=Difference_het), size=0.3) +
    geom_line(aes(x= pos_split/1000000, y=Difference_het), alpha=0.3)  +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=-1,ymax=1 , fill="red") +
    geom_vline(xintercept=boundary_strata_1_2/1000000, linetype="dashed", alpha=0.5) +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    ylab("Variation of SNPs Nber") +
    xlim(4.9,6.5)
  
  variation_het_plot_1_2 +
    ggtitle("Variation of heterosigosity between the haplotypes in the putative \n non-recombining region") +
    annotate("text",x=5.15,y=-0.006,label="mating-type locus",angle=90, color="red") 
  
  ggsave(paste0(path,"/Variation-heterozygosity-CBS815.71",i,"_strata_1_2",i,".pdf"),width=6,height=4)
  
  # {r Testing for 4 strata: boundary between strata 0 and 1, echo=TRUE}
  # selecting the orthogroups in the putative strata 1 and 2
  putative_strata_0_and_1 <- scaffold1 %>%
    filter(start_matp>=start_NRR,
           start_matp<=boundary_strata_1_2)
  
  # creating the database with the variation of heterozigosity
  var_het_strata_0_1 <- var_het(putative_strata_0_and_1)
  summary(var_het_strata_0_1)
  
  #position of the boundary between strata 1 and 2
  min_var <- max(abs(var_het_strata_0_1$Difference_het))
  pos_min<-(abs(var_het_strata_0_1$Difference_het)==min_var)
  boundary_strata_0_1 <- var_het_strata_0_1$pos_split[pos_min]
  boundary_strata_0_1
  
  # plotting the results
  variation_het_plot_0_1 <- ggplot(var_het_strata_0_1) +
    geom_point(aes(x= pos_split/1000000, y=Difference_het), size=0.3) +
    geom_line(aes(x= pos_split/1000000, y=Difference_het), alpha=0.3)  +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=-1,ymax=1 , fill="red") +
    geom_vline(xintercept=boundary_strata_0_1/1000000, linetype="dashed", alpha=0.5) +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    ylab("Variation of SNPs Nber") +
    xlim(4.9,6.5)
  
  variation_het_plot_0_1 +
    ggtitle("Variation of heterozygosity between the haplotypes in the putative \n non-recombining region") +
    annotate("text",x=5.15,y=-0.006,label="mating-type locus",angle=90, color="red") 
  
  ggsave(paste0(path,"/Variation-heterozygosity-CBS815.71",i,"strata_0_1",i,".pdf"), width=6,height=4)
  
  ## Summary figure about strata
  
  #{r Strata summary figure, echo=TRUE}
  
  het_plot_annotated <- het_plot +
    geom_vline(xintercept=boundary_strata_1_2/1000000, linetype="dashed", alpha=0.5) +
    geom_vline(xintercept=boundary_strata_2_3/1000000, linetype="dashed", alpha=0.5) +
    geom_vline(xintercept=boundary_strata_0_1/1000000, linetype="dashed", alpha=0.5) 
  
  
  final_plot <- grid.arrange(het_plot_annotated, variation_het_plot_2_3, variation_het_plot_1_2, variation_het_plot_0_1,  ncol=1) 
  
  ggsave(filename=paste0(path,"/Plot_4strata",i,".pdf"), plot=final_plot, width=6,height=8)
  
  print(boundary_strata_2_3)
  print(boundary_strata_1_2)
  print(boundary_strata_0_1)
  
  NRR_souche$limit1<- boundary_strata_2_3
  NRR_souche$limit2<- boundary_strata_1_2
  NRR_souche$limit3<- boundary_strata_0_1
  
  write.csv(NRR_souche, paste0(path,"/boundary_strata_strain_",i,".csv"),quote=F, row.names=F)
  
  write.table(NRR_souche, paste0(path,"/boundary_strata_strain_",i,".txt"),quote=F, row.names=F,col.names=F)
}

####Order_2#####
#strains_order_2=c("S2.CBS26270.h","S8.PSN644","S9.PSN707","S12.PSN663","S16.PSN1310.h")

for (i in strains_order_2)
{

  df <- read.table(paste0(path,"/VERIF_CORRECTION_souche_SNPsSyn_avecHET",i,".txt"),sep=" ",header=T)
  
  head(df)
  df$scaffold_matp=factor(df$CHROM)
  df$scaffold_matm=factor(df$CHROM)
  df$start_matp = as.numeric(df$POS_GENE)
  df <- df[order(df$POS_GENE),] 
  head(df)
  summary(df)
  dim(df)
  scaffold1 <- df 
  
  
  ### the limit of the putative non-recombining region
  NRR <- read.csv(paste0(path,"/SIZE_NRR.txt"),sep="")
  head(NRR)
  summary(NRR)
  dim(NRR)
  
  NRR_souche <- subset(NRR, SOUCHE == as.character(i)) 
  head(NRR_souche)
  summary(NRR_souche)
  dim(NRR_souche)
  
  start_NRR <- as.numeric(as.character(NRR_souche$MIN))
  end_NRR <-  as.numeric(as.character(NRR_souche$MAX))
  
  putative_NRR <- scaffold1 %>%
    filter(start_matp>=start_NRR,
           start_matp<=end_NRR)
  
  summary(putative_NRR)
  dim(putative_NRR)
  
  
  
  #{r plot: heterozygosity between MAT1_1 and MAT1_2 contig1, echo=TRUE}
  het_plot <- ggplot(putative_NRR) +
    geom_point(aes(x=start_matp/1000000, y=SNPs_par_pb ), 
               alpha=0.3, size=0.1) +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=0, ymax= 0.01 , fill="red") +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    annotate("text",x=5.15,y=0.01,label="mating-type locus",angle=90, color="red") +
    ylab("Nber of het SNPs per gene") +
    xlim(4.9,6.5)
  
  #
  
  ### Testing for four strata
  
  
  #`{r Testing for 4 strata: boundary between strata 2 and 3, echo=TRUE}
  # selecting the orthogroups in the putative strata 2 and 3
  putative_strata_2_and_3 <- scaffold1 %>%
    filter(start_matp>=5200000,
           start_matp<=end_NRR)
  
  # creating the database with the variation of heterozygosity
  var_het_strata_2_3 <- var_het(putative_strata_2_and_3)
  summary(var_het_strata_2_3)
  
  # position of the boundary between strata 2 and 3
  max_var <- max(abs(var_het_strata_2_3$Difference_het))
  pos_max<-(var_het_strata_2_3$Difference_het==max_var)
  boundary_strata_2_3 <- var_het_strata_2_3$pos_split[pos_max]
  boundary_strata_2_3
  
  # plotting the results
  variation_het_plot_2_3 <- ggplot(var_het_strata_2_3) +
    geom_point(aes(x= pos_split/1000000, y=Difference_het), size=0.3) +
    geom_line(aes(x= pos_split/1000000, y=Difference_het), alpha=0.3)  +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=-1,ymax=1 , fill="red") +
    geom_vline(xintercept=boundary_strata_2_3/1000000, linetype="dashed", alpha=0.5) +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    ylab("Variation of SNPs Nber") +
    xlim(4.9,6.5)
  
  variation_het_plot_2_3 +
    ggtitle("Variation of heterozygosity between the haplotypes in the putative \n non-recombining region") +
    annotate("text",x=5.15,y=0,label="mating-type locus",angle=90, color="red") 
  
  ggsave(paste0(path,"/Variation-het-CBS815.71_strata_2_3_",i,".pdf"), width=6,height=4)
  
  
  
  #
  
  #{r Testing for 4 strata: boundary between strata 1 and 2, echo=TRUE}
  # selecting the orthogroups in the putative strata 1 and 2
  putative_strata_1_and_2 <- scaffold1 %>%
    filter(start_matp>=start_NRR,
           start_matp<=boundary_strata_2_3)
  
  # creating the database with the variation of heterozygosity
  var_het_strata_1_2 <- var_het(putative_strata_1_and_2)
  summary(var_het_strata_1_2)
  
  # position of the boundary between strata 1 and 2
  min_var <- max(abs(var_het_strata_1_2$Difference_het))
  pos_min<-(abs(var_het_strata_1_2$Difference_het)==min_var)
  boundary_strata_1_2 <- var_het_strata_1_2$pos_split[pos_min]
  boundary_strata_1_2
  
  # plotting the results
  variation_het_plot_1_2 <- ggplot(var_het_strata_1_2) +
    geom_point(aes(x= pos_split/1000000, y=Difference_het), size=0.3) +
    geom_line(aes(x= pos_split/1000000, y=Difference_het), alpha=0.3)  +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=-1,ymax=1 , fill="red") +
    geom_vline(xintercept=boundary_strata_1_2/1000000, linetype="dashed", alpha=0.5) +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    ylab("Variation of SNPs Nber") +
    xlim(4.9,6.5)
  
  variation_het_plot_1_2 +
    ggtitle("Variation of heterozygosity between the haplotypes in the putative \n non-recombining region") +
    annotate("text",x=5.15,y=-0.006,label="mating-type locus",angle=90, color="red") 
  
  ggsave(paste0(path,"/Variation-heterozygosity-CBS815.71_strata_1_2_",i,".pdf"), width=6,height=4)
  
  
  
  
  # {r Testing for 4 strata: boundary between strata 0 and 1, echo=TRUE}
  # selecting the orthogroups in the putative strata 1 and 2
  
  
  
  putative_strata_0_and_1 <- scaffold1 %>%
    filter(start_matp>=boundary_strata_1_2,
           start_matp<=boundary_strata_2_3 )
  
  
  # creating the database with the variation of heterozygosity
  var_het_strata_0_1 <- var_het(putative_strata_0_and_1)
  summary(var_het_strata_0_1)
  
  # position of the boundary between strata 1 and 2
  min_var <- max(abs(var_het_strata_0_1$Difference_het))
  pos_min<-(abs(var_het_strata_0_1$Difference_het)==min_var)
  boundary_strata_0_1 <- min(var_het_strata_0_1$pos_split[pos_min])
  boundary_strata_0_1
  
  # plotting the results
  variation_het_plot_0_1 <- ggplot(var_het_strata_0_1) +
    geom_point(aes(x= pos_split/1000000, y=Difference_het), size=0.3) +
    geom_line(aes(x= pos_split/1000000, y=Difference_het), alpha=0.3)  +
    geom_rect(xmin=5179614/1000000, xmax=5180779/1000000,ymin=-1,ymax=1 , fill="red") +
    geom_vline(xintercept=boundary_strata_0_1/1000000, linetype="dashed", alpha=0.5) +
    xlab("Position on the contig 1 of the CBS815.71-sp3 assembly (Mb)") +
    ylab("Variation of SNPs Nber") +
    xlim(4.9,6.5)
  
  variation_het_plot_0_1 +
    ggtitle("Variation of heterozygosity between the haplotypes in the putative \n non-recombining region") +
    annotate("text",x=5.15,y=-0.006,label="mating-type locus",angle=90, color="red") 
  
  ggsave(paste0(path,"/Variation-heterozygosity-CBS815.71_strata_0_1_",i,".pdf"), width=6,height=4)
  
  
  ## Summary figure about strata
  

  het_plot_annotated <- het_plot +
    geom_vline(xintercept=boundary_strata_1_2/1000000, linetype="dashed", alpha=0.5) +
    geom_vline(xintercept=boundary_strata_2_3/1000000, linetype="dashed", alpha=0.5) +
    geom_vline(xintercept=boundary_strata_0_1/1000000, linetype="dashed", alpha=0.5) 
  
  
  final_plot <- grid.arrange(het_plot_annotated, variation_het_plot_2_3, variation_het_plot_1_2, variation_het_plot_0_1,  ncol=1) 
  
  ggsave(filename=paste0(path,"/Plot_4strata",i,".pdf"), plot=final_plot, width=6,height=8)
  
  print(boundary_strata_2_3)
  print(boundary_strata_1_2)
  print(boundary_strata_0_1)
  
  NRR_souche$limit1<- boundary_strata_2_3
  NRR_souche$limit2<- boundary_strata_1_2
  NRR_souche$limit3<- boundary_strata_0_1
  
  write.csv(NRR_souche, paste0(path,"/boundary_strata_souche_",i,".csv"),quote=F, row.names=F)
  
  write.table(NRR_souche, paste0(path,"/boundary_strata_souche_",i,".txt"),quote=F, row.names=F,col.names=F)
  
  
}
