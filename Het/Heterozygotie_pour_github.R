library(tidyverse)
library(ggplot2)
library(paletteer)
library(gridExtra)


#########Boucle hétérokaryons

souches = c("S13.CBS13294.h","S3.PSN1040.h","S10.PSN684.h","S3.PSN1070.h","S6.PSN1067.h","S2.CBS26270.h","S1.CBS31562.h","S1.PSN1217.h","S13.PSN1219.h","S3.PSN1236.h","S16.PSN1310.h","S3.PSN1433.h","S3.PSN1434.h","S3.PSN352.h","S3.PSN571.h","S3.PSN667.h","S1.PSN710.h","S3.PSN777.h","S1.PSN779.h","S3.PSN829.h","S4.PSN936.h","S13.PSN1042.h")

for (i in souches ) {
  
  #Ouvrir la liste des SNPs qui sont dans les gènes (ils sont déjà sélectionnés mais cela permet d'avoir la position du gène etc)
  
  SNP_GENES<- read.table("/Volumes/Elsa_Donnees/Octobre_2022/MAT+/VCF_final/VCF_final.strates_avecHET.Bedtools_intersect_SNPs_genes.bis.bed", header=F)
  colnames(SNP_GENES) <- c("CHROM", "POS", "VA","VB","DEBUT","FIN","VC","GENE","VD","VE")
  SNP_GENES_1=select(SNP_GENES, -c(3,4,7,9,10))
  SNP_GENES_2=select(SNP_GENES_1,-c(2))
  
  #Ouvrir les fichier de SNPs : fichiers avec tous les SNPs  
  
  SNP <- read.table(paste0("/Volumes/Elsa_Donnees/Octobre_2022/MAT+/VCF_final/souches_avecHET/Genotype.dikaryotic.Stetrasporum.",i,".VCF_final.Genes_avecHET.sansMAAF.txt"),header=T)
  colnames(SNP) <- c("CHROM", "POS", "TYPE", "REF", "VERSION")
  
  SNP_1=filter(SNP,VERSION!=".")
  SNP_2 = filter(SNP_1,VERSION!="./.")
  
  #Ouvir les fichiers avec SNPs synonymes
  
  SNPs_syn<- read.table("/Volumes/Elsa_Donnees/Octobre_2022/MAT+/VCF_final/VCF.Octobre2022.CBS815.71Sp3.final.diploid.VCF_final.genotyped.1.SNPs.filtered.geno0.9.only.biallelic.pol.renamed.strates_avecHET.genes.recode.annotated.effects.txt", header=F,col.names = paste0("V",seq_len(21)), fill = TRUE)
  colnames(SNPs_syn)=c("CHROM","POS","REF","ALT","EFFECT_1","EFFECT_2","EFFECT_3","EFFECT_4","EFFECT_5","EFFECT_6","EFFECT_7","EFFECT_8","EFFECT_9","EFFECT_10","EFFECT_11","EFFECT_12","EFFECT_13","EFFECT_14","EFFECT_15","EFFECT_16","EFFECT_17")
  SNPs_syn_1 = SNPs_syn %>% slice(-1)
  
  #on ne garde que les synonymous variant en colonne 1
  SNPs_syn_Eff1 = SNPs_syn_1  %>% filter(EFFECT_1=="synonymous_variant")
  dim(SNPs_syn_Eff1) #1451286
  
  #CHROMOSOME 1 
  SNP_syn_chrom1 = SNPs_syn_Eff1 %>% 
    filter(CHROM=="CBS815.71p_contig_1")
  
  SNP_3_chrom1=SNP_2 %>% 
    filter(CHROM == "CBS815.71p_contig_1")
  
  SNP_GENES_chrom1 = SNP_GENES_1 %>% 
    filter(CHROM=="CBS815.71p_contig_1")
  
  SNP_4_chrom1=merge(x = SNP_3_chrom1, y = SNP_GENES_1, all.x=TRUE)
  
  #On sélectionne uniquement les SNPs synonymes
  
  SNP_4bis_chrom1 = merge(x = SNP_4_chrom1, y = SNP_syn_chrom1, all.x=FALSE, all.y=FALSE)
  
  SNP_4ter_chrom1= SNP_4bis_chrom1 %>% arrange(POS,GENE)
  
  SNP_5_chrom1=SNP_4bis_chrom1 %>% distinct(POS,.keep_all = TRUE)
  
  #On ajoute une colonne score
  SNP_6_chrom1=SNP_5_chrom1 %>%
    mutate(score = case_when(
      VERSION=="A/A" ~ 0,
      VERSION=='A|A' ~ 0,
      VERSION== "C/C" ~ 0,
      VERSION=='C|C' ~ 0,
      VERSION=='G/G' ~ 0,
      VERSION=='G|G' ~ 0,
      VERSION=='T/T' ~ 0,
      VERSION=='T|T' ~ 0,
      VERSION=='A/C' ~ 1,
      VERSION=='A|C' ~ 1,
      VERSION=='A/G' ~ 1,
      VERSION=='A|G' ~ 1,
      VERSION=='A/T' ~ 1,
      VERSION=='A|T' ~ 1,
      VERSION=='C/A' ~ 1,
      VERSION=='C|A' ~ 1,
      VERSION=='C/G' ~ 1,
      VERSION=='C|G' ~ 1,
      VERSION=='C/T' ~ 1,
      VERSION=='C|T' ~ 1,
      VERSION=='G/A' ~ 1,
      VERSION=='G|A' ~ 1,
      VERSION=='G/C' ~ 1,
      VERSION=='G|C' ~ 1,
      VERSION=='G/T' ~ 1,
      VERSION=='G|T' ~ 1,
      VERSION=='T/A' ~ 1,
      VERSION=='T|A' ~ 1,
      VERSION=='T/C' ~ 1,
      VERSION=='T|C' ~ 1,
      VERSION=='T/G' ~ 1,
      VERSION=='T|G' ~ 1))
  
  #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
  
  SNP_7_chrom1 =  SNP_6_chrom1 %>%
    group_by(GENE) %>%
    summarise(score = sum(score))
  
  SNP_GENES_2_chrom1 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_1")
  
  #Remettre la position des gènes et leur chromosome
  
  SNP_8_chrom1 = merge(x = SNP_7_chrom1, SNP_GENES_2,all.x=T,all.y=F)
  SNP_9_chrom1 = SNP_8_chrom1 %>% distinct()
  
  #Ajouter une colonne position et taille du gène
  
  SNP_10_chrom1 = SNP_9_chrom1 %>% mutate((SNP_9_chrom1$FIN+SNP_9_chrom1$DEBUT)/2)
  SNP_11_chrom1=  SNP_10_chrom1 %>% mutate((SNP_10_chrom1$FIN-SNP_10_chrom1$DEBUT)+1)
  head(SNP_11_chrom1)
  
  colnames(SNP_11_chrom1) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
  
  #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
  #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
  SNP_12_chrom1 = SNP_11_chrom1 %>% 
    mutate(SNPs_par_pb=SCORE/TAILLE)
  head(SNP_12_chrom1)
  
  #
  
  #CHROMOSOME 2 
  
  SNP_syn_chrom2 = SNPs_syn_Eff1 %>% 
    filter(CHROM=="CBS815.71p_contig_2")
  
  SNP_3_chrom2=SNP_2 %>% 
    filter(CHROM == "CBS815.71p_contig_2")
  
  SNP_GENES_chrom2 = SNP_GENES_1 %>% 
    filter(CHROM=="CBS815.71p_contig_2")
  
  SNP_4_chrom2=merge(x = SNP_3_chrom2, y = SNP_GENES_1, all.x=TRUE)
  
  #On sélectionne uniquement les SNPs synonymes
  
  SNP_4bis_chrom2 = merge(x = SNP_4_chrom2, y = SNP_syn_chrom2, all.x=FALSE, all.y=FALSE)
  
  SNP_4ter_chrom2= SNP_4bis_chrom2 %>% arrange(POS,GENE)
  
  SNP_5_chrom2=SNP_4bis_chrom2 %>% distinct(POS,.keep_all = TRUE)
  #On ajoute une colonne score
  SNP_6_chrom2=SNP_5_chrom2 %>%
    mutate(score = case_when(
      VERSION=="A/A" ~ 0,
      VERSION=='A|A' ~ 0,
      VERSION== "C/C" ~ 0,
      VERSION=='C|C' ~ 0,
      VERSION=='G/G' ~ 0,
      VERSION=='G|G' ~ 0,
      VERSION=='T/T' ~ 0,
      VERSION=='T|T' ~ 0,
      VERSION=='A/C' ~ 1,
      VERSION=='A|C' ~ 1,
      VERSION=='A/G' ~ 1,
      VERSION=='A|G' ~ 1,
      VERSION=='A/T' ~ 1,
      VERSION=='A|T' ~ 1,
      VERSION=='C/A' ~ 1,
      VERSION=='C|A' ~ 1,
      VERSION=='C/G' ~ 1,
      VERSION=='C|G' ~ 1,
      VERSION=='C/T' ~ 1,
      VERSION=='C|T' ~ 1,
      VERSION=='G/A' ~ 1,
      VERSION=='G|A' ~ 1,
      VERSION=='G/C' ~ 1,
      VERSION=='G|C' ~ 1,
      VERSION=='G/T' ~ 1,
      VERSION=='G|T' ~ 1,
      VERSION=='T/A' ~ 1,
      VERSION=='T|A' ~ 1,
      VERSION=='T/C' ~ 1,
      VERSION=='T|C' ~ 1,
      VERSION=='T/G' ~ 1,
      VERSION=='T|G' ~ 1))
  
  #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
  
  SNP_7_chrom2 =  SNP_6_chrom2 %>%
    group_by(GENE) %>%
    summarise(score = sum(score))
  
  SNP_GENES_2_chrom2 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_2")
  
  #Remettre la position des gènes et leur chromosome
  
  SNP_8_chrom2 = merge(x = SNP_7_chrom2, SNP_GENES_2,all.x=T,all.y=F)
  SNP_9_chrom2 = SNP_8_chrom2 %>% distinct()
  
  #Ajouter une colonne position et taille du gène
  
  SNP_10_chrom2 = SNP_9_chrom2 %>% mutate((SNP_9_chrom2$FIN+SNP_9_chrom2$DEBUT)/2)
  SNP_11_chrom2=  SNP_10_chrom2 %>% mutate((SNP_10_chrom2$FIN-SNP_10_chrom2$DEBUT)+1)
  head(SNP_11_chrom2)
  
  colnames(SNP_11_chrom2) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
  
  #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
  #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
  SNP_12_chrom2 = SNP_11_chrom2 %>% 
    mutate(SNPs_par_pb=SCORE/TAILLE)
  
  #Chromosome 3
  
  SNP_syn_chrom3 = SNPs_syn_Eff1 %>% 
    filter(CHROM=="CBS815.71p_contig_3")
  
  SNP_3_chrom3=SNP_2 %>% 
    filter(CHROM == "CBS815.71p_contig_3")
  
  SNP_GENES_chrom3 = SNP_GENES_1 %>% 
    filter(CHROM=="CBS815.71p_contig_3")
  
  SNP_4_chrom3=merge(x = SNP_3_chrom3, y = SNP_GENES_1, all.x=TRUE)
  
  #On sélectionne uniquement les SNPs synonymes
  
  SNP_4bis_chrom3 = merge(x = SNP_4_chrom3, y = SNP_syn_chrom3, all.x=FALSE, all.y=FALSE)
  
  SNP_4ter_chrom3= SNP_4bis_chrom3 %>% arrange(POS,GENE)
  
  SNP_5_chrom3=SNP_4bis_chrom3 %>% distinct(POS,.keep_all = TRUE)
  #On ajoute une colonne score
  SNP_6_chrom3=SNP_5_chrom3 %>%
    mutate(score = case_when(
      VERSION=="A/A" ~ 0,
      VERSION=='A|A' ~ 0,
      VERSION== "C/C" ~ 0,
      VERSION=='C|C' ~ 0,
      VERSION=='G/G' ~ 0,
      VERSION=='G|G' ~ 0,
      VERSION=='T/T' ~ 0,
      VERSION=='T|T' ~ 0,
      VERSION=='A/C' ~ 1,
      VERSION=='A|C' ~ 1,
      VERSION=='A/G' ~ 1,
      VERSION=='A|G' ~ 1,
      VERSION=='A/T' ~ 1,
      VERSION=='A|T' ~ 1,
      VERSION=='C/A' ~ 1,
      VERSION=='C|A' ~ 1,
      VERSION=='C/G' ~ 1,
      VERSION=='C|G' ~ 1,
      VERSION=='C/T' ~ 1,
      VERSION=='C|T' ~ 1,
      VERSION=='G/A' ~ 1,
      VERSION=='G|A' ~ 1,
      VERSION=='G/C' ~ 1,
      VERSION=='G|C' ~ 1,
      VERSION=='G/T' ~ 1,
      VERSION=='G|T' ~ 1,
      VERSION=='T/A' ~ 1,
      VERSION=='T|A' ~ 1,
      VERSION=='T/C' ~ 1,
      VERSION=='T|C' ~ 1,
      VERSION=='T/G' ~ 1,
      VERSION=='T|G' ~ 1))
  
  #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
  
  SNP_7_chrom3 =  SNP_6_chrom3 %>%
    group_by(GENE) %>%
    summarise(score = sum(score))
  
  SNP_GENES_2_chrom3 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_3")
  
  #Remettre la position des gènes et leur chromosome
  
  SNP_8_chrom3 = merge(x = SNP_7_chrom3, SNP_GENES_2,all.x=T,all.y=F)
  SNP_9_chrom3 = SNP_8_chrom3 %>% distinct()
  
  #Ajouter une colonne position et taille du gène
  
  SNP_10_chrom3 = SNP_9_chrom3 %>% mutate((SNP_9_chrom3$FIN+SNP_9_chrom3$DEBUT)/2)
  SNP_11_chrom3 =  SNP_10_chrom3 %>% mutate((SNP_10_chrom3$FIN-SNP_10_chrom3$DEBUT)+1)
  
  colnames(SNP_11_chrom3) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
  
  #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
  #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
  SNP_12_chrom3 = SNP_11_chrom3 %>% 
    mutate(SNPs_par_pb=SCORE/TAILLE)
  
  #CHROMOSOME 4
  SNP_syn_chrom4 = SNPs_syn_Eff1 %>% 
    filter(CHROM=="CBS815.71p_contig_4")
  
  SNP_3_chrom4=SNP_2 %>% 
    filter(CHROM == "CBS815.71p_contig_4")
  
  SNP_GENES_chrom4 = SNP_GENES_1 %>% 
    filter(CHROM=="CBS815.71p_contig_1")
  
  SNP_4_chrom4=merge(x = SNP_3_chrom4, y = SNP_GENES_1, all.x=TRUE)
  
  #On sélectionne uniquement les SNPs synonymes
  
  SNP_4bis_chrom4 = merge(x = SNP_4_chrom4, y = SNP_syn_chrom4, all.x=FALSE, all.y=FALSE)
  
  SNP_4ter_chrom4= SNP_4bis_chrom4 %>% arrange(POS,GENE)
  
  SNP_5_chrom4=SNP_4bis_chrom4 %>% distinct(POS,.keep_all = TRUE)
  #On ajoute une colonne score
  SNP_6_chrom4=SNP_5_chrom4 %>%
    mutate(score = case_when(
      VERSION=="A/A" ~ 0,
      VERSION=='A|A' ~ 0,
      VERSION== "C/C" ~ 0,
      VERSION=='C|C' ~ 0,
      VERSION=='G/G' ~ 0,
      VERSION=='G|G' ~ 0,
      VERSION=='T/T' ~ 0,
      VERSION=='T|T' ~ 0,
      VERSION=='A/C' ~ 1,
      VERSION=='A|C' ~ 1,
      VERSION=='A/G' ~ 1,
      VERSION=='A|G' ~ 1,
      VERSION=='A/T' ~ 1,
      VERSION=='A|T' ~ 1,
      VERSION=='C/A' ~ 1,
      VERSION=='C|A' ~ 1,
      VERSION=='C/G' ~ 1,
      VERSION=='C|G' ~ 1,
      VERSION=='C/T' ~ 1,
      VERSION=='C|T' ~ 1,
      VERSION=='G/A' ~ 1,
      VERSION=='G|A' ~ 1,
      VERSION=='G/C' ~ 1,
      VERSION=='G|C' ~ 1,
      VERSION=='G/T' ~ 1,
      VERSION=='G|T' ~ 1,
      VERSION=='T/A' ~ 1,
      VERSION=='T|A' ~ 1,
      VERSION=='T/C' ~ 1,
      VERSION=='T|C' ~ 1,
      VERSION=='T/G' ~ 1,
      VERSION=='T|G' ~ 1))
  
  #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
  
  SNP_7_chrom4 =  SNP_6_chrom4 %>%
    group_by(GENE) %>%
    summarise(score = sum(score))
  
  SNP_GENES_2_chrom4 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_4")
  
  #Remettre la position des gènes et leur chromosome
  
  SNP_8_chrom4 = merge(x = SNP_7_chrom4, SNP_GENES_2,all.x=T,all.y=F)
  SNP_9_chrom4 = SNP_8_chrom4 %>% distinct()
  
  #Ajouter une colonne position et taille du gène
  
  SNP_10_chrom4 = SNP_9_chrom4 %>% mutate((SNP_9_chrom4$FIN+SNP_9_chrom4$DEBUT)/2)
  SNP_11_chrom4=  SNP_10_chrom4 %>% mutate((SNP_10_chrom4$FIN-SNP_10_chrom4$DEBUT)+1)
  head(SNP_11_chrom4)
  
  colnames(SNP_11_chrom4) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
  
  #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
  #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
  SNP_12_chrom4 = SNP_11_chrom4 %>% 
    mutate(SNPs_par_pb=SCORE/TAILLE)
  
  #CHROMOSOME 5
  SNP_syn_chrom5 = SNPs_syn_Eff1 %>% 
    filter(CHROM=="CBS815.71p_contig_5")
  
  SNP_3_chrom5=SNP_2 %>% 
    filter(CHROM == "CBS815.71p_contig_5")
  
  SNP_GENES_chrom5 = SNP_GENES_1 %>% 
    filter(CHROM=="CBS815.71p_contig_5")
  
  SNP_4_chrom5=merge(x = SNP_3_chrom5, y = SNP_GENES_1, all.x=TRUE)
  
  #On sélectionne uniquement les SNPs synonymes
  
  SNP_4bis_chrom5 = merge(x = SNP_4_chrom5, y = SNP_syn_chrom5, all.x=FALSE, all.y=FALSE)
  
  SNP_4ter_chrom5 = SNP_4bis_chrom5 %>% arrange(POS,GENE)
  
  SNP_5_chrom5=SNP_4bis_chrom5 %>% distinct(POS,.keep_all = TRUE)
  #On ajoute une colonne score
  SNP_6_chrom5=SNP_5_chrom5 %>%
    mutate(score = case_when(
      VERSION=="A/A" ~ 0,
      VERSION=='A|A' ~ 0,
      VERSION== "C/C" ~ 0,
      VERSION=='C|C' ~ 0,
      VERSION=='G/G' ~ 0,
      VERSION=='G|G' ~ 0,
      VERSION=='T/T' ~ 0,
      VERSION=='T|T' ~ 0,
      VERSION=='A/C' ~ 1,
      VERSION=='A|C' ~ 1,
      VERSION=='A/G' ~ 1,
      VERSION=='A|G' ~ 1,
      VERSION=='A/T' ~ 1,
      VERSION=='A|T' ~ 1,
      VERSION=='C/A' ~ 1,
      VERSION=='C|A' ~ 1,
      VERSION=='C/G' ~ 1,
      VERSION=='C|G' ~ 1,
      VERSION=='C/T' ~ 1,
      VERSION=='C|T' ~ 1,
      VERSION=='G/A' ~ 1,
      VERSION=='G|A' ~ 1,
      VERSION=='G/C' ~ 1,
      VERSION=='G|C' ~ 1,
      VERSION=='G/T' ~ 1,
      VERSION=='G|T' ~ 1,
      VERSION=='T/A' ~ 1,
      VERSION=='T|A' ~ 1,
      VERSION=='T/C' ~ 1,
      VERSION=='T|C' ~ 1,
      VERSION=='T/G' ~ 1,
      VERSION=='T|G' ~ 1))
  
  #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
  
  SNP_7_chrom5 =  SNP_6_chrom5 %>%
    group_by(GENE) %>%
    summarise(score = sum(score))
  
  SNP_GENES_2_chrom5 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_5")
  
  #Remettre la position des gènes et leur chromosome
  
  SNP_8_chrom5 = merge(x = SNP_7_chrom5, SNP_GENES_2,all.x=T,all.y=F)
  SNP_9_chrom5 = SNP_8_chrom5 %>% distinct()
  
  #Ajouter une colonne position et taille du gène
  
  SNP_10_chrom5 = SNP_9_chrom5 %>% mutate((SNP_9_chrom5$FIN+SNP_9_chrom5$DEBUT)/2)
  SNP_11_chrom5=  SNP_10_chrom5 %>% mutate((SNP_10_chrom5$FIN-SNP_10_chrom5$DEBUT)+1)
  
  colnames(SNP_11_chrom5) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
  
  #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
  #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
  SNP_12_chrom5 = SNP_11_chrom5 %>% 
    mutate(SNPs_par_pb=SCORE/TAILLE)
  
  #Chromosome 6
  SNP_syn_chrom6 = SNPs_syn_Eff1 %>% 
    filter(CHROM=="CBS815.71p_contig_6")
  
  SNP_3_chrom6=SNP_2 %>% 
    filter(CHROM == "CBS815.71p_contig_6")
  
  SNP_GENES_chrom6 = SNP_GENES_1 %>% 
    filter(CHROM=="CBS815.71p_contig_6")
  
  SNP_4_chrom6=merge(x = SNP_3_chrom6, y = SNP_GENES_1, all.x=TRUE)
  
  #On sélectionne uniquement les SNPs synonymes
  
  SNP_4bis_chrom6 = merge(x = SNP_4_chrom6, y = SNP_syn_chrom6, all.x=FALSE, all.y=FALSE)
  
  SNP_4ter_chrom6 = SNP_4bis_chrom6 %>% arrange(POS,GENE)
  
  SNP_5_chrom6 = SNP_4bis_chrom6 %>% distinct(POS,.keep_all = TRUE)
  
  #On ajoute une colonne score
  SNP_6_chrom6=SNP_5_chrom6 %>%
    mutate(score = case_when(
      VERSION=="A/A" ~ 0,
      VERSION=='A|A' ~ 0,
      VERSION== "C/C" ~ 0,
      VERSION=='C|C' ~ 0,
      VERSION=='G/G' ~ 0,
      VERSION=='G|G' ~ 0,
      VERSION=='T/T' ~ 0,
      VERSION=='T|T' ~ 0,
      VERSION=='A/C' ~ 1,
      VERSION=='A|C' ~ 1,
      VERSION=='A/G' ~ 1,
      VERSION=='A|G' ~ 1,
      VERSION=='A/T' ~ 1,
      VERSION=='A|T' ~ 1,
      VERSION=='C/A' ~ 1,
      VERSION=='C|A' ~ 1,
      VERSION=='C/G' ~ 1,
      VERSION=='C|G' ~ 1,
      VERSION=='C/T' ~ 1,
      VERSION=='C|T' ~ 1,
      VERSION=='G/A' ~ 1,
      VERSION=='G|A' ~ 1,
      VERSION=='G/C' ~ 1,
      VERSION=='G|C' ~ 1,
      VERSION=='G/T' ~ 1,
      VERSION=='G|T' ~ 1,
      VERSION=='T/A' ~ 1,
      VERSION=='T|A' ~ 1,
      VERSION=='T/C' ~ 1,
      VERSION=='T|C' ~ 1,
      VERSION=='T/G' ~ 1,
      VERSION=='T|G' ~ 1))
  
  #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
  
  SNP_7_chrom6 =  SNP_6_chrom6 %>%
    group_by(GENE) %>%
    summarise(score = sum(score))
  
  SNP_GENES_2_chrom6 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_6")
  
  #Remettre la position des gènes et leur chromosome
  
  SNP_8_chrom6 = merge(x = SNP_7_chrom6, SNP_GENES_2,all.x=T,all.y=F)
  SNP_9_chrom6 = SNP_8_chrom6 %>% distinct()
  
  #Ajouter une colonne position et taille du gène
  
  SNP_10_chrom6 = SNP_9_chrom6 %>% mutate((SNP_9_chrom6$FIN+SNP_9_chrom6$DEBUT)/2)
  SNP_11_chrom6=  SNP_10_chrom6 %>% mutate((SNP_10_chrom6$FIN-SNP_10_chrom6$DEBUT)+1)
  
  colnames(SNP_11_chrom6) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
  
  #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
  #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
  
  SNP_12_chrom6 = SNP_11_chrom6 %>% 
    mutate(SNPs_par_pb=SCORE/TAILLE)
  
  #Chromosome 7
  SNP_syn_chrom7 = SNPs_syn_Eff1 %>% 
    filter(CHROM=="CBS815.71p_contig_7")
  
  SNP_3_chrom7=SNP_2 %>% 
    filter(CHROM == "CBS815.71p_contig_7")
  
  SNP_GENES_chrom7 = SNP_GENES_1 %>% 
    filter(CHROM=="CBS815.71p_contig_7")
  
  SNP_4_chrom7=merge(x = SNP_3_chrom7, y = SNP_GENES_1, all.x=TRUE)
  
  #On sélectionne uniquement les SNPs synonymes
  
  SNP_4bis_chrom7 = merge(x = SNP_4_chrom7, y = SNP_syn_chrom7, all.x=FALSE, all.y=FALSE)
  
  SNP_4ter_chrom7= SNP_4bis_chrom7 %>% arrange(POS,GENE)
  
  SNP_5_chrom7=SNP_4bis_chrom7 %>% distinct(POS,.keep_all = TRUE)
  #On ajoute une colonne score
  SNP_6_chrom7=SNP_5_chrom7 %>%
    mutate(score = case_when(
      VERSION=="A/A" ~ 0,
      VERSION=='A|A' ~ 0,
      VERSION== "C/C" ~ 0,
      VERSION=='C|C' ~ 0,
      VERSION=='G/G' ~ 0,
      VERSION=='G|G' ~ 0,
      VERSION=='T/T' ~ 0,
      VERSION=='T|T' ~ 0,
      VERSION=='A/C' ~ 1,
      VERSION=='A|C' ~ 1,
      VERSION=='A/G' ~ 1,
      VERSION=='A|G' ~ 1,
      VERSION=='A/T' ~ 1,
      VERSION=='A|T' ~ 1,
      VERSION=='C/A' ~ 1,
      VERSION=='C|A' ~ 1,
      VERSION=='C/G' ~ 1,
      VERSION=='C|G' ~ 1,
      VERSION=='C/T' ~ 1,
      VERSION=='C|T' ~ 1,
      VERSION=='G/A' ~ 1,
      VERSION=='G|A' ~ 1,
      VERSION=='G/C' ~ 1,
      VERSION=='G|C' ~ 1,
      VERSION=='G/T' ~ 1,
      VERSION=='G|T' ~ 1,
      VERSION=='T/A' ~ 1,
      VERSION=='T|A' ~ 1,
      VERSION=='T/C' ~ 1,
      VERSION=='T|C' ~ 1,
      VERSION=='T/G' ~ 1,
      VERSION=='T|G' ~ 1))
  
  #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
  
  SNP_7_chrom7 =  SNP_6_chrom7 %>%
    group_by(GENE) %>%
    summarise(score = sum(score))
  
  SNP_GENES_2_chrom7 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_7")
  
  #Remettre la position des gènes et leur chromosome
  
  SNP_8_chrom7 = merge(x = SNP_7_chrom7, SNP_GENES_2,all.x=T,all.y=F)
  SNP_9_chrom7 = SNP_8_chrom7 %>% distinct()
  
  #Ajouter une colonne position et taille du gène
  
  SNP_10_chrom7 = SNP_9_chrom7 %>% mutate((SNP_9_chrom7$FIN+SNP_9_chrom7$DEBUT)/2)
  SNP_11_chrom7=  SNP_10_chrom7 %>% mutate((SNP_10_chrom7$FIN-SNP_10_chrom7$DEBUT)+1)
  head(SNP_11_chrom7)
  
  colnames(SNP_11_chrom7) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
  
  #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
  #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
  SNP_12_chrom7 = SNP_11_chrom7 %>% 
    mutate(SNPs_par_pb=SCORE/TAILLE)
  
  #On fait un seul tableau en concaténant tous les tableaux
  
  SNP_13 <- bind_rows(SNP_12_chrom1, SNP_12_chrom2, SNP_12_chrom3, SNP_12_chrom4, SNP_12_chrom5,SNP_12_chrom6,SNP_12_chrom7)
  
  #On représente l'hétérozygotie (Nombre de SNPs hétérozygote par gènes par paire de base)
  #sur l'ensemble de la souche 
  
  plot<-ggplot(SNP_13, aes(y=SNPs_par_pb , x=POS_GENE/1000)) +
    labs(x=" ", y="SNPs synonymes par gène moyennée par pb") +
    geom_point(aes(color=variable),size=0.5,alpha=0.8,colour = "black") + 
    facet_grid( CHROM ~ ., scale="free")+ theme(legend.position="none")
  plot 
  
  #On ajouter une variable nom de la souche
  
  SNP_14 = SNP_13 %>% 
    mutate(SOUCHE=i)
  
  #On enregistre une table avec uniquement le premier chromosome
  SNP_14_chrom1 = SNP_12_chrom1 %>% 
    mutate(SOUCHE=i)
  
  write.table(SNP_14_chrom1, file=paste0("VERIF_CORRECTION_souche_SNPsSyn_avecHET",i,".txt"),row.names=T,col.names=T)
  
  #on enregistre une table avec tous les chromosomes
  write.table(SNP_14, file=paste0("VERIF_CORRECTION_souche_ENTIERE_SNPsSyn_avecHET",i,".txt"),row.names=T,col.names=T)
  
}

##Boucle homokaryons
#Boucle homokaryons###########

#On reprend le script en créant une boucle pour toutes les souches

-
  homokaryons = c("S1.PSN881","S5.PSQ32","S1.CBS81571","S13.PSN1042","S3.PSN1070","S2.PSN589","S8.PSN644","S12.PSN663","S3.PSN705","S9.PSN707","S1.PSN775","S6.PSN831","S15.PSN377")
  
  for (i in homokaryons ) {
    
    SNP_GENES<- read.table("/Volumes/Elsa_Donnees/Octobre_2022/MAT+/VCF_final/VCF_final.strates_avecHET.Bedtools_intersect_SNPs_genes.bis.bed", header=F)
    colnames(SNP_GENES) <- c("CHROM", "POS", "VA","VB","DEBUT","FIN","VC","GENE","VD","VE")
    SNP_GENES_1=select(SNP_GENES, -c(3,4,7,9,10))
    SNP_GENES_2=select(SNP_GENES_1,-c(2))
    
    #Ouvrir les fichier de SNPs : fichiers avec tous les SNPs  
    
    SNP_m <- read.table(paste0("/Volumes/Elsa_Donnees/Octobre_2022/MAT+/VCF_final/souches_avecHET/Genotype.dikaryotic.Stetrasporum.",i,".m.VCF_final.Genes_avecHET.sansMAAF.txt"),header=T)
    SNP_p <- read.table(paste0("/Volumes/Elsa_Donnees/Octobre_2022/MAT+/VCF_final/souches_avecHET/Genotype.dikaryotic.Stetrasporum.",i,".p.VCF_final.Genes_avecHET.sansMAAF.txt"),header=T)
    colnames(SNP_m) <- c("CHROM", "POS", "TYPE", "REF", "VERSION_m")
    colnames(SNP_p) <- c("CHROM", "POS", "TYPE", "REF", "VERSION_p")
    
    head(SNP_m, 30)
    head(SNP_p, 30)
    
    SNPs_syn<- read.table("/Volumes/Elsa_Donnees/Octobre_2022/MAT+/VCF_final/VCF.Octobre2022.CBS815.71Sp3.final.diploid.VCF_final.genotyped.1.SNPs.filtered.geno0.9.only.biallelic.pol.renamed.strates_avecHET.genes.recode.annotated.effects.txt", header=F,col.names = paste0("V",seq_len(21)), fill = TRUE)
    colnames(SNPs_syn)=c("CHROM","POS","REF","ALT","EFFECT_1","EFFECT_2","EFFECT_3","EFFECT_4","EFFECT_5","EFFECT_6","EFFECT_7","EFFECT_8","EFFECT_9","EFFECT_10","EFFECT_11","EFFECT_12","EFFECT_13","EFFECT_14","EFFECT_15","EFFECT_16","EFFECT_17")
    head(SNPs_syn)
    SNPs_syn_1 = SNPs_syn %>% slice(-1)
    
    SNPs_syn_Eff1 = SNPs_syn_1  %>% filter(EFFECT_1=="synonymous_variant")
    
    #D'abord fusionner les fichiers
    
    #Filtrer les SNPS non déterminés
    SNP_1_m=filter(SNP_m,  VERSION_m!= ".")
    head(SNP_1_m,30)
    
    SNP_2_m=filter(SNP_1_m,  VERSION_m!= "./.")
    head(SNP_2_m,30)
    
    SNP_1_p=filter(SNP_p,  VERSION_p!= ".")
    head(SNP_1_p,30)
    
    SNP_2_p=filter(SNP_1_p,  VERSION_p!= "./.")
    head(SNP_2_p,30)
    
    #D'abord fusionner les fichiers
    
    SNP_2=merge(SNP_2_p,SNP_2_m,all.x=T,all.y=T)
    head(SNP_2)
    
    SNP_3= SNP_2 %>%
      mutate(SNP = paste0(VERSION_m,"/",VERSION_p))
    head(SNP_3)
    
    #Chromosome par chromosome
    
    ###CHROMOSOME 1 
    SNP_syn_chrom1 = SNPs_syn_Eff1 %>% 
      filter(CHROM=="CBS815.71p_contig_1")
    
    SNP_4_chrom1=SNP_3 %>% 
      filter(CHROM == "CBS815.71p_contig_1")
    
    #Ordonner en croissant par positions puis par gène quand 2 fois la même position pour que 
    #quand deux fois le même SNP dans deux gènes : celui du "premier" gène dans le chromosome au niveau de la position seulement soit gardé
    #Permet d'avoir un critère répétable
    
    SNP_5_chrom1=merge(x = SNP_4_chrom1, y = SNP_GENES_1, all.x=TRUE)
    
    SNP_5bis_chrom1 = merge(x = SNP_5_chrom1, y = SNP_syn_chrom1, all.x=FALSE, all.y=FALSE)
    
    SNP_6_chrom1= SNP_5bis_chrom1 %>% arrange(POS,GENE)
    SNP_7_chrom1=SNP_6_chrom1 %>% distinct(POS,.keep_all = TRUE)
    
    #On ajouter une colonne score 
    SNP_8_chrom1=SNP_7_chrom1 %>%
      mutate(score = case_when(
        SNP=="A|A/A|A" ~ 0,
        SNP=="A/A/A/A" ~ 0,
        SNP=="A/A/A|A" ~ 0,
        SNP=="A|A/A/A" ~ 0,
        SNP=="T|T/T|T" ~ 0,
        SNP=="T/T/T/T" ~ 0,
        SNP=="T/T/T|T" ~ 0,
        SNP=="T|T/T/T" ~ 0,
        SNP=="C|C/C|C" ~ 0,
        SNP=="C/C/C/C" ~ 0,
        SNP=="C/C/C|C" ~ 0,
        SNP=="C|C/C/C" ~ 0,
        SNP=="G|G/G|G" ~ 0,
        SNP=="G/G/G/G" ~ 0,
        SNP=="G/G/G|G" ~ 0,
        SNP=="G|G/G/G" ~ 0,
        SNP=='A|A/C|C' ~ 1,
        SNP=='A/A/C|C' ~ 1,
        SNP=='A|A/C/C' ~ 1,
        SNP=='A/A/C/C' ~ 1,
        SNP=="A|A/G|G" ~ 1,
        SNP=="A/A/G/G" ~ 1,
        SNP=="A/A/G|G" ~ 1,
        SNP=="A|A/G/G" ~ 1,
        SNP=="A|A/T/T" ~ 1,
        SNP=="A|A/T|T" ~ 1,
        SNP=="A/A/T/T" ~ 1,
        SNP=="A/A/T|T" ~ 1,
        SNP=='C|C/A|A' ~ 1,
        SNP=='C/C/A|A' ~ 1,
        SNP=='C|C/A/A' ~ 1,
        SNP=='C/C/A/A' ~ 1,
        SNP=='C|C/G|G' ~ 1,
        SNP=='C/C/G|G' ~ 1,
        SNP=='C|C/G/G' ~ 1,
        SNP=='C/C/G/G' ~ 1,
        SNP=='C|C/T|T' ~ 1,
        SNP=='C/C/T|T' ~ 1,
        SNP=='C|C/T/T' ~ 1,
        SNP=='C/C/T/T' ~ 1,
        SNP=='G|G/A|A' ~ 1,
        SNP=='G/G/A|A' ~ 1,
        SNP=='G|G/A/A' ~ 1,
        SNP=='G/G/A/A' ~ 1,
        SNP=='G|G/C|C' ~ 1,
        SNP=='G/G/C|C' ~ 1,
        SNP=='G|G/C/C' ~ 1,
        SNP=='G/G/C/C' ~ 1,
        SNP=='G|G/T|T' ~ 1,
        SNP=='G/G/T|T' ~ 1,
        SNP=='G|G/T/T' ~ 1,
        SNP=='G/G/T/T' ~ 1,
        SNP=='T|T/A|A' ~ 1,
        SNP=='T/T/A|A' ~ 1,
        SNP=='T|T/A/A' ~ 1,
        SNP=='T/T/A/A' ~ 1,
        SNP=='T|T/C|C' ~ 1,
        SNP=='T/T/C|C' ~ 1,
        SNP=='T|T/C/C' ~ 1,
        SNP=='T/T/C/C' ~ 1,
        SNP=='T|T/G|G' ~ 1,
        SNP=='T/T/G|G' ~ 1,
        SNP=='T|T/G/G' ~ 1,
        SNP=='T/T/G/G' ~ 1))
    
    SNP_9_chrom1 = filter(SNP_8_chrom1, SNP != "A|A/NA")
    SNP_10_chrom1 = filter(SNP_9_chrom1, SNP != "A/A/NA")
    SNP_11_chrom1 = filter(SNP_10_chrom1, SNP != "NA/A/A")
    SNP_12_chrom1 = filter(SNP_11_chrom1, SNP != "NA/A|A")
    
    SNP_13_chrom1 = filter(SNP_12_chrom1, SNP != "C|C/NA")
    SNP_14_chrom1 = filter(SNP_13_chrom1, SNP != "C/C/NA")
    SNP_15_chrom1 = filter(SNP_14_chrom1, SNP != "NA/C/C")
    SNP_16_chrom1 = filter(SNP_15_chrom1, SNP != "NA/C|C")
    
    SNP_17_chrom1 = filter(SNP_16_chrom1, SNP != "T|T/NA")
    SNP_18_chrom1 = filter(SNP_17_chrom1, SNP != "T/T/NA")
    SNP_19_chrom1 = filter(SNP_18_chrom1, SNP != "NA/T/T")
    SNP_20_chrom1 = filter(SNP_19_chrom1, SNP != "NA/T|T")
    
    SNP_21_chrom1 = filter(SNP_20_chrom1, SNP != "G|G/NA")
    SNP_22_chrom1 = filter(SNP_21_chrom1, SNP != "G/G/NA")
    SNP_23_chrom1 = filter(SNP_22_chrom1, SNP != "NA/G/G")
    SNP_24_chrom1 = filter(SNP_23_chrom1, SNP != "NA/G|G")
    
    head(SNP_24_chrom1,200)
    
    #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
    SNP_24bis_chrom1 = na.omit(SNP_24_chrom1)
    head(SNP_24_chrom1,200)
    
    
    SNP_25_chrom1 =  SNP_24bis_chrom1 %>%
      group_by(GENE) %>%
      summarise(score = sum(score))
    
    SNP_GENES_2_chrom1 = filter(SNP_GENES_2,CHROM == "CBS815.71p_contig_1")
    
    
    #Remettre la position des gènes et leur chromosome
    
    SNP_26_chrom1 = merge(x = SNP_25_chrom1, SNP_GENES_2,all.x=T,all.y=F)
    SNP_27_chrom1 = SNP_26_chrom1 %>% distinct()
    
    #Ajouter une colonne position et taille du gène
    
    SNP_28_chrom1 = SNP_27_chrom1 %>% mutate((SNP_27_chrom1$FIN+SNP_27_chrom1$DEBUT)/2)
    SNP_29_chrom1=  SNP_28_chrom1 %>% mutate((SNP_28_chrom1$FIN-SNP_28_chrom1$DEBUT)+1)
    head(SNP_29_chrom1)
    colnames(SNP_29_chrom1) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
    head(SNP_29_chrom1)
    
    #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
    #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
    SNP_30_chrom1 = SNP_29_chrom1 %>% 
      mutate(SNPs_par_pb=SCORE/TAILLE)
    head(SNP_30_chrom1)
    
    ###CHROMOSOME 2
    SNP_syn_chrom2 = SNPs_syn_Eff1 %>% 
      filter(CHROM=="CBS815.71p_contig_2")
    
    SNP_4_chrom2=SNP_3 %>% 
      filter(CHROM == "CBS815.71p_contig_2")
    
    #Ordonner en croissant par positions puis par gène quand 2 fois la même position pour que 
    #quand deux fois le même SNP dans deux gènes : celui du "premier" gène dans le chromosome au niveau de la position seulement soit gardé
    #Permet d'avoir un critère répétable
    
    SNP_5_chrom2=merge(x = SNP_4_chrom2, y = SNP_GENES_1, all.x=TRUE)
    
    SNP_5bis_chrom2 = merge(x = SNP_5_chrom2, y = SNP_syn_chrom2, all.x=FALSE, all.y=FALSE)
    
    SNP_6_chrom2= SNP_5bis_chrom2 %>% arrange(POS,GENE)
    
    SNP_7_chrom2=SNP_6_chrom2 %>% distinct(POS,.keep_all = TRUE)
    
    #On ajouter une colonne score 
    SNP_8_chrom2=SNP_7_chrom2 %>%
      mutate(score = case_when(
        SNP=="A|A/A|A" ~ 0,
        SNP=="A/A/A/A" ~ 0,
        SNP=="A/A/A|A" ~ 0,
        SNP=="A|A/A/A" ~ 0,
        SNP=="T|T/T|T" ~ 0,
        SNP=="T/T/T/T" ~ 0,
        SNP=="T/T/T|T" ~ 0,
        SNP=="T|T/T/T" ~ 0,
        SNP=="C|C/C|C" ~ 0,
        SNP=="C/C/C/C" ~ 0,
        SNP=="C/C/C|C" ~ 0,
        SNP=="C|C/C/C" ~ 0,
        SNP=="G|G/G|G" ~ 0,
        SNP=="G/G/G/G" ~ 0,
        SNP=="G/G/G|G" ~ 0,
        SNP=="G|G/G/G" ~ 0,
        SNP=='A|A/C|C' ~ 1,
        SNP=='A/A/C|C' ~ 1,
        SNP=='A|A/C/C' ~ 1,
        SNP=='A/A/C/C' ~ 1,
        SNP=="A|A/G|G" ~ 1,
        SNP=="A/A/G/G" ~ 1,
        SNP=="A/A/G|G" ~ 1,
        SNP=="A|A/G/G" ~ 1,
        SNP=="A|A/T/T" ~ 1,
        SNP=="A|A/T|T" ~ 1,
        SNP=="A/A/T/T" ~ 1,
        SNP=="A/A/T|T" ~ 1,
        SNP=='C|C/A|A' ~ 1,
        SNP=='C/C/A|A' ~ 1,
        SNP=='C|C/A/A' ~ 1,
        SNP=='C/C/A/A' ~ 1,
        SNP=='C|C/G|G' ~ 1,
        SNP=='C/C/G|G' ~ 1,
        SNP=='C|C/G/G' ~ 1,
        SNP=='C/C/G/G' ~ 1,
        SNP=='C|C/T|T' ~ 1,
        SNP=='C/C/T|T' ~ 1,
        SNP=='C|C/T/T' ~ 1,
        SNP=='C/C/T/T' ~ 1,
        SNP=='G|G/A|A' ~ 1,
        SNP=='G/G/A|A' ~ 1,
        SNP=='G|G/A/A' ~ 1,
        SNP=='G/G/A/A' ~ 1,
        SNP=='G|G/C|C' ~ 1,
        SNP=='G/G/C|C' ~ 1,
        SNP=='G|G/C/C' ~ 1,
        SNP=='G/G/C/C' ~ 1,
        SNP=='G|G/T|T' ~ 1,
        SNP=='G/G/T|T' ~ 1,
        SNP=='G|G/T/T' ~ 1,
        SNP=='G/G/T/T' ~ 1,
        SNP=='T|T/A|A' ~ 1,
        SNP=='T/T/A|A' ~ 1,
        SNP=='T|T/A/A' ~ 1,
        SNP=='T/T/A/A' ~ 1,
        SNP=='T|T/C|C' ~ 1,
        SNP=='T/T/C|C' ~ 1,
        SNP=='T|T/C/C' ~ 1,
        SNP=='T/T/C/C' ~ 1,
        SNP=='T|T/G|G' ~ 1,
        SNP=='T/T/G|G' ~ 1,
        SNP=='T|T/G/G' ~ 1,
        SNP=='T/T/G/G' ~ 1))
    
    SNP_9_chrom2 = filter(SNP_8_chrom2, SNP != "A|A/NA")
    SNP_10_chrom2 = filter(SNP_9_chrom2, SNP != "A/A/NA")
    SNP_11_chrom2 = filter(SNP_10_chrom2, SNP != "NA/A/A")
    SNP_12_chrom2 = filter(SNP_11_chrom2, SNP != "NA/A|A")
    
    SNP_13_chrom2 = filter(SNP_12_chrom2, SNP != "C|C/NA")
    SNP_14_chrom2 = filter(SNP_13_chrom2, SNP != "C/C/NA")
    SNP_15_chrom2 = filter(SNP_14_chrom2, SNP != "NA/C/C")
    SNP_16_chrom2 = filter(SNP_15_chrom2, SNP != "NA/C|C")
    
    SNP_17_chrom2 = filter(SNP_16_chrom2, SNP != "T|T/NA")
    SNP_18_chrom2 = filter(SNP_17_chrom2, SNP != "T/T/NA")
    SNP_19_chrom2 = filter(SNP_18_chrom2, SNP != "NA/T/T")
    SNP_20_chrom2 = filter(SNP_19_chrom2, SNP != "NA/T|T")
    
    SNP_21_chrom2 = filter(SNP_20_chrom2, SNP != "G|G/NA")
    SNP_22_chrom2 = filter(SNP_21_chrom2, SNP != "G/G/NA")
    SNP_23_chrom2 = filter(SNP_22_chrom2, SNP != "NA/G/G")
    SNP_24_chrom2 = filter(SNP_23_chrom2, SNP != "NA/G|G")
    
    head(SNP_24_chrom2,200)
    
    SNP_24bis_chrom2 = na.omit(SNP_24_chrom2)
    head(SNP_24_chrom2,200)
    
    #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
    
    SNP_25_chrom2 =  SNP_24bis_chrom2 %>%
      group_by(GENE) %>%
      summarise(score = sum(score))
    #Remettre la position des gènes et leur chromosome
    
    SNP_26_chrom2 = merge(x = SNP_25_chrom2, SNP_GENES_2,all.x=T,all.y=F)
    SNP_27_chrom2 = SNP_26_chrom2 %>% distinct()
    
    #Ajouter une colonne position et taille du gène
    
    SNP_28_chrom2 = SNP_27_chrom2 %>% mutate((SNP_27_chrom2$FIN+SNP_27_chrom2$DEBUT)/2)
    SNP_29_chrom2=  SNP_28_chrom2 %>% mutate((SNP_28_chrom2$FIN-SNP_28_chrom2$DEBUT)+1)
    head(SNP_29_chrom2)
    colnames(SNP_29_chrom2) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
    head(SNP_29_chrom2)
    
    #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
    #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
    SNP_30_chrom2 = SNP_29_chrom2 %>% 
      mutate(SNPs_par_pb=SCORE/TAILLE)
    head(SNP_30_chrom2)
    
    ####CHROMOSOME 3
    SNP_syn_chrom3 = SNPs_syn_Eff1 %>% 
      filter(CHROM=="CBS815.71p_contig_3")
    
    SNP_4_chrom3=SNP_3 %>% 
      filter(CHROM == "CBS815.71p_contig_3")
    
    #Ordonner en croissant par positions puis par gène quand 2 fois la même position pour que 
    #quand deux fois le même SNP dans deux gènes : celui du "premier" gène dans le chromosome au niveau de la position seulement soit gardé
    #Permet d'avoir un critère répétable
    
    SNP_5_chrom3=merge(x = SNP_4_chrom3, y = SNP_GENES_1, all.x=TRUE)
    
    SNP_5bis_chrom3 = merge(x = SNP_5_chrom3, y = SNP_syn_chrom3, all.x=FALSE, all.y=FALSE)
    
    SNP_6_chrom3= SNP_5bis_chrom3 %>% arrange(POS,GENE)
    SNP_7_chrom3=SNP_6_chrom3 %>% distinct(POS,.keep_all = TRUE)
    
    #On ajouter une colonne score 
    SNP_8_chrom3=SNP_7_chrom3 %>%
      mutate(score = case_when(
        SNP=="A|A/A|A" ~ 0,
        SNP=="A/A/A/A" ~ 0,
        SNP=="A/A/A|A" ~ 0,
        SNP=="A|A/A/A" ~ 0,
        SNP=="T|T/T|T" ~ 0,
        SNP=="T/T/T/T" ~ 0,
        SNP=="T/T/T|T" ~ 0,
        SNP=="T|T/T/T" ~ 0,
        SNP=="C|C/C|C" ~ 0,
        SNP=="C/C/C/C" ~ 0,
        SNP=="C/C/C|C" ~ 0,
        SNP=="C|C/C/C" ~ 0,
        SNP=="G|G/G|G" ~ 0,
        SNP=="G/G/G/G" ~ 0,
        SNP=="G/G/G|G" ~ 0,
        SNP=="G|G/G/G" ~ 0,
        SNP=='A|A/C|C' ~ 1,
        SNP=='A/A/C|C' ~ 1,
        SNP=='A|A/C/C' ~ 1,
        SNP=='A/A/C/C' ~ 1,
        SNP=="A|A/G|G" ~ 1,
        SNP=="A/A/G/G" ~ 1,
        SNP=="A/A/G|G" ~ 1,
        SNP=="A|A/G/G" ~ 1,
        SNP=="A|A/T/T" ~ 1,
        SNP=="A|A/T|T" ~ 1,
        SNP=="A/A/T/T" ~ 1,
        SNP=="A/A/T|T" ~ 1,
        SNP=='C|C/A|A' ~ 1,
        SNP=='C/C/A|A' ~ 1,
        SNP=='C|C/A/A' ~ 1,
        SNP=='C/C/A/A' ~ 1,
        SNP=='C|C/G|G' ~ 1,
        SNP=='C/C/G|G' ~ 1,
        SNP=='C|C/G/G' ~ 1,
        SNP=='C/C/G/G' ~ 1,
        SNP=='C|C/T|T' ~ 1,
        SNP=='C/C/T|T' ~ 1,
        SNP=='C|C/T/T' ~ 1,
        SNP=='C/C/T/T' ~ 1,
        SNP=='G|G/A|A' ~ 1,
        SNP=='G/G/A|A' ~ 1,
        SNP=='G|G/A/A' ~ 1,
        SNP=='G/G/A/A' ~ 1,
        SNP=='G|G/C|C' ~ 1,
        SNP=='G/G/C|C' ~ 1,
        SNP=='G|G/C/C' ~ 1,
        SNP=='G/G/C/C' ~ 1,
        SNP=='G|G/T|T' ~ 1,
        SNP=='G/G/T|T' ~ 1,
        SNP=='G|G/T/T' ~ 1,
        SNP=='G/G/T/T' ~ 1,
        SNP=='T|T/A|A' ~ 1,
        SNP=='T/T/A|A' ~ 1,
        SNP=='T|T/A/A' ~ 1,
        SNP=='T/T/A/A' ~ 1,
        SNP=='T|T/C|C' ~ 1,
        SNP=='T/T/C|C' ~ 1,
        SNP=='T|T/C/C' ~ 1,
        SNP=='T/T/C/C' ~ 1,
        SNP=='T|T/G|G' ~ 1,
        SNP=='T/T/G|G' ~ 1,
        SNP=='T|T/G/G' ~ 1,
        SNP=='T/T/G/G' ~ 1))
    
    SNP_9_chrom3 = filter(SNP_8_chrom3, SNP != "A|A/NA")
    SNP_10_chrom3 = filter(SNP_9_chrom3, SNP != "A/A/NA")
    SNP_11_chrom3 = filter(SNP_10_chrom3, SNP != "NA/A/A")
    SNP_12_chrom3 = filter(SNP_11_chrom3, SNP != "NA/A|A")
    
    SNP_13_chrom3 = filter(SNP_12_chrom3, SNP != "C|C/NA")
    SNP_14_chrom3 = filter(SNP_13_chrom3, SNP != "C/C/NA")
    SNP_15_chrom3 = filter(SNP_14_chrom3, SNP != "NA/C/C")
    SNP_16_chrom3 = filter(SNP_15_chrom3, SNP != "NA/C|C")
    
    SNP_17_chrom3 = filter(SNP_16_chrom3, SNP != "T|T/NA")
    SNP_18_chrom3 = filter(SNP_17_chrom3, SNP != "T/T/NA")
    SNP_19_chrom3 = filter(SNP_18_chrom3, SNP != "NA/T/T")
    SNP_20_chrom3 = filter(SNP_19_chrom3, SNP != "NA/T|T")
    
    SNP_21_chrom3 = filter(SNP_20_chrom3, SNP != "G|G/NA")
    SNP_22_chrom3 = filter(SNP_21_chrom3, SNP != "G/G/NA")
    SNP_23_chrom3 = filter(SNP_22_chrom3, SNP != "NA/G/G")
    SNP_24_chrom3 = filter(SNP_23_chrom3, SNP != "NA/G|G")
    
    head(SNP_24_chrom3,200)
    
    SNP_24bis_chrom3 = na.omit(SNP_24_chrom3)
    head(SNP_24_chrom3,200)
    
    #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
    
    SNP_25_chrom3 =  SNP_24bis_chrom3 %>%
      group_by(GENE) %>%
      summarise(score = sum(score))
    #Remettre la position des gènes et leur chromosome
    
    SNP_26_chrom3 = merge(x = SNP_25_chrom3, SNP_GENES_2,all.x=T,all.y=F)
    SNP_27_chrom3 = SNP_26_chrom3 %>% distinct()
    
    #Ajouter une colonne position et taille du gène
    
    SNP_28_chrom3 = SNP_27_chrom3 %>% mutate((SNP_27_chrom3$FIN+SNP_27_chrom3$DEBUT)/2)
    SNP_29_chrom3=  SNP_28_chrom3 %>% mutate((SNP_28_chrom3$FIN-SNP_28_chrom3$DEBUT)+1)
    head(SNP_29_chrom3)
    colnames(SNP_29_chrom3) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
    head(SNP_29_chrom3)
    
    #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
    #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
    SNP_30_chrom3 = SNP_29_chrom3 %>% 
      mutate(SNPs_par_pb=SCORE/TAILLE)
    head(SNP_30_chrom3)
    
    ###CHROMOSOME 4
    SNP_syn_chrom4 = SNPs_syn_Eff1 %>% 
      filter(CHROM=="CBS815.71p_contig_4")
    
    SNP_4_chrom4=SNP_3 %>% 
      filter(CHROM == "CBS815.71p_contig_4")
    
    #Ordonner en croissant par positions puis par gène quand 2 fois la même position pour que 
    #quand deux fois le même SNP dans deux gènes : celui du "premier" gène dans le chromosome au niveau de la position seulement soit gardé
    #Permet d'avoir un critère répétable
    
    SNP_5_chrom4=merge(x = SNP_4_chrom4, y = SNP_GENES_1, all.x=TRUE)
    
    SNP_5bis_chrom4 = merge(x = SNP_5_chrom4, y = SNP_syn_chrom4, all.x=FALSE, all.y=FALSE)
    
    SNP_6_chrom4 = SNP_5bis_chrom4 %>% arrange(POS,GENE)
    
    SNP_7_chrom4=SNP_6_chrom4 %>% distinct(POS,.keep_all = TRUE)
    
    #On ajouter une colonne score 
    SNP_8_chrom4=SNP_7_chrom4 %>%
      mutate(score = case_when(
        SNP=="A|A/A|A" ~ 0,
        SNP=="A/A/A/A" ~ 0,
        SNP=="A/A/A|A" ~ 0,
        SNP=="A|A/A/A" ~ 0,
        SNP=="T|T/T|T" ~ 0,
        SNP=="T/T/T/T" ~ 0,
        SNP=="T/T/T|T" ~ 0,
        SNP=="T|T/T/T" ~ 0,
        SNP=="C|C/C|C" ~ 0,
        SNP=="C/C/C/C" ~ 0,
        SNP=="C/C/C|C" ~ 0,
        SNP=="C|C/C/C" ~ 0,
        SNP=="G|G/G|G" ~ 0,
        SNP=="G/G/G/G" ~ 0,
        SNP=="G/G/G|G" ~ 0,
        SNP=="G|G/G/G" ~ 0,
        SNP=='A|A/C|C' ~ 1,
        SNP=='A/A/C|C' ~ 1,
        SNP=='A|A/C/C' ~ 1,
        SNP=='A/A/C/C' ~ 1,
        SNP=="A|A/G|G" ~ 1,
        SNP=="A/A/G/G" ~ 1,
        SNP=="A/A/G|G" ~ 1,
        SNP=="A|A/G/G" ~ 1,
        SNP=="A|A/T/T" ~ 1,
        SNP=="A|A/T|T" ~ 1,
        SNP=="A/A/T/T" ~ 1,
        SNP=="A/A/T|T" ~ 1,
        SNP=='C|C/A|A' ~ 1,
        SNP=='C/C/A|A' ~ 1,
        SNP=='C|C/A/A' ~ 1,
        SNP=='C/C/A/A' ~ 1,
        SNP=='C|C/G|G' ~ 1,
        SNP=='C/C/G|G' ~ 1,
        SNP=='C|C/G/G' ~ 1,
        SNP=='C/C/G/G' ~ 1,
        SNP=='C|C/T|T' ~ 1,
        SNP=='C/C/T|T' ~ 1,
        SNP=='C|C/T/T' ~ 1,
        SNP=='C/C/T/T' ~ 1,
        SNP=='G|G/A|A' ~ 1,
        SNP=='G/G/A|A' ~ 1,
        SNP=='G|G/A/A' ~ 1,
        SNP=='G/G/A/A' ~ 1,
        SNP=='G|G/C|C' ~ 1,
        SNP=='G/G/C|C' ~ 1,
        SNP=='G|G/C/C' ~ 1,
        SNP=='G/G/C/C' ~ 1,
        SNP=='G|G/T|T' ~ 1,
        SNP=='G/G/T|T' ~ 1,
        SNP=='G|G/T/T' ~ 1,
        SNP=='G/G/T/T' ~ 1,
        SNP=='T|T/A|A' ~ 1,
        SNP=='T/T/A|A' ~ 1,
        SNP=='T|T/A/A' ~ 1,
        SNP=='T/T/A/A' ~ 1,
        SNP=='T|T/C|C' ~ 1,
        SNP=='T/T/C|C' ~ 1,
        SNP=='T|T/C/C' ~ 1,
        SNP=='T/T/C/C' ~ 1,
        SNP=='T|T/G|G' ~ 1,
        SNP=='T/T/G|G' ~ 1,
        SNP=='T|T/G/G' ~ 1,
        SNP=='T/T/G/G' ~ 1))
    
    SNP_9_chrom4 = filter(SNP_8_chrom4, SNP != "A|A/NA")
    SNP_10_chrom4 = filter(SNP_9_chrom4, SNP != "A/A/NA")
    SNP_11_chrom4 = filter(SNP_10_chrom4, SNP != "NA/A/A")
    SNP_12_chrom4 = filter(SNP_11_chrom4, SNP != "NA/A|A")
    
    SNP_13_chrom4 = filter(SNP_12_chrom4, SNP != "C|C/NA")
    SNP_14_chrom4 = filter(SNP_13_chrom4, SNP != "C/C/NA")
    SNP_15_chrom4 = filter(SNP_14_chrom4, SNP != "NA/C/C")
    SNP_16_chrom4 = filter(SNP_15_chrom4, SNP != "NA/C|C")
    
    SNP_17_chrom4 = filter(SNP_16_chrom4, SNP != "T|T/NA")
    SNP_18_chrom4 = filter(SNP_17_chrom4, SNP != "T/T/NA")
    SNP_19_chrom4 = filter(SNP_18_chrom4, SNP != "NA/T/T")
    SNP_20_chrom4 = filter(SNP_19_chrom4, SNP != "NA/T|T")
    
    SNP_21_chrom4 = filter(SNP_20_chrom4, SNP != "G|G/NA")
    SNP_22_chrom4 = filter(SNP_21_chrom4, SNP != "G/G/NA")
    SNP_23_chrom4 = filter(SNP_22_chrom4, SNP != "NA/G/G")
    SNP_24_chrom4 = filter(SNP_23_chrom4, SNP != "NA/G|G")
    
    head(SNP_24_chrom4,200)
    
    SNP_24bis_chrom4 = na.omit(SNP_24_chrom4)
    head(SNP_24_chrom4,200)
    
    #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
    
    SNP_25_chrom4 =  SNP_24bis_chrom4 %>%
      group_by(GENE) %>%
      summarise(score = sum(score))
    #Remettre la position des gènes et leur chromosome
    
    SNP_26_chrom4 = merge(x = SNP_25_chrom4, SNP_GENES_2,all.x=T,all.y=F)
    SNP_27_chrom4 = SNP_26_chrom4 %>% distinct()
    
    #Ajouter une colonne position et taille du gène
    
    SNP_28_chrom4 = SNP_27_chrom4 %>% mutate((SNP_27_chrom4$FIN+SNP_27_chrom4$DEBUT)/2)
    SNP_29_chrom4=  SNP_28_chrom4 %>% mutate((SNP_28_chrom4$FIN-SNP_28_chrom4$DEBUT)+1)
    head(SNP_29_chrom4)
    colnames(SNP_29_chrom4) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
    head(SNP_29_chrom4)
    
    #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
    #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
    SNP_30_chrom4 = SNP_29_chrom4 %>% 
      mutate(SNPs_par_pb=SCORE/TAILLE)
    
    
    ####CHROM5
    SNP_syn_chrom5 = SNPs_syn_Eff1 %>% 
      filter(CHROM=="CBS815.71p_contig_5")
    
    
    SNP_4_chrom5=SNP_3 %>% 
      filter(CHROM == "CBS815.71p_contig_5")
    
    #Ordonner en croissant par positions puis par gène quand 2 fois la même position pour que 
    #quand deux fois le même SNP dans deux gènes : celui du "premier" gène dans le chromosome au niveau de la position seulement soit gardé
    #Permet d'avoir un critère répétable
    
    SNP_5_chrom5=merge(x = SNP_4_chrom5, y = SNP_GENES_1, all.x=TRUE)
    
    SNP_5bis_chrom5 = merge(x = SNP_5_chrom5, y = SNP_syn_chrom5, all.x=FALSE, all.y=FALSE)
    
    SNP_6_chrom5= SNP_5bis_chrom5 %>% arrange(POS,GENE)
    SNP_7_chrom5=SNP_6_chrom5 %>% distinct(POS,.keep_all = TRUE)
    
    #On ajouter une colonne score 
    SNP_8_chrom5=SNP_7_chrom5 %>%
      mutate(score = case_when(
        SNP=="A|A/A|A" ~ 0,
        SNP=="A/A/A/A" ~ 0,
        SNP=="A/A/A|A" ~ 0,
        SNP=="A|A/A/A" ~ 0,
        SNP=="T|T/T|T" ~ 0,
        SNP=="T/T/T/T" ~ 0,
        SNP=="T/T/T|T" ~ 0,
        SNP=="T|T/T/T" ~ 0,
        SNP=="C|C/C|C" ~ 0,
        SNP=="C/C/C/C" ~ 0,
        SNP=="C/C/C|C" ~ 0,
        SNP=="C|C/C/C" ~ 0,
        SNP=="G|G/G|G" ~ 0,
        SNP=="G/G/G/G" ~ 0,
        SNP=="G/G/G|G" ~ 0,
        SNP=="G|G/G/G" ~ 0,
        SNP=='A|A/C|C' ~ 1,
        SNP=='A/A/C|C' ~ 1,
        SNP=='A|A/C/C' ~ 1,
        SNP=='A/A/C/C' ~ 1,
        SNP=="A|A/G|G" ~ 1,
        SNP=="A/A/G/G" ~ 1,
        SNP=="A/A/G|G" ~ 1,
        SNP=="A|A/G/G" ~ 1,
        SNP=="A|A/T/T" ~ 1,
        SNP=="A|A/T|T" ~ 1,
        SNP=="A/A/T/T" ~ 1,
        SNP=="A/A/T|T" ~ 1,
        SNP=='C|C/A|A' ~ 1,
        SNP=='C/C/A|A' ~ 1,
        SNP=='C|C/A/A' ~ 1,
        SNP=='C/C/A/A' ~ 1,
        SNP=='C|C/G|G' ~ 1,
        SNP=='C/C/G|G' ~ 1,
        SNP=='C|C/G/G' ~ 1,
        SNP=='C/C/G/G' ~ 1,
        SNP=='C|C/T|T' ~ 1,
        SNP=='C/C/T|T' ~ 1,
        SNP=='C|C/T/T' ~ 1,
        SNP=='C/C/T/T' ~ 1,
        SNP=='G|G/A|A' ~ 1,
        SNP=='G/G/A|A' ~ 1,
        SNP=='G|G/A/A' ~ 1,
        SNP=='G/G/A/A' ~ 1,
        SNP=='G|G/C|C' ~ 1,
        SNP=='G/G/C|C' ~ 1,
        SNP=='G|G/C/C' ~ 1,
        SNP=='G/G/C/C' ~ 1,
        SNP=='G|G/T|T' ~ 1,
        SNP=='G/G/T|T' ~ 1,
        SNP=='G|G/T/T' ~ 1,
        SNP=='G/G/T/T' ~ 1,
        SNP=='T|T/A|A' ~ 1,
        SNP=='T/T/A|A' ~ 1,
        SNP=='T|T/A/A' ~ 1,
        SNP=='T/T/A/A' ~ 1,
        SNP=='T|T/C|C' ~ 1,
        SNP=='T/T/C|C' ~ 1,
        SNP=='T|T/C/C' ~ 1,
        SNP=='T/T/C/C' ~ 1,
        SNP=='T|T/G|G' ~ 1,
        SNP=='T/T/G|G' ~ 1,
        SNP=='T|T/G/G' ~ 1,
        SNP=='T/T/G/G' ~ 1))
    
    SNP_9_chrom5 = filter(SNP_8_chrom5, SNP != "A|A/NA")
    SNP_10_chrom5 = filter(SNP_9_chrom5, SNP != "A/A/NA")
    SNP_11_chrom5 = filter(SNP_10_chrom5, SNP != "NA/A/A")
    SNP_12_chrom5 = filter(SNP_11_chrom5, SNP != "NA/A|A")
    
    SNP_13_chrom5 = filter(SNP_12_chrom5, SNP != "C|C/NA")
    SNP_14_chrom5 = filter(SNP_13_chrom5, SNP != "C/C/NA")
    SNP_15_chrom5 = filter(SNP_14_chrom5, SNP != "NA/C/C")
    SNP_16_chrom5 = filter(SNP_15_chrom5, SNP != "NA/C|C")
    
    SNP_17_chrom5 = filter(SNP_16_chrom5, SNP != "T|T/NA")
    SNP_18_chrom5 = filter(SNP_17_chrom5, SNP != "T/T/NA")
    SNP_19_chrom5 = filter(SNP_18_chrom5, SNP != "NA/T/T")
    SNP_20_chrom5 = filter(SNP_19_chrom5, SNP != "NA/T|T")
    
    SNP_21_chrom5 = filter(SNP_20_chrom5, SNP != "G|G/NA")
    SNP_22_chrom5 = filter(SNP_21_chrom5, SNP != "G/G/NA")
    SNP_23_chrom5 = filter(SNP_22_chrom5, SNP != "NA/G/G")
    SNP_24_chrom5 = filter(SNP_23_chrom5, SNP != "NA/G|G")
    
    head(SNP_24_chrom5,200)
    
    SNP_24bis_chrom5 = na.omit(SNP_24_chrom5)
    head(SNP_24_chrom5,200)
    
    #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
    
    SNP_25_chrom5 =  SNP_24bis_chrom5 %>%
      group_by(GENE) %>%
      summarise(score = sum(score))
    #Remettre la position des gènes et leur chromosome
    
    SNP_26_chrom5 = merge(x = SNP_25_chrom5, SNP_GENES_2,all.x=T,all.y=F)
    SNP_27_chrom5 = SNP_26_chrom5 %>% distinct()
    
    #Ajouter une colonne position et taille du gène
    
    SNP_28_chrom5 = SNP_27_chrom5 %>% mutate((SNP_27_chrom5$FIN+SNP_27_chrom5$DEBUT)/2)
    SNP_29_chrom5=  SNP_28_chrom5 %>% mutate((SNP_28_chrom5$FIN-SNP_28_chrom5$DEBUT)+1)
    head(SNP_29_chrom5)
    colnames(SNP_29_chrom5) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
    head(SNP_29_chrom5)
    
    #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
    #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
    SNP_30_chrom5 = SNP_29_chrom5 %>% 
      mutate(SNPs_par_pb=SCORE/TAILLE)
    
    ###CHROMOSOME 6
    SNP_syn_chrom6 = SNPs_syn_Eff1 %>% 
      filter(CHROM=="CBS815.71p_contig_6")
    
    SNP_4_chrom6=SNP_3 %>% 
      filter(CHROM == "CBS815.71p_contig_6")
    
    #Ordonner en croissant par positions puis par gène quand 2 fois la même position pour que 
    #quand deux fois le même SNP dans deux gènes : celui du "premier" gène dans le chromosome au niveau de la position seulement soit gardé
    #Permet d'avoir un critère répétable
    SNP_5_chrom6=merge(x = SNP_4_chrom6, y = SNP_GENES_1, all.x=TRUE)
    
    SNP_5bis_chrom6 = merge(x = SNP_5_chrom6, y = SNP_syn_chrom6, all.x=FALSE, all.y=FALSE)
    
    SNP_6_chrom6= SNP_5bis_chrom6 %>% arrange(POS,GENE)
    SNP_7_chrom6=SNP_6_chrom6 %>% distinct(POS,.keep_all = TRUE)
    
    #On ajouter une colonne score 
    SNP_8_chrom6=SNP_7_chrom6 %>%
      mutate(score = case_when(
        SNP=="A|A/A|A" ~ 0,
        SNP=="A/A/A/A" ~ 0,
        SNP=="A/A/A|A" ~ 0,
        SNP=="A|A/A/A" ~ 0,
        SNP=="T|T/T|T" ~ 0,
        SNP=="T/T/T/T" ~ 0,
        SNP=="T/T/T|T" ~ 0,
        SNP=="T|T/T/T" ~ 0,
        SNP=="C|C/C|C" ~ 0,
        SNP=="C/C/C/C" ~ 0,
        SNP=="C/C/C|C" ~ 0,
        SNP=="C|C/C/C" ~ 0,
        SNP=="G|G/G|G" ~ 0,
        SNP=="G/G/G/G" ~ 0,
        SNP=="G/G/G|G" ~ 0,
        SNP=="G|G/G/G" ~ 0,
        SNP=='A|A/C|C' ~ 1,
        SNP=='A/A/C|C' ~ 1,
        SNP=='A|A/C/C' ~ 1,
        SNP=='A/A/C/C' ~ 1,
        SNP=="A|A/G|G" ~ 1,
        SNP=="A/A/G/G" ~ 1,
        SNP=="A/A/G|G" ~ 1,
        SNP=="A|A/G/G" ~ 1,
        SNP=="A|A/T/T" ~ 1,
        SNP=="A|A/T|T" ~ 1,
        SNP=="A/A/T/T" ~ 1,
        SNP=="A/A/T|T" ~ 1,
        SNP=='C|C/A|A' ~ 1,
        SNP=='C/C/A|A' ~ 1,
        SNP=='C|C/A/A' ~ 1,
        SNP=='C/C/A/A' ~ 1,
        SNP=='C|C/G|G' ~ 1,
        SNP=='C/C/G|G' ~ 1,
        SNP=='C|C/G/G' ~ 1,
        SNP=='C/C/G/G' ~ 1,
        SNP=='C|C/T|T' ~ 1,
        SNP=='C/C/T|T' ~ 1,
        SNP=='C|C/T/T' ~ 1,
        SNP=='C/C/T/T' ~ 1,
        SNP=='G|G/A|A' ~ 1,
        SNP=='G/G/A|A' ~ 1,
        SNP=='G|G/A/A' ~ 1,
        SNP=='G/G/A/A' ~ 1,
        SNP=='G|G/C|C' ~ 1,
        SNP=='G/G/C|C' ~ 1,
        SNP=='G|G/C/C' ~ 1,
        SNP=='G/G/C/C' ~ 1,
        SNP=='G|G/T|T' ~ 1,
        SNP=='G/G/T|T' ~ 1,
        SNP=='G|G/T/T' ~ 1,
        SNP=='G/G/T/T' ~ 1,
        SNP=='T|T/A|A' ~ 1,
        SNP=='T/T/A|A' ~ 1,
        SNP=='T|T/A/A' ~ 1,
        SNP=='T/T/A/A' ~ 1,
        SNP=='T|T/C|C' ~ 1,
        SNP=='T/T/C|C' ~ 1,
        SNP=='T|T/C/C' ~ 1,
        SNP=='T/T/C/C' ~ 1,
        SNP=='T|T/G|G' ~ 1,
        SNP=='T/T/G|G' ~ 1,
        SNP=='T|T/G/G' ~ 1,
        SNP=='T/T/G/G' ~ 1))
    
    SNP_9_chrom6 = filter(SNP_8_chrom6, SNP != "A|A/NA")
    SNP_10_chrom6 = filter(SNP_9_chrom6, SNP != "A/A/NA")
    SNP_11_chrom6 = filter(SNP_10_chrom6, SNP != "NA/A/A")
    SNP_12_chrom6 = filter(SNP_11_chrom6, SNP != "NA/A|A")
    
    SNP_13_chrom6 = filter(SNP_12_chrom6, SNP != "C|C/NA")
    SNP_14_chrom6 = filter(SNP_13_chrom6, SNP != "C/C/NA")
    SNP_15_chrom6 = filter(SNP_14_chrom6, SNP != "NA/C/C")
    SNP_16_chrom6 = filter(SNP_15_chrom6, SNP != "NA/C|C")
    
    SNP_17_chrom6 = filter(SNP_16_chrom6, SNP != "T|T/NA")
    SNP_18_chrom6 = filter(SNP_17_chrom6, SNP != "T/T/NA")
    SNP_19_chrom6 = filter(SNP_18_chrom6, SNP != "NA/T/T")
    SNP_20_chrom6 = filter(SNP_19_chrom6, SNP != "NA/T|T")
    
    SNP_21_chrom6 = filter(SNP_20_chrom6, SNP != "G|G/NA")
    SNP_22_chrom6 = filter(SNP_21_chrom6, SNP != "G/G/NA")
    SNP_23_chrom6 = filter(SNP_22_chrom6, SNP != "NA/G/G")
    SNP_24_chrom6 = filter(SNP_23_chrom6, SNP != "NA/G|G")
    
    head(SNP_24_chrom6,200)
    
    SNP_24bis_chrom6 = na.omit(SNP_24_chrom6)
    head(SNP_24_chrom6,200)
    
    #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
    
    SNP_25_chrom6 =  SNP_24bis_chrom6 %>%
      group_by(GENE) %>%
      summarise(score = sum(score))
    #Remettre la position des gènes et leur chromosome
    
    SNP_26_chrom6 = merge(x = SNP_25_chrom6, SNP_GENES_2,all.x=T,all.y=F)
    SNP_27_chrom6 = SNP_26_chrom6 %>% distinct()
    
    #Ajouter une colonne position et taille du gène
    
    SNP_28_chrom6 = SNP_27_chrom6 %>% mutate((SNP_27_chrom6$FIN+SNP_27_chrom6$DEBUT)/2)
    SNP_29_chrom6=  SNP_28_chrom6 %>% mutate((SNP_28_chrom6$FIN-SNP_28_chrom6$DEBUT)+1)
    head(SNP_29_chrom6)
    colnames(SNP_29_chrom6) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
    head(SNP_29_chrom6)
    
    #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
    #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
    SNP_30_chrom6 = SNP_29_chrom6 %>% 
      mutate(SNPs_par_pb=SCORE/TAILLE)
    
    ##CHROMOSOME 7
    SNP_syn_chrom7 = SNPs_syn_Eff1 %>% 
      filter(CHROM=="CBS815.71p_contig_7")
    
    SNP_4_chrom7=SNP_3 %>% 
      filter(CHROM == "CBS815.71p_contig_7")
    
    #Ordonner en croissant par positions puis par gène quand 2 fois la même position pour que 
    #quand deux fois le même SNP dans deux gènes : celui du "premier" gène dans le chromosome au niveau de la position seulement soit gardé
    #Permet d'avoir un critère répétable
    
    SNP_5_chrom7=merge(x = SNP_4_chrom7, y = SNP_GENES_1, all.x=TRUE)
    
    SNP_5bis_chrom7 = merge(x = SNP_5_chrom7, y = SNP_syn_chrom7, all.x=FALSE, all.y=FALSE)
    
    SNP_6_chrom7= SNP_5bis_chrom7 %>% arrange(POS,GENE)
    
    SNP_7_chrom7=SNP_6_chrom7 %>% distinct(POS,.keep_all = TRUE)
    
    #On ajouter une colonne score 
    SNP_8_chrom7=SNP_7_chrom7 %>%
      mutate(score = case_when(
        SNP=="A|A/A|A" ~ 0,
        SNP=="A/A/A/A" ~ 0,
        SNP=="A/A/A|A" ~ 0,
        SNP=="A|A/A/A" ~ 0,
        SNP=="T|T/T|T" ~ 0,
        SNP=="T/T/T/T" ~ 0,
        SNP=="T/T/T|T" ~ 0,
        SNP=="T|T/T/T" ~ 0,
        SNP=="C|C/C|C" ~ 0,
        SNP=="C/C/C/C" ~ 0,
        SNP=="C/C/C|C" ~ 0,
        SNP=="C|C/C/C" ~ 0,
        SNP=="G|G/G|G" ~ 0,
        SNP=="G/G/G/G" ~ 0,
        SNP=="G/G/G|G" ~ 0,
        SNP=="G|G/G/G" ~ 0,
        SNP=='A|A/C|C' ~ 1,
        SNP=='A/A/C|C' ~ 1,
        SNP=='A|A/C/C' ~ 1,
        SNP=='A/A/C/C' ~ 1,
        SNP=="A|A/G|G" ~ 1,
        SNP=="A/A/G/G" ~ 1,
        SNP=="A/A/G|G" ~ 1,
        SNP=="A|A/G/G" ~ 1,
        SNP=="A|A/T/T" ~ 1,
        SNP=="A|A/T|T" ~ 1,
        SNP=="A/A/T/T" ~ 1,
        SNP=="A/A/T|T" ~ 1,
        SNP=='C|C/A|A' ~ 1,
        SNP=='C/C/A|A' ~ 1,
        SNP=='C|C/A/A' ~ 1,
        SNP=='C/C/A/A' ~ 1,
        SNP=='C|C/G|G' ~ 1,
        SNP=='C/C/G|G' ~ 1,
        SNP=='C|C/G/G' ~ 1,
        SNP=='C/C/G/G' ~ 1,
        SNP=='C|C/T|T' ~ 1,
        SNP=='C/C/T|T' ~ 1,
        SNP=='C|C/T/T' ~ 1,
        SNP=='C/C/T/T' ~ 1,
        SNP=='G|G/A|A' ~ 1,
        SNP=='G/G/A|A' ~ 1,
        SNP=='G|G/A/A' ~ 1,
        SNP=='G/G/A/A' ~ 1,
        SNP=='G|G/C|C' ~ 1,
        SNP=='G/G/C|C' ~ 1,
        SNP=='G|G/C/C' ~ 1,
        SNP=='G/G/C/C' ~ 1,
        SNP=='G|G/T|T' ~ 1,
        SNP=='G/G/T|T' ~ 1,
        SNP=='G|G/T/T' ~ 1,
        SNP=='G/G/T/T' ~ 1,
        SNP=='T|T/A|A' ~ 1,
        SNP=='T/T/A|A' ~ 1,
        SNP=='T|T/A/A' ~ 1,
        SNP=='T/T/A/A' ~ 1,
        SNP=='T|T/C|C' ~ 1,
        SNP=='T/T/C|C' ~ 1,
        SNP=='T|T/C/C' ~ 1,
        SNP=='T/T/C/C' ~ 1,
        SNP=='T|T/G|G' ~ 1,
        SNP=='T/T/G|G' ~ 1,
        SNP=='T|T/G/G' ~ 1,
        SNP=='T/T/G/G' ~ 1))
    
    SNP_9_chrom7 = filter(SNP_8_chrom7, SNP != "A|A/NA")
    SNP_10_chrom7 = filter(SNP_9_chrom7, SNP != "A/A/NA")
    SNP_11_chrom7 = filter(SNP_10_chrom7, SNP != "NA/A/A")
    SNP_12_chrom7 = filter(SNP_11_chrom7, SNP != "NA/A|A")
    
    SNP_13_chrom7 = filter(SNP_12_chrom7, SNP != "C|C/NA")
    SNP_14_chrom7 = filter(SNP_13_chrom7, SNP != "C/C/NA")
    SNP_15_chrom7 = filter(SNP_14_chrom7, SNP != "NA/C/C")
    SNP_16_chrom7 = filter(SNP_15_chrom7, SNP != "NA/C|C")
    
    SNP_17_chrom7 = filter(SNP_16_chrom7, SNP != "T|T/NA")
    SNP_18_chrom7 = filter(SNP_17_chrom7, SNP != "T/T/NA")
    SNP_19_chrom7 = filter(SNP_18_chrom7, SNP != "NA/T/T")
    SNP_20_chrom7 = filter(SNP_19_chrom7, SNP != "NA/T|T")
    
    SNP_21_chrom7 = filter(SNP_20_chrom7, SNP != "G|G/NA")
    SNP_22_chrom7 = filter(SNP_21_chrom7, SNP != "G/G/NA")
    SNP_23_chrom7 = filter(SNP_22_chrom7, SNP != "NA/G/G")
    SNP_24_chrom7 = filter(SNP_23_chrom7, SNP != "NA/G|G")
    
    head(SNP_24_chrom7,200)
    
    SNP_24bis_chrom7 = na.omit(SNP_24_chrom7)
    head(SNP_24_chrom7,200)
    
    #Faire un fichier avec une seule ligne par gène, et dans "somme", l'addition des "scores"
    
    SNP_25_chrom7 =  SNP_24bis_chrom7 %>%
      group_by(GENE) %>%
      summarise(score = sum(score))
    #Remettre la position des gènes et leur chromosome
    
    SNP_26_chrom7 = merge(x = SNP_25_chrom7, SNP_GENES_2,all.x=T,all.y=F)
    SNP_27_chrom7 = SNP_26_chrom7 %>% distinct()
    
    #Ajouter une colonne position et taille du gène
    
    SNP_28_chrom7 = SNP_27_chrom7 %>% mutate((SNP_27_chrom7$FIN+SNP_27_chrom7$DEBUT)/2)
    SNP_29_chrom7=  SNP_28_chrom7 %>% mutate((SNP_28_chrom7$FIN-SNP_28_chrom7$DEBUT)+1)
    head(SNP_29_chrom7)
    colnames(SNP_29_chrom7) <- c("GENE","SCORE","CHROM","DEBUT","FIN","POS_GENE","TAILLE")
    head(SNP_29_chrom7)
    
    #On ajouter une colonne dans laquelle on obtient le nombre de SNPs / la taille du gène pour obtenir
    #pour chaque gène le nombre de SNPs moyenné par la taille du gène, donc un nombre entre 0 et 1
    SNP_30_chrom7 = SNP_29_chrom7 %>% 
      mutate(SNPs_par_pb=SCORE/TAILLE)
    
    #On fait un seul tableau en concaténant tous les tableaux
    
    SNP_31 <- bind_rows(SNP_30_chrom1, SNP_30_chrom2, SNP_30_chrom3, SNP_30_chrom4, SNP_30_chrom5,SNP_30_chrom6,SNP_30_chrom7)
    head(SNP_31)
    
    #On représente l'hétérozygotie (Nombre de SNPs hétérozygote par gènes par paire de base)
    #sur l'ensemble de la souche 
    
    plot<-ggplot(SNP_31, aes(y=SNPs_par_pb , x=POS_GENE/1000)) +
      labs(x=" ", y="SNPs synonymes par gène moyennée par pb") +
      geom_point(aes(color=variable),size=0.5,alpha=0.8,colour = "black") + 
      facet_grid( CHROM ~ ., scale="free")+ theme(legend.position="none")
    plot
    
    #On enregistre une table avec uniquement le premier chromosome
    SNP_32 = SNP_31 %>% 
      mutate(SOUCHE=i)
    
    write.table(SNP_32, file=paste0("VERIF_CORRECTION_souche_ENTIERE_SNPsSyn_avecHET",i,".txt"),row.names=T,col.names=T)
    
    #On enregistre une table avec tous les chromosomes
    SNP_32_chrom1 = SNP_30_chrom1 %>% 
      mutate(SOUCHE=i)
    
    write.table(SNP_32_chrom1, file=paste0("VERIF_CORRECTION_souche_SNPsSyn_avecHET",i,".txt"),row.names=T,col.names=T)
    
    
  }
