#Strata detection based on change-point method#

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(magrittr)
library(mcp)
library(rjags)

strains=c("S1.PSN779.h","S1.CBS81571","S2.CBS26270.h","S5.PSQ32","S6.PSN1067.h","S8.PSN644","S9.PSN707","S12.PSN663","S13.PSN1042","S3.PSN1070","S16.PSN1310.h","S10.PSN684.h","S15.PSN377","S4.PSN936.h")
path <- "$path_to_strains"

for (i in strains) {
  
  d2=read.table(paste0(path,"/VERIF_CORRECTION_souche_SNPsSyn_avecHET",i,".txt"), header=T, as.is=T)
  NRR <- read.csv(paste0(path,"/SIZE_NRR.txt"),sep="")
  NRR_strain <- subset(NRR, SOUCHE == as.character(i))
  
  start_NRR <- as.numeric(as.character(NRR_strain$MIN))
  end_NRR <-  as.numeric(as.character(NRR_strain$MAX))
  
  #Selecting columns for gene position and for heterozygous SNPs per bp
  df_1= d2 %>% select("POS_GENE","SNPs_par_pb")
  
  #Selecting the non-recombining region
  df_2 = df_1 %>% filter(POS_GENE>start_NRR-1)
  #disparition virgule des position genes à investiguer
  df_3 = df_2 %>% filter(POS_GENE<=end_NRR)
  
  #Computing models for 2, 3 and 4 strata
  model2 = list(SNPs_par_pb ~ 1, 1~ 1) # 2 strata, 1 change-point
  fit_mcp2 = mcp(model2, data = df_3, par_x = "POS_GENE",chains=4)
  summary(fit_mcp2)
  fit_mcp=fit_mcp2
  
  model3 = list(SNPs_par_pb ~ 1, 1~ 1, 1~ 1) # 3 strata, 2 change-points
  fit_mcp3 = mcp(model3, data = df_3, par_x = "POS_GENE",chains=4)
  summary(fit_mcp3)
  fit_mcp=fit_mcp3
  
  model4 = list(SNPs_par_pb ~ 1, 1~ 1, 1~ 1, 1~ 1) # 4 strata, 3 change-points
  fit_mcp4 = mcp(model4, data = df_3, par_x = "POS_GENE",chains=4)
  summary(fit_mcp4)
  fit_mcp=fit_mcp4
  
  #opening the strata limits detected with the Vittorelli et al method
  limits <- read.csv(paste0(path,"/boundary_strata_strain_",i,".csv"), header=TRUE)
  
  #plotting each model, indicating the limits of the strata detected with the Vittorelli et al method
  MCP=plot(fit_mcp2, q_fit=TRUE) + ylab(expression(italic(SNPs_par_pb))) + theme_bw() +
    xlab('')+
    ylim(-1,0.0075)+
    geom_point(size = 1, shape = 3)  +
    geom_vline(xintercept=limits$limit1, color="orange") + 
    geom_vline(xintercept=limits$limit2, color="orange")+ 
    geom_vline(xintercept=limits$limit3, color="orange") 
  MCP
  sum2=summary(fit_mcp2)
  
  write.table(sum2,paste0(path,"/summary_mcp2_",i,".txt"),sep=";",row.names = FALSE)
  ggsave(paste0(path,"/mcp2_",i,"_limits_Vit.png"),width =12, height =5, scale=2, dpi=200)
  
  MCP=plot(fit_mcp3, q_fit=TRUE) + ylab(expression(italic(SNPs_par_pb))) + theme_bw() +
    xlab('')+
    ylim(-1,0.0075)+
    geom_point(size = 1, shape = 3)  +
    geom_vline(xintercept=limits$limit1, color="orange") + 
    geom_vline(xintercept=limits$limit2, color="orange")+ 
    geom_vline(xintercept=limits$limit3, color="orange") 
  MCP
  sum3=summary(fit_mcp3)
  
  write.table(sum3,paste0(path,"/summary_mcp3_",i,".txt"),sep=";",row.names = FALSE)
  ggsave(paste0(path,"/mcp3_",i,"_limits_Vit.png"),width =12, height =5, scale=2, dpi=200)
  
  
  MCP=plot(fit_mcp4, q_fit=TRUE) + ylab(expression(italic(SNPs_par_pb))) + theme_bw() +
    xlab('')+
    ylim(-1,0.0075)+
    geom_point(size = 1, shape = 3)  +
    geom_vline(xintercept=limits$limit1, color="orange") + 
    geom_vline(xintercept=limits$limit2, color="orange")+ 
    geom_vline(xintercept=limits$limit3, color="orange") 
  MCP
  sum4=summary(fit_mcp4)
  
  write.table(sum4,paste0(path,"/summary_mcp4_",i,".txt"),sep=";",row.names = FALSE)
  ggsave(paste0(path,"/mcp4_",i,"_limits_Vit.png"),width =12, height =5, scale=2, dpi=200)
  
}
