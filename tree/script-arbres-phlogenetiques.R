Script arbre : 
library(ggtree)
library(tidyverse)
library(reshape2)
genes = c("1625","1774","1497","1495")for (i in genes ) {  mytree <- (paste0("Gene",i,"/rooted/Arbre_TransPol_39g_rooted_Gene",i,".contree"))  tree<-read.tree(mytree)  plot(tree)    #Arbre avec uniquement bootstraps + noms des souches, branches significantes  
ggtree(tree,layout ='rectangular')+    theme(legend.position='none') +    scale_size_discrete(range = c(3,0)) +    scale_shape_manual(values=c(23,27))+    geom_text2(aes(label=label, subset=isTip), hjust=-.1)+    geom_label2(aes(subset =! isTip ,label=label,label.padding = 2,parse=T),size=3)+    geom_treescale()+    ggtitle(paste0("Gene",i,"ROOTED"))    ggsave(paste0("Gene",i,"/rooted/VCF_final_Transpol_Gene",i,"_brut_phylogramme_rooted_legende.pdf"),width =15, height =15, scale=1, dpi=200)      }

