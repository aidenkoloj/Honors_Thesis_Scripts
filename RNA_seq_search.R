library(tidyverse)
library(readr)
library(dplyr)
library(readxl)

## Created by Aiden Kolodziej for querying and graphing of RNA-seq data. 

### Load in data
data= read_csv("RNA_seq.csv")
data = data %>% rename(method = sync_meth)


## Function: collect data for a given gene
## input should be a string corresponding to a gene name

gene = function(name_of_gene){
  data %>% 
    #Select columns of interest
    select(strain,ext_gene,normalized_counts,method) %>% 
    #Filter data according to ugt_name
    filter(ext_gene == name_of_gene)
}


## Function: graph gene of interest function
## input should be a string corresponding to a gene name

graph_gene = function(gene_name){
  df = gene(gene_name)
  title = df$ext_gene[1]
  ggplot(df, aes_string('strain','normalized_counts'))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    geom_bar(stat = "summary",fill = c("#cf6090", "#3853a4", "#78c6cf","#faa41a"), width =.7)+
    geom_point(aes(color = method, shape = method), size =2)+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 5, face = "bold"))+
    labs(y= "Normalized Counts", x = "Strain")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(axis.title.y = element_text(size =18, vjust =2))+
    theme(axis.title.x = element_text(size =20))+
    theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}



### RUN LINE TO GRAPH GENE OF INTEREST ###

graph_gene("hach-1")

