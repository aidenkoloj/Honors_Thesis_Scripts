library(tidyverse)
library(readr)
library(dplyr)
library(readxl)
library(lsa)

## Created by Aiden Kolodziej for comparison of RNA transcripts by cosine similarity 

### Load in data

data= read_csv("RNA_seq_NIDDK.csv")

data = data %>% rename(gene_name = ext_gene, counts = normalized_counts, meth = sync_meth)

group_means <- data %>%
  select(strain, gene_name, counts)%>%
  group_by(gene_name,strain)%>%
  summarise(counts = mean(counts))%>%
  ungroup(gene_name,strain)
  #mutate(another_column= 6)

gene_names = 
  data %>%
  select(gene_name)%>%
  distinct()

gene_names = gene_names$gene_name
  

### takes some vector as input, returns list of hits based on a set cosine score.

same_pattern = function(pattern){
  hits = c()
  z = nrow(group_means)/4
  q = 0
  w = 1
  
  for (i in c(1:z)){
    q = q+4
    
    b = group_means$counts[w:q]
    c = cosine(pattern, b)
    
    w = w+4
    
    if (c  > 0.998){
      hits <- c(hits,gene_names[i],c)
    }
  }
  
  hits = tibble(hits)
  return(hits)
}



# Use this function to create a vector for a gene to compare expression to other genes

vec_interest = function(a_gene_name){
  vec_data <- data %>%
    select(strain, gene_name, counts)%>%
    group_by(gene_name,strain)%>%
    summarise(counts = mean(counts))%>%
    ungroup(gene_name,strain)%>%
    filter(gene_name == a_gene_name)%>%
    select(counts)%>%
    pull(counts)
  
  return(vec_data)
  
}



# Aquiring vector for compound of interest --------------------------------
### for ascr3-glucose


ascr3glu = read_excel("ascr_glu2.xlsx")

###
#vec_interest = read_csv("iglas.csv")
#vec_interest <- subset(vec_interest, select = -c(iglos1))

first = colnames(ascr3glu)[1]
second = colnames(ascr3glu)[2]

ascr3glu_vec = ascr3glu %>% 
  rename("groups" = first) %>%
  rename("expression" = second)%>%
  select(groups,expression)%>%
  group_by(groups) %>%
  dplyr::summarize(mean = mean(expression, na.rm=TRUE))%>%
  pull(mean)




# Running Program ---------------------------------------------------------


#vec = vec_interest("ugt-53")

a_pattern = c(5201.334, 5142.230, 2284.991, 5029.740)

hits = same_pattern(a_pattern)
hits






# Graphing Function ----------------------------------------------------------------



plotty <- 
  data %>% 
  #Select columns of interest
  select(strain,gene_name,counts, meth) %>% 
  #Filter data according to ugt_name
  #filter(gene_name == hits$hits[2])
  filter(gene_name == "molo-1") 

ggplot(plotty, aes(strain,counts))+
  
  geom_bar(stat = "summary", width =.7, fill = c("#cf6090", "#3853a4", "#78c6cf","#faa41a"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(aes(shape =meth, color = meth))+
  ggtitle(plotty$gene_name[1])+
  theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 1, face = "bold"))







  

