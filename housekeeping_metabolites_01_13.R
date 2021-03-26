library(tidyverse)
library(readr)
library(dplyr)
library(readxl)
library(lsa)
library(reshape2)
library(stringr)

## takes a csv with a mz_rt_comment column and n number of abundances columns
df <- read.csv("NIDDK_Soup_Neg_npeaks5_rt180-1400_FPQ95_v01_tidy.csv")
df_abundances <- df[,3:ncol(df)]
names <- colnames(df_abundances)

mz_rt_comments <- str_split_fixed(df$mz_rt_comment, "_", 3)
df_mz_rt_comments <- cbind(mz_rt_comments, df_abundances)
colnames(df_mz_rt_comments)<- c("mz", "rt", "comment", names)

df_comments_only <- df_mz_rt_comments[!(is.na(df_mz_rt_comments$comment) | df_mz_rt_comments$comment==""), ]
df_comments_only_abundances <- df_comments_only[,4:ncol(df_comments_only)]


compiled_m_value_2 <- tibble(mz = c(), rt = c(), feature = c(), m = numeric())


# -------------------------------------------------------------------------

#for (i in 1:nrow(df_comments_only_abundances))

for (i in 1:nrow(df_comments_only_abundances)){
  #print(df_comments_only_abundances[1,1:4])
  xx <- mapply('/', df_abundances, as.numeric(df_comments_only_abundances[i,1:ncol(df_comments_only_abundances)]))
  #print(xx[i,1:4])
  xx <- as.tibble(xx) %>% mutate_all(~log2(.))
  #print(xx[i,1:4])
  xx<- xx %>% rowwise() %>% mutate(sd = sd(c_across(names[1]:names[ncol(df_abundances)])))
  #print(xx$sd[1])
  m_compute <- (sum(xx$sd)/(nrow(df_abundances)-1))
  print(m_compute)
  compiled_m_value_2 <- compiled_m_value_2 %>% add_row(
    mz = as.character(df_comments_only[i,1]),
    rt = as.character(df_comments_only[i,2]),feature = as.character(df_comments_only[i,3]), m = m_compute)
  print(i)
}


# -------------------------------------------------------------------------



compiled_m_value_2 <- compiled_m_value_2 %>% arrange(m)

write.csv(compiled_m_value_2,
          "/Users/aidenkoloj/Desktop/R_Projects👾/R_Schroeder_Lab/R_Genetic_Background_4/housekeeping_metabolites/m_values_log2trans_NIDDK_Soup_Neg_npeaks5_rt180-1400_FPQ95_v01.csv")


# -------------------------------------------------------------------------
# 213.1126649_597.686_C11H17O4- present hb101
# 301.1656548_520.374_ascr#3

df_tidy = df %>% filter(grepl("ascr", mz_rt_comments)) 
 
norm = df %>% filter( mz_rt_comments == "359.2439627_718.135_ascr#22")
norm = norm[,2:ncol(norm)]

met = df %>% filter( mz_rt_comments == "446.2184264_736.657_icas#10")
met = met[,2:ncol(met)]


norm_melt = melt(norm)
met_melt = melt(met)

met_tidy = met_melt %>% 
  mutate(strains = c(rep("BRC20067", 6), rep("CB4856", 6), rep("DL238",6), rep("N2",6)))
colnames(met_tidy)<-c("feature", "strain", "abundance", "strains")


plotie <- met_melt %>%
  mutate(norm_data = value/norm_melt$value) %>%
  select(mz_rt_comments, variable, norm_data) %>%
  mutate(strains = c(rep("BRC20067", 6), rep("CB4856", 6), rep("DL238",6), rep("N2",6)))

colnames(plotie) <- c("feature", "strain", "abundance", "strains")



graph_gene = function(df){
  ggplot(df, aes(strains,abundance))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    geom_bar(stat = "summary",fill = c("#cf6090", "#3853a4", "#78c6cf","#faa41a"), width =.4)+
    #geom_bar(stat = "summary",fill = c( "#78c6cf","#faa41a"), width =.4)+
    geom_point(color="black", shape=18, size =3)+
    ggtitle(paste0(df[1,1]))+
    theme(plot.title = element_text(hjust = 0.5, size =10, vjust = 5, face = "bold"))+
    labs(y= "Abundance", x = "Strain")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(axis.title.y = element_text(size =18, vjust =2))+
    theme(axis.title.x = element_text(size =20))+
    theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}

graph_gene(plotie)
graph_gene(met_tidy)


# Example of normalization method step by step ----------------------------

df = data.frame( gene = paste0(rep("gene", 15000), rep("_", 15000), rep(1:15000)), 
                 Sample_A1 = c(floor(runif(15000, min=1, max=101))), 
                 Sample_B1 = c(floor(runif(15000, min=1, max=101))),
                 Sample_C1 = c(floor(runif(15000, min=1, max=101))), 
                 Sample_D1 = c(floor(runif(15000, min=1, max=101))),
                 Sample_A2 = c(floor(runif(15000, min=1, max=101))), 
                 Sample_B2 = c(floor(runif(15000, min=1, max=101))),
                 Sample_C2 = c(floor(runif(15000, min=1, max=101))), 
                 Sample_D2 = c(floor(runif(15000, min=1, max=101))),
                 Sample_A3 = c(floor(runif(15000, min=1, max=101))), 
                 Sample_B3 = c(floor(runif(15000, min=1, max=101))),
                 Sample_C3 = c(floor(runif(15000, min=1, max=101))), 
                 Sample_D3 = c(floor(runif(15000, min=1, max=101))))

df_abundances <- df[,2:ncol(df)]

gene_names <- c(as.character(df$gene))

compiled_m_value_x<- tibble(gene = c(), m = c())

for (i in 1:nrow(df_abundances)){
 
dfx <- mapply('/', df_abundances, as.numeric(df_abundances[i,1:ncol(df_abundances)]))

dfx <- as.tibble(dfx) %>% mutate_all(~log2(.))

dfx<- dfx %>% rowwise() %>% mutate(sd = sd(c_across(Sample_A1:Sample_D3)))

m_compute <- (sum(dfx$sd)/(nrow(dfx)-1))

compiled_m_value_x <- compiled_m_value_x %>% 
  add_row(gene = as.character(gene_names[i]), m = m_compute)

print(i)
}


# END ---------------------------------------------------------------------




graph_gene_x(df,"gene_535")

graph_gene_x_ratio(df, "gene_1210","gene_6205")








# x is a dataframe, y is a character of a gene name 

graph_gene_x = function(x, y){
  
  #data tidying 
  x <- melt(x, value.name = "abundance", variable.name = "replicate")
  x1 <- x %>% filter(gene == y) 
  df <- x1 %>% mutate(sample = x1$replicate)
  df$sample <- substr(as.character(x1$replicate), start = 1,stop = nchar(as.character(x1$replicate))-1)
  
  ggplot(df, aes(sample,abundance))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    geom_bar(stat = "summary", fill = c("#cf6090", "#3853a4", "#78c6cf","#faa41a"), width =.4)+
    geom_point(color="black", shape=18, size =3)+
    ggtitle(paste0(y))+
    theme(plot.title = element_text(hjust = 0.5, size =10, vjust = 5, face = "bold"))+
    labs(x= paste0(colnames(df)[2]), y = paste0(colnames(df)[3]))+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(axis.title.y = element_text(size =18, vjust =2))+
    theme(axis.title.x = element_text(size =20))+
    theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}


# x is a dataframe, y is a gene name, z is the gene name you are normalizing to

graph_gene_x_ratio = function(x, y, z){
  
  #data tidying 
  x <- melt(x, value.name = "abundance", variable.name = "replicate")
  x1 <- x %>% filter(gene == y) 
  df <- x1 %>% mutate(sample = x1$replicate)
  df$sample <- substr(as.character(x1$replicate), start = 1,stop = nchar(as.character(x1$replicate))-1)
  
  x2 <- x %>% filter(gene == z)
  df2 <- x2 %>% mutate(sample = x2$replicate)
  df2$sample <- substr(as.character(x2$replicate), start = 1,stop = nchar(as.character(x2$replicate))-1)
  print(df2)
  
  df <- df %>% mutate(abundance_ratio = df$abundance/df2$abundance)

  ggplot(df, aes(sample,abundance_ratio))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    geom_bar(stat = "summary", fill = c("#cf6090", "#3853a4", "#78c6cf","#faa41a"), width =.4)+
    geom_point(color="black", shape=18, size =3)+
    ggtitle(paste0(y,"/",z))+
    theme(plot.title = element_text(hjust = 0.5, size =10, vjust = 5, face = "bold"))+
    labs(x= paste0(colnames(df)[2]), y = paste0(colnames(df)[3]))+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(axis.title.y = element_text(size =18, vjust =2))+
    theme(axis.title.x = element_text(size =20))+
    theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}

