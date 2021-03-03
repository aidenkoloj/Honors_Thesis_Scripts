library(tidyverse)
library(readr)
library(dplyr)
library(readxl)
library(lsa)
library(reshape2)


# Created by Aiden Kolodziej 

metra_data = read.csv("NIDDK_covametra_data_all_soup_neg_features_includes_BWF_annotations.csv")
ascr3 =  metra_data %>% filter(feature == "ascr#3_metabolite") %>% rename( abundance_ascr3 = abundance, replicate_ascr3 = replicate) %>% select(replicate_ascr3, abundance_ascr3)
ascr3_rep = bind_rows(replicate((nrow(metra_data)/24), ascr3, simplify = FALSE))
metra_data_ascr3 = cbind(metra_data,ascr3_rep)
metra_data_ascr3_norm = metra_data_ascr3 %>% mutate(abundance_ascr3_norm = abundance/abundance_ascr3) %>% select(replicate,feature, abundance_ascr3_norm)

# CovaMeTra -------------------------------------------------------------------------

#Takes 'x' a feature name and returns highest correlated features. 
#Metabolite features must have "_metabolite" on the end 

#-------------------------------------------------------------------------

covametra = function(x){
  
  input_data <- 
    metra_data_ascr3_norm %>% 
    #Select columns of interest
    filter(feature == x) 
  
  df <- data.frame(a = c(input_data$feature), b = c(input_data$abundance_ascr3_norm)) 
  df <- bind_rows(replicate((nrow(metra_data_ascr3_norm)/24), df, simplify = FALSE))
  combined <- cbind(metra_data_ascr3_norm, df) 
  combined_r2 <- combined %>% rename(input_names = a, input_abundance = b) %>% group_by(feature) %>% 
    summarise(r_2 = cor(abundance_ascr3_norm,input_abundance), method = "pearson") 
  combined_r2$r_2 = combined_r2$r_2*combined_r2$r_2
  combined_r2 <- combined_r2 %>% 
    arrange(desc(r_2)) 
  return(combined_r2)
  
  
}
#-------------------------------------------------------------------------
# Plot 2 features x,y graph
# Data and colors
# -------------------------------------------------------------------------
cbPalette <- c(rep("#cf6090",6), rep("#3853a4",6), rep("#78c6cf",6), rep("#faa41a",6))

plot_2features = function(f1, f2){
  
  title = paste(f1, f2, sep ="/")
  
  df1 <- metra_data_ascr3_norm %>% filter(feature == f1) %>% rename( df1_value = abundance_ascr3_norm)
  df2 <- metra_data_ascr3_norm %>% filter(feature == f2) %>% rename( df2_value = abundance_ascr3_norm) %>% select(df2_value)
  df_combined <- cbind(df1, df2) 
  print(df_combined)
  contains_r2_value = lm(df1_value ~ df2_value, data=df_combined)
  
  ggplot(df_combined, aes(df1_value, df2_value))+
    geom_point(aes(color = replicate), size = 2)+
    ggtitle(title)+
    #ggtitle("ascr.3/ascr73frag")+
    theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 5, face = "bold"))+
    theme_bw() + theme(panel.border = element_blank(),
                       text=element_text(size=10,  family="Arial"), legend.position = "none")+
    theme(plot.title = element_text(hjust = 0.5, size =10, vjust = 5, face = "bold"))+
    theme(plot.margin=unit(c(1,1.5,1.5,1.2),"cm"))+
    geom_smooth(method = "lm", col = "gray", se = F)+
    labs(x = as.character(f1), y = as.character(f2))+
    #labs( x = "ascr#3", y = "ascr73frag")+
    geom_text(x = -Inf, y = Inf,hjust = 0, vjust = 1, label = paste("R2 = ",signif(summary(contains_r2_value)$r.squared, 5)))+
    scale_colour_manual(values=cbPalette)
}


# Plot single feature -------------------------------------------------------------------------

df1 = dplyr::filter(metra_data, feature == "cest-3") %>% mutate(strain = c(rep("BRC20067",6), rep("CB4856",6), rep("DL238", 6), rep("N2", 6)))
plot1feature <- function(x){
  ggplot(x, aes(strain,abundance))+
    labs(y= "Normalized Counts", x = "Strain")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="black", width=0.2)+
    geom_point(aes(color = replicate), size = 2)+
    theme_bw() + theme(panel.border = element_blank(), axis.text.x=element_blank())+
    theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 5, face = "bold"))+
    theme(plot.margin=unit(c(1,1.5,1.5,1.2),"cm"))+
    theme(plot.title = element_text(hjust = 0.5, size =20, vjust = 5, face = "bold"))+
    ggtitle("indole glucoside")+
    scale_colour_manual(values=cbPalette)
}

